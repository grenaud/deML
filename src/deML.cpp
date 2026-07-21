// vim:ts=8
#include <iostream>
#include <string>
#include <cstring>
#include <sys/stat.h>
#include <sstream>
#include <map>
#include <queue>
#include <vector>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <chrono>
#include <atomic>
#include <limits>
#include <cstdlib>
#include <gzstream.h>


#include <htslib/sam.h>
#include <htslib/hts.h>

#include "RGAssign.h"
#include "FastQParser.h"

#include "JSON.h"
#include "libgab.h"

using namespace std;

typedef struct { 

    ogzstream pairr1;
    ogzstream pairr2;
    ogzstream pairi1;
    ogzstream pairi2;

    ogzstream pairr1f;
    ogzstream pairr2f;
    ogzstream pairi1f;
    ogzstream pairi2f;


} fqwriters;

/****************************************/
/*                                      */
/*          STATIC VARIABLES            */
/*                                      */
/****************************************/


struct tallyForRG{    
    unsigned int assigned;
    unsigned int conflict;
    unsigned int unknown;
    unsigned int wrong;    
};

int qualOffset=33;

static string tagIndex1Seq  = "XI";
static string tagIndex1Qual = "YI";
static string tagIndex2Seq  = "XJ";
static string tagIndex2Qual = "YJ";

static double rgScoreCutoff  = 80 ;             // two HQ bases can mismatch
static double fracConflict   = 20 ;             // top shall be 100x more likely
static double wrongness      = 60 ;             // flag wrong if the wrong pair is a lot more likely
static int    mismatchesTrie = 2;
static int    maxErrorHits   = 20;
PrefixTree<string> * trieKnownString;
static string dashes = "--------------------------------";


/****************************************/
/*                                      */
/*          MULTI-THREADING             */
/*                                      */
/****************************************/

// Number of worker threads. 1 (the default) means the program runs
// single-threaded on the main thread, byte-for-byte identical to the
// original serial implementation. A value >1 activates the ordered
// pipeline below.
static int numberOfThreads = 1;

// Number of input records ("work units") bundled into a single chunk
// that is handed to a worker thread.
static const size_t chunkSizeUnits = 20000;

// Abort on an unrecoverable error from code that may run on a worker thread.
// std::exit() would run atexit handlers and destroy namespace-scope objects
// (e.g. the shared index trie) while other workers still read them, risking a
// crash during shutdown; when threads are live, terminate immediately without
// that teardown instead. The caller must already have printed a flushed
// message (cerr<<...<<endl). Single-threaded behaviour is unchanged.
[[noreturn]] static void fatalExit(){
    if(numberOfThreads > 1)
	std::_Exit(1);
    exit(1);
}

// Protects all shared statistics: the global namesMap tally and the
// per-run unknownSeq/wrongSeq/conflictSeq maps. Only ever contended
// when numberOfThreads>1; uncontended (hence negligible) otherwise.
static std::mutex statsMutex;


// A chunk of work items carrying its position (rank) in the input so
// that output can be re-ordered to exactly match input order.
template<typename T>
struct PipelineChunk {
    long          rank;
    std::vector<T> items;
};

template<typename T>
struct PipelineChunkGreater {
    // top() of the priority_queue is the SMALLEST rank (next to write)
    bool operator()(PipelineChunk<T> * a, PipelineChunk<T> * b) const {
        return a->rank > b->rank;
    }
};


// Ordered parallel map: reads items with produce() (single threaded),
// transforms each item with work() (concurrently on numThreads
// workers), and emits each item with consume() (single threaded, in
// the exact order produce() returned them).
//
//   produce(T&)  -> bool   fill the next item, return false at EOF
//   work(T&)               heavy, thread-safe transformation, in place
//   consume(T&)            write the item (called in input order)
//
// Ordering guarantee: every chunk keeps the rank it was read at; the
// writer only ever emits the chunk whose rank == lastWritten+1, so the
// output order is independent of which worker finished first.
template<typename T, typename ProduceF, typename WorkF, typename ConsumeF>
void runOrderedPipeline(int numThreads, size_t chunkSize,
                        ProduceF produce, WorkF work, ConsumeF consume){

    // ---- serial fast path: preserves the original behaviour exactly ----
    if(numThreads <= 1){
        T item;
        while( produce(item) ){
            work(item);
            consume(item);
        }
        return;
    }

    typedef PipelineChunk<T> Chunk;

    std::queue<Chunk *>                                                         procQ;   // input, FIFO
    std::priority_queue<Chunk *, std::vector<Chunk *>, PipelineChunkGreater<T> > writeQ;  // output, ordered by rank

    std::mutex              procM;
    std::mutex              writeM;
    std::mutex              readerM;
    std::condition_variable procCV;    // procQ gained work OR readDone
    std::condition_variable writeCV;   // writeQ gained a chunk OR the run finished
    std::condition_variable readerCV;  // the writer advanced lastWritten (memory freed)

    std::atomic<bool> readDone(false);
    std::atomic<long> lastWritten(-1);
    std::atomic<long> finalRank( std::numeric_limits<long>::max() ); // real value set when reader finishes
    // Bound how many chunks the reader may run ahead of the writer, which caps
    // total memory (procQ + in-flight + writeQ) without ever blocking a worker.
    const long        maxOutstanding = std::max(2, 3*numThreads);

    // ---- workers: procQ -> work() -> writeQ (heavy assignment, concurrent) ----
    auto workerFn = [&](){
        while(true){
            Chunk * c = NULL;
            {
                std::unique_lock<std::mutex> lk(procM);
                procCV.wait(lk, [&]{ return !procQ.empty() || readDone.load(); });
                if( procQ.empty() )   // readDone && nothing left
                    return;
                c = procQ.front();
                procQ.pop();
            }
            for(size_t i=0;i<c->items.size();i++)
                work(c->items[i]);

            {
                std::lock_guard<std::mutex> lk(writeM);
                writeQ.push(c);       // write queue is unbounded: workers never block
            }
            writeCV.notify_one();      // wake the writer thread
        }
    };

    // ---- writer: drain writeQ strictly in rank order, on its own thread so
    //      output compression overlaps the reader's input decompression ----
    auto writerFn = [&](){
        while(true){
            Chunk * c = NULL;
            {
                std::unique_lock<std::mutex> lk(writeM);
                writeCV.wait_for(lk, std::chrono::milliseconds(50), [&]{
                    return ( !writeQ.empty() && writeQ.top()->rank == lastWritten.load()+1 )
                        || ( lastWritten.load() == finalRank.load() );
                });
                if( !( !writeQ.empty() && writeQ.top()->rank == lastWritten.load()+1 ) ){
                    if( lastWritten.load() == finalRank.load() )
                        return;       // everything has been written
                    continue;         // next chunk not ready yet
                }
                c = writeQ.top();
                writeQ.pop();
            }
            for(size_t i=0;i<c->items.size();i++)
                consume(c->items[i]);
            lastWritten.store(c->rank);
            delete c;
            readerCV.notify_one();     // reader may resume producing
        }
    };

    std::vector<std::thread> workers;
    for(int i=0;i<numThreads;i++)
        workers.push_back( std::thread(workerFn) );
    std::thread writer(writerFn);

    // Hand a finished chunk to the workers, throttling the reader so it never
    // runs more than maxOutstanding chunks ahead of the writer.
    auto pushChunk = [&](Chunk * c){
        while( (c->rank - lastWritten.load()) > maxOutstanding ){
            std::unique_lock<std::mutex> lk(readerM);
            readerCV.wait_for(lk, std::chrono::milliseconds(2),
                              [&]{ return (c->rank - lastWritten.load()) <= maxOutstanding; });
        }
        {
            std::lock_guard<std::mutex> lk(procM);
            procQ.push(c);
        }
        procCV.notify_one();
    };

    // ---- reader (this thread): produce chunks ----
    long    rank = 0;
    Chunk * cur  = new Chunk();
    cur->rank = rank;

    T item;
    while( produce(item) ){
        cur->items.push_back( std::move(item) );
        item = T();
        if( cur->items.size() >= chunkSize ){
            pushChunk(cur);
            rank++;
            cur = new Chunk();
            cur->rank = rank;
        }
    }
    pushChunk(cur);                    // last chunk (may be empty; carries final rank)

    finalRank.store(rank);
    {
        std::lock_guard<std::mutex> lk(procM);   // publish readDone to workers
        readDone.store(true);
    }
    procCV.notify_all();
    { std::lock_guard<std::mutex> lk(writeM); }  // publish finalRank to the writer
    writeCV.notify_one();

    for(size_t i=0;i<workers.size();i++)
        workers[i].join();
    writer.join();
}


map<string,string> p7id2seq;
map<string,string> p5id2seq;


// XXX really, a global?!
map<string,tallyForRG> namesMap; //a map name of RG to count of how many observed

// string getCWD(){
//    char temp[1000];
//    return ( getcwd(temp, 1000) ? string( temp ) : string("") );
// }

struct compareNameRG {
    bool operator() (pair<string,int> i,pair<string,int> j) {
        return (i.second>j.second);
    }
};

struct compareNameTally {
    bool operator() (pair<string,tallyForRG> i,pair<string,tallyForRG> j) { 
        return ( (i.second.assigned+i.second.unknown+i.second.conflict+i.second.wrong)
		 >
		 (j.second.assigned+j.second.unknown+j.second.conflict+j.second.wrong) );
    }
};


// Fetch a Z/H string tag into 'out'; returns false if the tag is absent
// (or is not a string tag).
static bool getZTag( bam1_t *b, const string & tag, string & out ){
    uint8_t *p = bam_aux_get(b, tag.c_str());
    if(!p) return false;
    if( *p != 'Z' && *p != 'H' ) return false;
    char *s = bam_aux2Z(p);
    if(!s) return false;
    out = s;
    return true;
}

// Delete a tag if it is present.
static inline void delTag( bam1_t *b, const char *tag ){
    uint8_t *p = bam_aux_get(b, tag);
    if(p) bam_aux_del(b, p);
}

static string get_string_field( bam1_t *b, const char* name )
{
    uint8_t *p = bam_aux_get(b, name);
    if(!p) return string();
    char ttype = *p;
    if( ttype=='Z' || ttype=='H' ){
        char *s = bam_aux2Z(p);
        return s ? string(s) : string();
    } else if( ttype=='A' ){
        return string(1, (char)bam_aux2A(p));
    }
    return string() ;
}

static bool getIndices( bam1_t *al,string & index1,string & index1Q,string & index2,string & index2Q,bool suppressErrors=false){
    if(!getZTag(al,tagIndex1Seq,index1) ){
        if (!suppressErrors) {
	    cerr << "Cannot retrieve "<<tagIndex1Seq<<" field  "<<bam_get_qname(al) << endl;
	    fatalExit();
        } else
            return false;
    }

    if(!getZTag(al,tagIndex1Qual,index1Q)){
        if (!suppressErrors) {
	    cerr << "Cannot retrieve "<<tagIndex1Qual<<" field  "<<bam_get_qname(al) << endl;
	    fatalExit();
        } else
            return false;
    }

    if(!getZTag(al,tagIndex2Seq,index2)) {
        index2 ="";
        index2Q="";
    }
    else if(!getZTag(al,tagIndex2Qual,index2Q)) { //we got XJ, double indexed
        if (!suppressErrors) {
    	    cerr << "Cannot retrieve "<<tagIndex2Qual<<" field  "<<bam_get_qname(al) << endl;
    	    fatalExit();
        } else
            return false;
    }

    return true;
}


void updateRecord( bam1_t *al, const rgAssignment &rg , bool failBAM)
{
    // get old ZQ field, remove "ICW"
    string zq = get_string_field(al, "ZQ");
    string::iterator p = zq.begin(), q = zq.begin(), e = zq.end() ;
    while( p != e ) {
        if( *p != 'I' && *p != 'C' && *p != 'W' ) 
        {
            *q = *p ;
            ++q ;
        }
        ++p ;
    }
    zq.erase( q, e ) ;

    string predictedGroup=rg.predictedGroup;
    //will overwrite the RG
    bool assigned=true;
    if( rg.predictedGroup.empty() ) {
	predictedGroup="unknown";
    }
    if( rg.conflict ){ zq += 'C' ; assigned=false; }
    if( rg.wrong    ){ zq += 'W' ; assigned=false; }
    if( rg.unknown  ){ zq += 'I' ; assigned=false; }
    {
	// namesMap is shared across worker threads (see runOrderedPipeline)
	std::lock_guard<std::mutex> lk(statsMutex);
	bool incrInTally=false;
	if( rg.conflict ){ if(!incrInTally){ namesMap[predictedGroup].conflict++; incrInTally=true;}  }
	if( rg.wrong    ){ if(!incrInTally){ namesMap[predictedGroup].wrong++;    incrInTally=true;}  }
	if( rg.unknown  ){ if(!incrInTally){ namesMap[predictedGroup].unknown++;  incrInTally=true;}  }
	if(assigned){
	    namesMap[predictedGroup].assigned++;
	}
    }


    if( rg.predictedGroup.empty() ) {
        delTag(al,"RG");
        delTag(al,"Z0");
        delTag(al,"Z1");
        delTag(al,"Z2");
        zq += 'I';
    } else    {
        bam_aux_update_str(al,"RG",rg.predictedGroup.size()+1,rg.predictedGroup.c_str());

        bam_aux_update_int(al,"Z0",(int)round(-10 * rg.logLikelihoodScore));

        if( rg.logRatioTopToSecond <= 0 )
            bam_aux_update_int(al,"Z1",(int)round(-10 * rg.logRatioTopToSecond));
        else
            delTag(al,"Z1") ;



        if( (rg.topWrongToTopCorrect) >= 0)  //if the log odds ratio is greater than 0, the ratio is greater than 1 and the mispairing scenario is more likely than the correct pairing one
            bam_aux_update_int(al,"Z2",(int)round(10.0 * rg.topWrongToTopCorrect));
        else
            delTag(al,"Z2") ;

    }

    // store new ZQ field and set FailedQC flag if it isn't empty if failBAM
    if(failBAM){
	if( !zq.empty() ) al->core.flag |=  BAM_FQCFAIL;
	else              al->core.flag &= ~BAM_FQCFAIL;

	if( zq.empty() )
	    delTag(al,"ZQ") ;
	else
	    bam_aux_update_str(al,"ZQ",zq.size()+1,zq.c_str()) ;
    }
}

inline bool containsNoNs(const string & sN){
    return (sN.find("N") == string::npos);
}

void initializeKnownIndices(PrefixTree<string> * trieKnownString,string configFile){
    //p7 300
    string line;
    ifstream myFile;
    string content="";
    myFile.open(configFile.c_str(), ios::in);

    if (myFile.is_open()){
	while ( getline (myFile,line)){
	    content+=line;
	}
	myFile.close();
    }else{
	cerr << "Unable to open config file "<<configFile<<endl;
	exit(1);
    }


    JSONValue *value = JSON::Parse(content.c_str());
    if (value == NULL){
	cerr<<"Failed to parse JSON file"<<endl;
	exit(1);
    }

    JSONObject root;
    root = value->AsObject();
    if(root.find(L"indices") == root.end()){
	cerr<<"Failed to parse JSON file, needs a indices field"<<endl;
	exit(1);
    }
    

    JSONValue *jsonIndices    = root.at(L"indices");
    JSONObject jsonIndicesObj = jsonIndices->AsObject();


    //p7 indices
    if(jsonIndicesObj.find(L"p7indices") == jsonIndicesObj.end()){
	cerr<<"Failed to parse JSON file, needs a p7indices field"<<endl;
	exit(1);
    }
        
    JSONValue *jsonIndicesp7    = jsonIndicesObj.at(L"p7indices");
    JSONObject jsonIndicesp7Obj = jsonIndicesp7->AsObject();
    JSONArray arrayp7 = jsonIndicesp7Obj[L"p7index"]->AsArray();
    for (unsigned int i = 0; i < arrayp7.size(); i++){
	JSONObject temp=	arrayp7[i]->AsObject();
	
	if(temp.find(L"seq") == temp.end()){
	    cerr<<"Failed to parse JSON file, needs a seq in the p7indices field"<<endl;
	    exit(1);
	}
	if(temp.find(L"id") == temp.end()){
	    cerr<<"Failed to parse JSON file, needs a id in the p7indices field"<<endl;
	    exit(1);
	}


	string tempSeq (temp[L"seq"]->AsString().begin(),
			temp[L"seq"]->AsString().end());
	string tempID  (temp[L"id"]->AsString().begin(),
			temp[L"id"]->AsString().end());
	trieKnownString->insertIntoTree( tempSeq.c_str() , "p7#"    +tempID);
	trieKnownString->insertIntoTree( tempSeq.c_str() , "p7REVC#"+tempID);
    }


    //p5 indices
    if(jsonIndicesObj.find(L"p5indices") == jsonIndicesObj.end()){
	cerr<<"Failed to parse JSON file, needs a p5indices field"<<endl;
	exit(1);
    }
        
    JSONValue *jsonIndicesp5    = jsonIndicesObj.at(L"p5indices");
    JSONObject jsonIndicesp5Obj = jsonIndicesp5->AsObject();
    JSONArray arrayp5 = jsonIndicesp5Obj[L"p5index"]->AsArray();
    for (unsigned int i = 0; i < arrayp5.size(); i++){
	JSONObject temp=	arrayp5[i]->AsObject();
	
	if(temp.find(L"seq") == temp.end()){
	    cerr<<"Failed to parse JSON file, needs a seq in the p5indices field"<<endl;
	    exit(1);
	}
	if(temp.find(L"id") == temp.end()){
	    cerr<<"Failed to parse JSON file, needs a id in the p5indices field"<<endl;
	    exit(1);
	}


	string tempSeq (temp[L"seq"]->AsString().begin(),
			temp[L"seq"]->AsString().end());
	string tempID  (temp[L"id"]->AsString().begin(),
			temp[L"id"]->AsString().end());

	trieKnownString->insertIntoTree( tempSeq.c_str() , "p5#"    +tempID);
	trieKnownString->insertIntoTree( tempSeq.c_str() , "p5REVC#"+tempID);
    }




    //Other sequences after the primming site
    string IS4="AGATCTC";
    trieKnownString->insertIntoTree( IS4.c_str() ,                    "IS4");
    trieKnownString->insertIntoTree( reverseComplement(IS4).c_str() , "REVC#IS4");


}


//todo
//add fastq

void readIndexNumbers(string configFile){


    //p7 300
    string line;
    ifstream myFile;
    string content="";
    myFile.open(configFile.c_str(), ios::in);

    if (myFile.is_open()){
	while ( getline (myFile,line)){
	    content+=line;
	}
	myFile.close();
    }else{
	cerr << "Unable to open config file "<<configFile<<endl;
	exit(1);
    }


    JSONValue *value = JSON::Parse(content.c_str());
    if (value == NULL){
	cerr<<"Failed to parse JSON file"<<endl;
	exit(1);
    }

    JSONObject root;
    root = value->AsObject();
    if(root.find(L"indices") == root.end()){
	cerr<<"Failed to parse JSON file, needs a indices field"<<endl;
	exit(1);
    }
    

    JSONValue *jsonIndices    = root.at(L"indices");
    JSONObject jsonIndicesObj = jsonIndices->AsObject();


    //p7 indices
    if(jsonIndicesObj.find(L"p7indices") == jsonIndicesObj.end()){
	cerr<<"Failed to parse JSON file, needs a p7indices field"<<endl;
	exit(1);
    }
        
    JSONValue *jsonIndicesp7    = jsonIndicesObj.at(L"p7indices");
    JSONObject jsonIndicesp7Obj = jsonIndicesp7->AsObject();
    JSONArray arrayp7 = jsonIndicesp7Obj[L"p7index"]->AsArray();
    for (unsigned int i = 0; i < arrayp7.size(); i++){
	JSONObject temp=	arrayp7[i]->AsObject();
	
	if(temp.find(L"seq") == temp.end()){
	    cerr<<"Failed to parse JSON file, needs a seq in the p7indices field"<<endl;
	    exit(1);
	}
	if(temp.find(L"id") == temp.end()){
	    cerr<<"Failed to parse JSON file, needs a id in the p7indices field"<<endl;
	    exit(1);
	}


	string tempSeq (temp[L"seq"]->AsString().begin(),
			temp[L"seq"]->AsString().end());
	string tempID  (temp[L"id"]->AsString().begin(),
			temp[L"id"]->AsString().end());
	//cerr<<"p7 #"<<tempID<<"#\t#"<<tempSeq<<"#"<<endl;
	if(p7id2seq.find(tempID) == p7id2seq.end()){
	    p7id2seq[tempID] = tempSeq;
	}else{
	    cerr<<"Error while parsing JSON file, ID: "<<tempID<<" was found twice"<<endl;
	    exit(1);
	}
	
    }


    //p5 indices
    if(jsonIndicesObj.find(L"p5indices") == jsonIndicesObj.end()){
	cerr<<"Failed to parse JSON file, needs a p5indices field"<<endl;
	exit(1);
    }
        
    JSONValue *jsonIndicesp5    = jsonIndicesObj.at(L"p5indices");
    JSONObject jsonIndicesp5Obj = jsonIndicesp5->AsObject();
    JSONArray arrayp5 = jsonIndicesp5Obj[L"p5index"]->AsArray();
    for (unsigned int i = 0; i < arrayp5.size(); i++){
	JSONObject temp=	arrayp5[i]->AsObject();
	
	if(temp.find(L"seq") == temp.end()){
	    cerr<<"Failed to parse JSON file, needs a seq in the p5indices field"<<endl;
	    exit(1);
	}
	if(temp.find(L"id") == temp.end()){
	    cerr<<"Failed to parse JSON file, needs a id in the p5indices field"<<endl;
	    exit(1);
	}


	string tempSeq (temp[L"seq"]->AsString().begin(),
			temp[L"seq"]->AsString().end());
	string tempID  (temp[L"id"]->AsString().begin(),
			temp[L"id"]->AsString().end());

	if(p5id2seq.find(tempID) == p5id2seq.end()){
	    p5id2seq[tempID] = tempSeq;
	}else{
	    cerr<<"Error while parsing JSON file, ID: "<<tempID<<" was found twice"<<endl;
	    exit(1);
	}


    }






}


void printUnfoundToFile(vector< pair<string,int> > * unfound,ofstream & fileError,const bool  printErrorJSONabsent){

    for(int i=0;i<min(int(unfound->size()),maxErrorHits);i++){	       
	//Searching in known strings
	vector<string> temp = allTokens((*unfound)[i].first,'#');
	vector<string> temp2;
	for(unsigned int j=0;j<temp.size();j++){
	    vector<string> * temp3=new vector<string>();
	    vector<string> * temp4=new vector<string>();
	    if(!printErrorJSONabsent){
		trieKnownString->searchMismatch(temp[j].c_str(),temp3,0);
		trieKnownString->searchMismatch( ( "N"+temp[j].substr(0, temp[j].size() -1)     ).c_str(),temp4,1);
		trieKnownString->searchMismatch( (     temp[j].substr(1, temp[j].size() -1)+"N" ).c_str(),temp4,1);

		//adding a tag before the shifted ones
		for(unsigned int k=0;k<temp4->size();k++)
		    (*temp4)[k]="SHFT#"+(*temp4)[k];
		temp3->insert( temp3->end(), temp4->begin(), temp4->end() );
	    }
	    if(temp3->size() == 0 && temp4->size() == 0){
		temp2.push_back( "?");
	    }else{
		temp2.push_back( vectorToString(*temp3,","));  
	    }


	    delete temp3;
	    delete temp4;

	    
	}
	
	fileError<< 
	    vectorToString( temp,"\t" )<<"\t"<<
	    (*unfound)[i].second<<"\t"<<
	    vectorToString( temp2,"\t" )<< endl;
    }
}

void check_thresholds( rgAssignment &rg ) {
    rg.unknown  =    (-10.0 * rg.logLikelihoodScore)   > rgScoreCutoff ;
    rg.conflict =  (( -10.0 * rg.logRatioTopToSecond)  < fracConflict) && (rg.logRatioTopToSecond < 0.0) ;
    rg.wrong    =  (   10.0 * rg.topWrongToTopCorrect) > wrongness ;
}

// One BAM "work unit": either a single read (paired==false, only a1 is
// used) or a mate pair (paired==true, a1 & a2). The reader groups reads
// into units so that a pair is never split across a chunk boundary.
// Owns its htslib records; move-only so ownership transfers cleanly as the
// unit flows reader -> worker -> writer through the ordered pipeline.
struct ReadUnit {
    bam1_t *a1;
    bam1_t *a2;
    bool    paired;

    ReadUnit() : a1(NULL), a2(NULL), paired(false) {}
    ~ReadUnit(){ if(a1) bam_destroy1(a1); if(a2) bam_destroy1(a2); }

    ReadUnit(ReadUnit&& o) noexcept : a1(o.a1), a2(o.a2), paired(o.paired){ o.a1=NULL; o.a2=NULL; }
    ReadUnit& operator=(ReadUnit&& o) noexcept {
        if(this!=&o){
            if(a1) bam_destroy1(a1);
            if(a2) bam_destroy1(a2);
            a1=o.a1; a2=o.a2; paired=o.paired;
            o.a1=NULL; o.a2=NULL;
        }
        return *this;
    }
    ReadUnit(const ReadUnit&)            = delete;
    ReadUnit& operator=(const ReadUnit&) = delete;
};

// Assigns a read group and edits the tags of 'al' in place. Does NOT
// write the record: the caller (or the ordered pipeline's consume step)
// is responsible for saving it, so that output stays in input order.
void processSingleEndReads( bam1_t *al, bool printError, map<string,int> &unknownSeq, map<string,int> &wrongSeq, map<string,int> &conflictSeq,bool failBAM)
{
    string index1;
    string index1Q;
    string index2;
    string index2Q;

    getIndices(al,index1,index1Q,index2,index2Q);

    rgAssignment rgReturn=assignReadGroup(index1,index1Q,index2,index2Q,rgScoreCutoff,fracConflict,mismatchesTrie,qualOffset);
    check_thresholds( rgReturn ) ;

    updateRecord(al,rgReturn,failBAM);

    //record unresolved indices
    if(printError){
	string keyIndex;
	if(index2.empty()){
	    keyIndex=index1;
	}else{
	    keyIndex=index1+"#"+index2;
	}

	std::lock_guard<std::mutex> lk(statsMutex);
	if( rgReturn.conflict ) conflictSeq[ keyIndex ] ++;
	if( rgReturn.unknown  ) unknownSeq [ keyIndex ] ++;
	if( rgReturn.wrong    ) wrongSeq   [ keyIndex ] ++;
    }

}

// Same as above for a read pair; edits both 'al' and 'al2' in place.
void processPairedEndReads( bam1_t *al, bam1_t *al2, bool printError, map<string,int> &unknownSeq, map<string,int> &wrongSeq, map<string,int> &conflictSeq,bool failBAM)
{
    string index1;
    string index1Q;
    string index2;
    string index2Q;

    string sindex1;
    string sindex1Q;
    string sindex2;
    string sindex2Q;

    //retrieve indices
    getIndices(al,index1,index1Q,index2,index2Q);
    //check to see if the other indices are the same just for fun
    bool al2HasIndex=getIndices(al2,sindex1,sindex1Q,sindex2,sindex2Q,true);

    if(al2HasIndex && (index1 !=sindex1 )){cerr<<"Seq#1 has a different index 1 than seq #2, exiting "        <<bam_get_qname(al)<<" vs "<<bam_get_qname(al2)<< endl; fatalExit();}
    if(al2HasIndex && (index1Q!=sindex1Q)){cerr<<"Seq#1 has a different index 1 quality than seq #2, exiting "<<bam_get_qname(al)<<" vs "<<bam_get_qname(al2)<< endl; fatalExit();}
    if(al2HasIndex && (index2 !=sindex2 )){cerr<<"Seq#1 has a different index 2 than seq #2, exiting "        <<bam_get_qname(al)<<" vs "<<bam_get_qname(al2)<< endl; fatalExit();}
    if(al2HasIndex && (index2Q!=sindex2Q)){cerr<<"Seq#1 has a different index 2 quality than seq #2, exiting "<<bam_get_qname(al)<<" vs "<<bam_get_qname(al2)<< endl; fatalExit();}

    rgAssignment rgReturn = assignReadGroup(index1,index1Q,index2,index2Q,rgScoreCutoff,fracConflict,mismatchesTrie,qualOffset);
    check_thresholds( rgReturn ) ;

    updateRecord(al, rgReturn,failBAM);
    updateRecord(al2,rgReturn,failBAM);

    // record unresolved indices
    if(printError){
        string keyIndex;
        if(index2.empty()){
            keyIndex=index1;
        }else{
            keyIndex=index1+"#"+index2;
        }

        std::lock_guard<std::mutex> lk(statsMutex);
        if( rgReturn.conflict ) conflictSeq[ keyIndex ] += 1;
        if( rgReturn.unknown  ) unknownSeq [ keyIndex ] += 1;
        if( rgReturn.wrong    ) wrongSeq   [ keyIndex ] += 1;
    }
}


void checkFD(){
    int rofd  = returnOpenFileDescriptors();
    int rofdM = returnOpenFileDescriptorsMax();
    if(rofd == rofdM){
	cerr<<"WARNING: deML has detected that you have "<<rofd<<" open file descriptors out of a max. of "<<rofd<<endl;
	cerr<<"         therefore you have reached the maximum. If the information is correct, certain files might be empty"<<endl;
	cerr<<"         Either:"<<endl;
	cerr<<"           1) Use BAM as input/output"<<endl;
	cerr<<"           2) Check \"ulimit -n\" and put a higher number e.g. \"ulimit -n 1024\""<<endl;
	cerr<<"              If you already at the limit, increase the system limits"<<endl;
    }

}

// One fastq "work unit": the raw records read for a single template
// (forward, index1, and optionally index2 and reverse), plus the read
// group predicted for it by a worker thread.
struct FqUnit {
    string fId,  fSeq,  fQual;   // forward read
    string i1Id, i1Seq, i1Qual;  // first index
    string i2Id, i2Seq, i2Qual;  // second index  (only if hasId2Bool)
    string rId,  rSeq,  rQual;    // reverse read  (only if hasRevBool)

    string predictedGroup;        // filled by work()
    bool   fail;                  // filled by work(): failed a threshold
};

void processFastq(string           forwardfq,
		  string           reversefq,
		  string           index1fq,
		  string           index2fq,
		  string           prefixOut,
		  map<string,int> &unknownSeq,
		  map<string,int> &wrongSeq,
		  map<string,int> &conflictSeq,
		  const bool       printSummary,
		  const bool       printError){

    FastQParser * fqpf=0;
    FastQParser * fqpr=0;
    FastQParser * fqpi1=0;
    FastQParser * fqpi2=0;

    map<string,fqwriters *> rg2FqWriters;
    
    
    fqpf  = new FastQParser (forwardfq);
    fqpi1 = new FastQParser (index1fq);

    bool hasRevBool = (!reversefq.empty());
    bool hasId2Bool = (!index2fq.empty());


    if(hasRevBool)
	fqpr  = new FastQParser (reversefq);

    if(hasId2Bool)	
	fqpi2 = new FastQParser (index2fq);

    // ---- reader: pull one template's records from the parsers (single threaded) ----
    auto produce = [&](FqUnit &u) -> bool {
        if(!fqpf->hasData())
            return false;

        FastQObj * ffo = fqpf->getData();

        vector<string> deff = allTokens( *(ffo->getID()), ' '  );
        string deffs        = deff[0];

        if(!fqpi1->hasData()){
            cerr << "ERROR: Discrepancy between fastq files at record with first index " <<  *(ffo->getID()) <<endl;
            exit(1);
        }

        FastQObj * i1fo = fqpi1->getData();

        if(strEndsWith(deffs,  "/1")){
            deffs=deffs.substr(0,deffs.size()-2);
        }

        vector<string> defi1 = allTokens( *(i1fo->getID()), ' '  );
        string defi1s        = defi1[0];

        if( (deffs != defi1s ) ){
            cerr << "ERROR: Discrepancy between fastq files, different names with first index " <<deffs <<" and "<< defi1s <<endl;
            exit(1);
        }

        u.fId  = *(ffo->getID());   u.fSeq  = *(ffo->getSeq());   u.fQual  = *(ffo->getQual());
        u.i1Id = *(i1fo->getID());  u.i1Seq = *(i1fo->getSeq());  u.i1Qual = *(i1fo->getQual());

        if(hasId2Bool){
            if(!fqpi2->hasData()){
                cerr << "ERROR: Discrepancy between fastq files at record with second index " <<  *(ffo->getID()) <<endl;
                exit(1);
            }

            FastQObj * i2fo = fqpi2->getData();
            vector<string> defi2 = allTokens( *(i2fo->getID()), ' '  );
            string defi2s        = defi2[0];

            if( (deffs != defi2s ) ){
                cerr << "ERROR: Discrepancy between fastq files, different names with second index " <<deffs <<" and "<< defi2s <<endl;
                exit(1);
            }

            u.i2Id = *(i2fo->getID()); u.i2Seq = *(i2fo->getSeq()); u.i2Qual = *(i2fo->getQual());
        }

        if(hasRevBool){
            if(!fqpr->hasData()){
                cerr << "ERROR: Discrepency between fastq files at record with second read " <<  *(ffo->getID()) <<endl;
                exit(1);
            }

            FastQObj * rfo = fqpr->getData();
            vector<string> defr = allTokens( *(rfo->getID()), ' '  );
            string defrs        = defr[0];

            if(strEndsWith(defrs,  "/2")){
                defrs=defrs.substr(0,defrs.size()-2);
            }

            if( (deffs != defrs) ){
                cerr << "ERROR: Discrepency between fastq files, different names " <<deffs <<" and "<< defrs <<endl;
                exit(1);
            }

            u.rId = *(rfo->getID()); u.rSeq = *(rfo->getSeq()); u.rQual = *(rfo->getQual());
        }

        return true;
    };

    // ---- worker: assign the read group (heavy work, runs concurrently) ----
    auto work = [&](FqUnit &u){
        // assignReadGroup takes non-const references, so hand it local copies
        string i1s=u.i1Seq, i1q=u.i1Qual, i2s=u.i2Seq, i2q=u.i2Qual;
        rgAssignment rgReturn = assignReadGroup(i1s,i1q,i2s,i2q,rgScoreCutoff,fracConflict,mismatchesTrie,qualOffset);
        check_thresholds( rgReturn ) ;

        string predictedGroup = rgReturn.predictedGroup.empty() ? string("unknown") : rgReturn.predictedGroup;

        //////////////////////////
        //BEGIN update counters //
        //////////////////////////
        if(printSummary || printError){
            std::lock_guard<std::mutex> lk(statsMutex);

            if(printSummary){
                bool assigned=true;
                bool incrInTally=false;
                if( rgReturn.conflict ){  assigned=false; if(!incrInTally){ namesMap[predictedGroup].conflict++; incrInTally=true;}  }
                if( rgReturn.wrong    ){  assigned=false; if(!incrInTally){ namesMap[predictedGroup].wrong++;    incrInTally=true;}  }
                if( rgReturn.unknown  ){  assigned=false; if(!incrInTally){ namesMap[predictedGroup].unknown++;  incrInTally=true;}  }
                if(assigned){
                    namesMap[predictedGroup].assigned++;
                }
            }

            if(printError){
                string keyIndex;
                if(!hasId2Bool){
                    keyIndex=u.i1Seq;
                }else{
                    keyIndex=u.i1Seq+"#"+u.i2Seq;
                }

                if( rgReturn.conflict ) conflictSeq[ keyIndex ] += 1;
                if( rgReturn.unknown  ) unknownSeq [ keyIndex ] += 1;
                if( rgReturn.wrong    ) wrongSeq   [ keyIndex ] += 1;
            }
        }
        //////////////////////////
        //  END update counters //
        //////////////////////////

        u.predictedGroup = predictedGroup;
        u.fail           = ( rgReturn.conflict || rgReturn.wrong || rgReturn.unknown );
    };

    // ---- writer: open a read group's files on first sight, emit in input order ----
    auto consume = [&](FqUnit &u){
        const string & predictedGroup = u.predictedGroup;

        if(rg2FqWriters.find(predictedGroup) == rg2FqWriters.end()){ //new
            rg2FqWriters[predictedGroup] = new fqwriters();

            string outpairr1   = prefixOut+"_"+predictedGroup+"_"+"r1.fq.gz";
            string outpairr2   = prefixOut+"_"+predictedGroup+"_"+"r2.fq.gz";
            string outpairi1   = prefixOut+"_"+predictedGroup+"_"+"i1.fq.gz";
            string outpairi2   = prefixOut+"_"+predictedGroup+"_"+"i2.fq.gz";

            string outpairr1f   = prefixOut+"_"+predictedGroup+"_"+"r1.fail.fq.gz";
            string outpairr2f   = prefixOut+"_"+predictedGroup+"_"+"r2.fail.fq.gz";
            string outpairi1f   = prefixOut+"_"+predictedGroup+"_"+"i1.fail.fq.gz";
            string outpairi2f   = prefixOut+"_"+predictedGroup+"_"+"i2.fail.fq.gz";

            rg2FqWriters[predictedGroup]->pairr1.open( outpairr1.c_str(),  ios::out);
            rg2FqWriters[predictedGroup]->pairr1f.open(outpairr1f.c_str(), ios::out);

            rg2FqWriters[predictedGroup]->pairi1.open( outpairi1.c_str(),  ios::out);
            rg2FqWriters[predictedGroup]->pairi1f.open(outpairi1f.c_str(), ios::out);

            if(!rg2FqWriters[predictedGroup]->pairr1.good()){          checkFD(); cerr<<"Cannot write to file "<< outpairr1 <<" either you do not have permissions or you have too many read groups, in that case, convert your input data to a single BAM file and demultiplex it"<<endl; fatalExit(); }
            if(!rg2FqWriters[predictedGroup]->pairr1f.good()){         checkFD(); cerr<<"Cannot write to file "<< outpairr1f<<" either you do not have permissions or you have too many read groups, in that case, convert your input data to a single BAM file and demultiplex it"<<endl; fatalExit(); }

            if(!rg2FqWriters[predictedGroup]->pairi1.good()){          checkFD(); cerr<<"Cannot write to file "<< outpairi1 <<" either you do not have permissions or you have too many read groups, in that case, convert your input data to a single BAM file and demultiplex it"<<endl; fatalExit(); }
            if(!rg2FqWriters[predictedGroup]->pairi1f.good()){         checkFD(); cerr<<"Cannot write to file "<< outpairi1f<<" either you do not have permissions or you have too many read groups, in that case, convert your input data to a single BAM file and demultiplex it"<<endl; fatalExit(); }


            if(hasId2Bool){
                rg2FqWriters[predictedGroup]->pairi2.open( outpairi2.c_str(),  ios::out);
                rg2FqWriters[predictedGroup]->pairi2f.open(outpairi2f.c_str(), ios::out);

                if(!rg2FqWriters[predictedGroup]->pairi2.good()){      checkFD(); cerr<<"Cannot write to file "<< outpairi2 <<" either you do not have permissions or you have too many read groups, in that case, convert your input data to a single BAM file and demultiplex it"<<endl; fatalExit(); }
                if(!rg2FqWriters[predictedGroup]->pairi2f.good()){     checkFD(); cerr<<"Cannot write to file "<< outpairi2f<<" either you do not have permissions or you have too many read groups, in that case, convert your input data to a single BAM file and demultiplex it"<<endl; fatalExit(); }
            }

            if(hasRevBool){
                rg2FqWriters[predictedGroup]->pairr2.open( outpairr2.c_str(),  ios::out);
                rg2FqWriters[predictedGroup]->pairr2f.open(outpairr2f.c_str(), ios::out);

                if(!rg2FqWriters[predictedGroup]->pairr2.good()){      checkFD(); cerr<<"Cannot write to file "<<outpairr2 <<" either you do not have permissions or you have too many read groups, in that case, convert your input data to a single BAM file and demultiplex it"<<endl; fatalExit(); }
                if(!rg2FqWriters[predictedGroup]->pairr2f.good()){     checkFD(); cerr<<"Cannot write to file "<<outpairr2f<<" either you do not have permissions or you have too many read groups, in that case, convert your input data to a single BAM file and demultiplex it"<<endl; fatalExit(); }
            }
        }//end new rg


        if(u.fail){ //if the reads fail to meet thresholds

            rg2FqWriters[predictedGroup]->pairr1f     << u.fId  << "\n" << u.fSeq  << "\n+\n" << u.fQual  << "\n";
            rg2FqWriters[predictedGroup]->pairi1f     << u.i1Id << "\n" << u.i1Seq << "\n+\n" << u.i1Qual << "\n";
            if(hasId2Bool){
                rg2FqWriters[predictedGroup]->pairi2f << u.i2Id << "\n" << u.i2Seq << "\n+\n" << u.i2Qual << "\n";
            }

            if(hasRevBool){
                rg2FqWriters[predictedGroup]->pairr2f << u.rId  << "\n" << u.rSeq  << "\n+\n" << u.rQual  << "\n";
            }

        }else{
            rg2FqWriters[predictedGroup]->pairr1     << u.fId  << "\n" << u.fSeq  << "\n+\n" << u.fQual  << "\n";
            rg2FqWriters[predictedGroup]->pairi1     << u.i1Id << "\n" << u.i1Seq << "\n+\n" << u.i1Qual << "\n";
            if(hasId2Bool){
                rg2FqWriters[predictedGroup]->pairi2 << u.i2Id << "\n" << u.i2Seq << "\n+\n" << u.i2Qual << "\n";
            }

            if(hasRevBool){
                rg2FqWriters[predictedGroup]->pairr2 << u.rId  << "\n" << u.rSeq  << "\n+\n" << u.rQual  << "\n";
            }
        }
    };

    runOrderedPipeline<FqUnit>(numberOfThreads, chunkSizeUnits, produce, work, consume);
    
    //map<string,fqwriters *>::iterator rg2FqwritersIt;
    vector<string> allRgsFound  = allKeysMap(rg2FqWriters);
    for (unsigned int i=0; i<allRgsFound.size();i++){
	string rgTag = allRgsFound[i];

	rg2FqWriters[rgTag]->pairr1.close( );
	rg2FqWriters[rgTag]->pairr1f.close();
	
	rg2FqWriters[rgTag]->pairi1.close( );
	rg2FqWriters[rgTag]->pairi1f.close();
	
	if(hasId2Bool){
	    rg2FqWriters[rgTag]->pairi2.close( );
	    rg2FqWriters[rgTag]->pairi2f.close();
	}
	
	if(hasRevBool){
	    rg2FqWriters[rgTag]->pairr2.close( );
	    rg2FqWriters[rgTag]->pairr2f.close();	   
	}
    }
    
    delete fqpf;
    delete fqpi1;
 
   if(hasId2Bool)
	delete fqpi2;

    if(hasRevBool)
	delete fqpr;

}



int main (int argc, char *argv[]) {

    string bamFile;
    string bamFileOUT="";

    string index="";
    string outfile;
    bool   printSummary=false;
    string filenameSummary;

    bool   printError=false;
    bool   printErrorJSONabsent=false;

    string filenameError;

    ofstream ratioValuesOS;
    ofstream rgqualOS;

    bool ratioValuesFlag = false;
    bool rgqualFlag      = false;
    bool shiftByOne      = false;
    bool failBAM         = true;

    bool useFastq=false;
    string forwardfq;
    string reversefq;
    string index1fq;
    string index2fq;

    bool produceUnCompressedBAM=false; 

    const string usage=string(string(argv[0])+
			      " [options] BAMfile"+"\n\n"+

			      "\t"+"Mandatory parameters:"+"\n"+
			      "\t\t"+"-i"+","+"--index"+"\t[index]"+"\t\t\t"+"File describing index sequences used"+"\n"+
			      "\t\t"+"-o"+","+"--outfile"+"\t[outfile]"+"\t\t"+"Specify output file"+"\n\n"+

			      "\t"+"BAM input options:"+"\n"+
			      "\t\t"+"--bamtags"+"\t[Idx1Seq],[Idx1Qual],[Idx2Seq],[Idx2Qual]"+"\t"+"BAM Tags containing barcodes\n"+

			      "\t"+"Fastq input/output (optional):"+"\n"+
			      "\t\t"+"You can specify fastq as input/output, in which case the -o option will be"+"\n"+
			      "\t\t"+"treated as an output prefix"+"\n"+
			      
			      "\t\t"+"-f"  +"\t[forward fastq]"+"\t\t\t"+""+"Forward reads in fastq\n"+
			      "\t\t"+"-r"  +"\t[reverse fastq]"+"\t\t\t"+""+"Reverse reads in fastq\n"+
			      "\t\t"+"-if1"  +"\t[index1 file fastq]"+"\t\t"+""+"First index sequences in fastq\n"+
			      "\t\t"+"-if2"  +"\t[index2 file fastq]"+"\t\t"+""+"Second index sequences in fastq\n"+
			      
			      "\tCutoffs options:"+"\n"
			      "\t\t"+"--phred64"  +"\t"+"\t\t\t"+""+"Use PHRED 64 as the offset for QC scores (default : PHRED33)"+"\n"+

			      "\t\t"+"--rgqual"  +"\t[quality]"+"\t\t"+""+"Worst quality before flagging as unknown ["+stringify(rgScoreCutoff)+"]\n"+
			      "\t\t"+"--fracconf"+"\t[quality]"+"\t\t"+""+"Maximum quality difference before being considered a conflict ["+stringify(fracConflict)+"] \n"+
			      "\t\t"+"--wrongness"+"\t[quality]"+"\t\t"+""+"Mininum quality difference to flag as wrongly paired ["+stringify(wrongness)+"] \n"+
			      "\t\t"+"--mm"+"\t\t[mismatches]"+"\t\t"+""+"Maximum # of tolerated mismatches ["+stringify(mismatchesTrie)+"] \n"+

			      "\n\tRG assignment options:"+"\n"+
			      "\t\t"+"" +""+"--shift"+"\t"+"\t\t\t\t"+"Try shifting the index right by one at the cost of a mismatch"+"\n"+

			      
			      //"\t"+"\tOptional:"+"\n"+
			      "\n\tPerformance options:"+"\n"+
			      "\t\t"+"-t"+"\t\t[threads]"+"\t\t"+""+"Number of worker threads, output stays in input order (default : "+stringify(numberOfThreads)+")\n"+

			      "\n\tOutput options:"+"\n"+
			      "\t\t"+"--maxerr"+"\t[max err]"+"\t\t"+""+"Print  # wrongly of assigned RG in the error log (--error) ["+stringify(maxErrorHits)+"] \n"+
                              "\t\t"+"-u" +"\t\t\t\t\t"           +"Produce uncompressed bam (good for pipe)"+"\n"+ 
			      "\t\t"+"-s"+","+"--summary"+"\t[summary file]"+"\t\t"+"Summarize the RG tally in this file"+"\n"+
			      "\t\t"+"-e"+","+"--error"  +"\t[error file]"+"\t\t"+"Summarize the indices that were not assigned to a RG"+"\n"+
			      "\t\t"+""+""+"--rgval"  +"\t[file]"+"\t\t\t\t"+"Write the rg qualities as a binary file"+"\n"+
			      "\t\t"+""+""+"--ratio"   +"\t\t[file]"+"\t\t\t"+"Write the likelihood ratios as a binary file"+"\n"
			      "\t\t"+""+""+"--nofail"   +"\t\t\t"+"\t"+"Do not set the QC fail in the output BAM file"+"\n"

			      );
			      

    if( (argc== 1) ||
    	(argc== 2 && string(argv[1]) == "-h") ||
    	(argc== 2 && string(argv[1]) == "-help") ||
    	(argc== 2 && string(argv[1]) == "--help") ){
    	cout<<"Usage:"<<endl;
    	cout<<""<<endl;
    	cout<<usage<<endl;
    	return 1;
    }

    int lastIndexArgc = argc-1;
    for(int i=1;i<(lastIndexArgc);i++){ //all but the last
	
	if(string(argv[i]) == "--nofail"  ){
            failBAM=false;
            continue;
        }

	if(strcmp(argv[i],"-f") == 0 ){
	    forwardfq = string(argv[i+1]);
	    useFastq=true;
	    lastIndexArgc = argc;
	    i++;
	    continue;
	}

	if(strcmp(argv[i],"-r") == 0 ){
	    reversefq = string(argv[i+1]);
	    useFastq=true;
	    lastIndexArgc = argc;
	    i++;
	    continue;
	}

	if(strcmp(argv[i],"--bamtags") == 0 ){
	    istringstream iss(argv[i+1]);
	    string tok;
	    int p=0;
	    while(getline(iss, tok, ',')) {
	        switch (p) {
	            case 0: tagIndex1Seq = tok; break;
	            case 1: tagIndex1Qual = tok; break;
	            case 2: tagIndex2Seq = tok; break;
	            case 3: tagIndex2Qual = tok; break;
	            default:
	                cerr<<"too many commas in "<<argv[i+1]<<endl;
	                exit(1);
	        }
	        ++p;
	    }
	    i++;
	    continue;
	}

	if(strcmp(argv[i],"-if1") == 0 ){
	    index1fq = string(argv[i+1]);
	    useFastq=true;
	    lastIndexArgc = argc;
	    i++;
	    continue;
	}

	if(string(argv[i]) == "--phred64"  ){
            qualOffset=64;
            continue;
        }

	if(strcmp(argv[i],"-if2") == 0 ){
	    index2fq = string(argv[i+1]);
	    useFastq=true;
	    lastIndexArgc = argc;
	    i++;
	    continue;
	}
	

	if(strcmp(argv[i],"--shift") == 0 ){
	    shiftByOne      = true;
	    continue;
	}

	if(strcmp(argv[i],"--rgval") == 0 ){
	    string temp =string(argv[i+1]);
	    rgqualOS.open(temp.c_str(), ios::out | ios::binary);
	    rgqualFlag      = true;
	    if (!rgqualOS){
		cerr<<"Cannot print to file "<<temp<<endl;
		exit(1);
	    }
	    setFileForRGQual(&rgqualOS);
	    i++;
	    continue;
	}


	if(strcmp(argv[i],"-u") == 0  ){ 
	    produceUnCompressedBAM=true; 
	    continue; 
	} 
	

	if(strcmp(argv[i],"--ratio") == 0 ){
	    string temp =string(argv[i+1]);
	    ratioValuesOS.open(temp.c_str(), ios::out | ios::binary);
	    ratioValuesFlag = true;

	    if (!ratioValuesOS){
		cerr<<"Cannot print to file "<<temp<<endl;
		exit(1);
	    }
	    setFileForRatio(&ratioValuesOS);
	    i++;
	    continue;
	}

	if(strcmp(argv[i],"-e") == 0 || strcmp(argv[i],"--error") == 0 ){
	    printError=true;
	    filenameError =string(argv[i+1]);
	    i++;
	    continue;
	}



	if(strcmp(argv[i],"-s") == 0 || strcmp(argv[i],"--summary") == 0 ){
	    printSummary=true;
	    filenameSummary =string(argv[i+1]);
	    i++;
	    continue;
	}





	if(strcmp(argv[i],"--maxerr") == 0 ){
	    maxErrorHits =destringify<int>(argv[i+1]);
	    i++;
	    continue;
	}


	if(strcmp(argv[i],"--rgqual") == 0 ){
	    rgScoreCutoff =destringify<double>(argv[i+1]);
	    i++;
	    continue;
	}



	if(strcmp(argv[i],"--fracconf") == 0 ){
	    fracConflict =destringify<double>(argv[i+1]);
	    i++;
	    continue;
	}

	if(strcmp(argv[i],"--wrongness") == 0 ){
	    wrongness =destringify<double>(argv[i+1]);
	    i++;
	    continue;
	}

	if(strcmp(argv[i],"--mm") == 0 ){
	    mismatchesTrie =destringify<int>(argv[i+1]);
	    i++;
	    continue;
	}

	if(strcmp(argv[i],"-t") == 0 ){
	    numberOfThreads =destringify<int>(argv[i+1]);
	    if(numberOfThreads < 1){
		cerr<<"The number of threads (-t) must be at least 1"<<endl;
		return 1;
	    }
	    i++;
	    continue;
	}



	if(strcmp(argv[i],"-i") == 0 || strcmp(argv[i],"--index") == 0 ){
	    index =string(argv[i+1]);
	    i++;
	    continue;
	}

	if(strcmp(argv[i],"-o") == 0 || strcmp(argv[i],"--outfile") == 0 ){
	    outfile=string(argv[i+1]);
	    i++;
	    continue;
	}


	cerr<<"Unknown option "<<argv[i] <<" exiting"<<endl;
	return 1;             
    }



    if(useFastq){

	if(forwardfq.empty() ){
	    cerr<<"If fastq is used, the forward read must be specified"<<endl;
	    return 1;             
	}

	if(index1fq.empty() ){
	    cerr<<"If fastq is used, the forward read must be specified"<<endl;
	    return 1;             
	}

	if(produceUnCompressedBAM){
	    cerr<<"If fastq is used, cannot use -u flag"<<endl;
	    return 1;             
	}
	
    }
    
    
    
    if(outfile.size() == 0){
	cerr<<"The field -o is mandatory exiting"<<endl;
	return 1;             
    }

    if(index.size() == 0){
	cerr<<"The field -i is mandatory exiting"<<endl;
	return 1;
    }

    // The --ratio / --rgval binary dumps are written from inside
    // assignReadGroup() as reads are processed. Under multiple threads
    // their order would be non-deterministic (and the shared stream would
    // be corrupted), so fall back to a single thread when they are used.
    if( (ratioValuesFlag || rgqualFlag) && numberOfThreads > 1 ){
	cerr<<"WARNING: --ratio/--rgval write per-read values in input order; forcing -t 1"<<endl;
	numberOfThreads = 1;
    }


    if(printError){
	trieKnownString = new PrefixTree<string>();

	string configFile = "/../../webForm/config.json";

	if(!isFile( getCWD(argv[0])+configFile )){
	    configFile = "/../webForm/config.json";
	    if(!isFile( getCWD(argv[0])+configFile )){
		cerr<<"ERROR: file "<<(getCWD(argv[0])+"/../../webForm/config.json")<<" or "<<(getCWD(argv[0])+"/../webForm/config.json")<<" were not found, this file is required with the -e or --error option"<<endl;
		//return 1;             
	    }
	    printErrorJSONabsent=true;
	}else{
	    initializeKnownIndices(trieKnownString,getCWD(argv[0])+configFile);
	}
	// //debug
	// vector<string> * temp3=new vector<string>();
	// vector<string> * temp4=new vector<string>();
	// string t="AAGGTCT";
	// trieKnownString->searchMismatch(t.c_str(),temp3,0);
	// trieKnownString->searchMismatch( ( "N"+t.substr(0, t.size() -1)     ).c_str(),temp4,1);
	// trieKnownString->searchMismatch( (     t.substr(1, t.size() -1)+"N" ).c_str(),temp4,1);
	// cout<<"t3 "<<vectorToString(*temp3)<<endl;
	// cout<<"t4 "<<vectorToString(*temp4)<<endl;

	// return 1;
	// //end debug
    }

    bamFile=argv[argc-1];
    ifstream myIndexFile;


    //
    //  BEGIN : reading index file
    //    
    string indexStringFile="";
    bool firstLine=true;
    bool fileContainsSeq=true;
    bool fileContainsSeqDouble=false;

    myIndexFile.open(index.c_str(), ios::in);
    if (myIndexFile.is_open()){
	string line;
	while ( getline (myIndexFile,line)){
	    indexStringFile+=line+"\n";

	    if(firstLine){
		firstLine=false;
	    }else{
		vector<string> tempf = splitWS(line);

		if(tempf.size() == 2){ //contains 2 fields, single index
		    fileContainsSeqDouble=false;

		    if(isStringNatNumber(    tempf[0]) || 
		       isStringNatNumber("t"+tempf[0])){ //for truseq
			fileContainsSeq=false;
		    }

		}else{
		    if(tempf.size() == 3){ //contains 3 fields, double index
			fileContainsSeqDouble=true;

			if(isStringNatNumber(    tempf[0]) || 
			   isStringNatNumber("t"+tempf[0])){ //for truseq
			    fileContainsSeq=false;

			    if( isStringNatNumber(    tempf[1]) || 
			        isStringNatNumber("t"+tempf[1]) ){ //for truseq
				//fine
			    }else{
				cerr << "Error: line "<<line<<" cannot have a mixture of index sequence and index numbers"<<endl;
				exit(1);				
			    }

			}

		    }else{
			cerr << "Error: line "<<line<<" does not have 2 or 3 fields"<<endl;
			exit(1);
		    }

		}

	    }
	}
	myIndexFile.close();
    }else{ 
	cerr << "Unable to open file "<<index<<endl;
	exit(1);
    }
    //
    //  END : reading index file
    //    


    //
    //  BEGIN : Processing index with numbers
    //
    if(!fileContainsSeq){


	string configFile = "/../../webForm/config.json";

	if(!isFile( getCWD(argv[0])+configFile )){
	    configFile = "/../webForm/config.json";
	    if(!isFile( getCWD(argv[0])+configFile )){
		cerr<<"ERROR: file "<<(getCWD(argv[0])+"/../../webForm/config.json")<<" or "<<(getCWD(argv[0])+"/../webForm/config.json")<<" were not found, this file is required with the -e or --error option"<<endl;
		return 1;             
	    }
	}
	//initializeKnownIndices(trieKnownString,getCWD(argv[0])+configFile);
	
	readIndexNumbers( getCWD(argv[0])+configFile);	

	firstLine=true;
	vector<string> allLinesIndexFile = allTokens(indexStringFile,'\n');
	string tempIndexStringFile="";
	//while ( getline (myFile,line)){
	for(unsigned int i=0;i<allLinesIndexFile.size();i++){
	    string line = allLinesIndexFile[i];
	    if(line.empty())
		continue;
	    vector<string> tempf = splitWS(line);

	    if( ((tempf.size() == 2) && !fileContainsSeqDouble) || 
		((tempf.size() == 3) && fileContainsSeqDouble) ){
	    }else{
		cerr << "Error: line "<<line<<" does not have 2 or 3 fields"<<endl;
		exit(1);
	    }

	    if(firstLine){
		tempIndexStringFile += line;
		firstLine=false;
	    }else{

		if(p7id2seq.find(tempf[0]) != p7id2seq.end()){
		    tempIndexStringFile +=  p7id2seq[ tempf[0] ] ;
		}else{
		    cerr<<"Error p7 ID: #"<<tempf[0]<<"# not found"<<endl;
		    exit(1);
		}
	    
		if(fileContainsSeqDouble){

		    if(p5id2seq.find(tempf[1]) != p5id2seq.end()){
			tempIndexStringFile +=  "\t"+p5id2seq[ tempf[1] ] ;
		    }else{
			cerr<<"Error p5 ID: #"<<tempf[1]<<"# not found"<<endl;
			//cerr<<"Error while parsing JSON file, ID: "<<tempID<<" was found twice"<<endl;
			exit(1);
		    }
		
		    tempIndexStringFile +=  "\t"+tempf[2]  ;
		}else{
		    tempIndexStringFile +=  "\t"+tempf[1]  ;
		}
	    }
	    
	    tempIndexStringFile +=  "\n";
	    //cerr<<"temp "<<tempIndexStringFile<<endl;
	}//for all lines
	
	indexStringFile=tempIndexStringFile;
	cerr<<"Using the following as index sequences"<<endl;
	cerr<<"----------------------------------"<<endl;
	cerr<<indexStringFile<<endl;
	cerr<<"----------------------------------"<<endl;
    }//end file contains seq
    //
    //  END : Processing index with numbers
    //



    map<string,string> rgs =readIndexFile(indexStringFile,mismatchesTrie,shiftByOne);
    map<string,int> unknownSeq;
    map<string,int> wrongSeq;
    map<string,int> conflictSeq;
    
    if(useFastq){

	processFastq(forwardfq,
		     reversefq,
		     index1fq,
		     index2fq,
		     outfile,
		     unknownSeq,
		     wrongSeq,
		     conflictSeq,
		     printSummary,
		     printError);

    }else{


	samFile *in = sam_open(bamFile.c_str(), "r");
	if ( !in ) {
	    cerr << "Could not open input BAM file  "<<bamFile << endl;
	    return 1;
	}
	// multithreaded BGZF decompression
	if( numberOfThreads > 1 ) hts_set_threads(in, numberOfThreads);

	sam_hdr_t *hdr = sam_hdr_read(in);
	if( !hdr ){
	    cerr << "Could not read header from input BAM file  "<<bamFile << endl;
	    return 1;
	}

	string pCommandLine = "";
	for(int i=0;i<(argc);i++){
	    pCommandLine += (string(argv[i])+" ");
	}
	// add our @PG line (htslib chains PP to the previous program automatically)
	if( sam_hdr_add_pg(hdr, "deML",
			   "PN", "deML",
			   "CL", pCommandLine.c_str(),
			   "VN", returnGitHubVersion(string(argv[0]),"..").c_str(),
			   NULL) < 0 ){
	    cerr << "Warning: could not add @PG line to the BAM header" << endl;
	}

	// Replace any read groups already present in the input header with the
	// ones from the index file (mirrors the old bamtools "ReadGroups=srgd").
	sam_hdr_remove_except(hdr, "RG", NULL, NULL);

	map<string,string>::const_iterator itRG;
	for ( itRG=rgs.begin(); itRG != rgs.end(); itRG++ ){

	    if(itRG->first == "conflict" || itRG->first == "unknown" || itRG->first == "wrong" ){
		cerr<<"ERROR: The RG names cannot contain the words: \"conflict\" or \"unknown\" or \"wrong\""<<endl;
		return 1;
	    }

	    // add @RG  ID:<name>  [DS:<description>]
	    int rgrc;
	    if( itRG->second.empty() )
		rgrc = sam_hdr_add_line(hdr, "RG", "ID", itRG->first.c_str(), NULL);
	    else
		rgrc = sam_hdr_add_line(hdr, "RG", "ID", itRG->first.c_str(), "DS", itRG->second.c_str(), NULL);
	    if( rgrc < 0 )
		cerr<<"Warning: could not add read group "<<itRG->first<<" to the BAM header"<<endl;

	    namesMap[ itRG->first ].assigned =0;
	    namesMap[ itRG->first ].unknown  =0;
	    namesMap[ itRG->first ].conflict =0;
	    namesMap[ itRG->first ].wrong    =0;
	}

	namesMap[ "unknown" ].assigned=0;
	namesMap[ "unknown" ].unknown=0;
	namesMap[ "unknown" ].conflict=0;
	namesMap[ "unknown" ].wrong=0;

	samFile *out = sam_open(outfile.c_str(), produceUnCompressedBAM ? "wb0" : "wb");
	if( !out ) {
	    cerr << "Could not open output BAM file  "<<outfile << endl;
	    return 1;
	}
	// multithreaded BGZF compression
	if( numberOfThreads > 1 ) hts_set_threads(out, numberOfThreads);

	if( sam_hdr_write(out, hdr) < 0 ){
	    cerr << "Could not write header to output BAM file  "<<outfile << endl;
	    return 1;
	}


	// Reader state for the look-ahead pairing below. 'pending' holds a
	// record that was read while probing for a mate but turned out to
	// belong to the next unit. Records are recycled by swapping pointers,
	// so no BAM record is ever deep-copied.
	bam1_t *pending     = bam_init1();
	bool    havePending = false;
	bool    readerEOF   = false;

	// produce(): reproduce the original convoluted pairing logic, but
	// emit one ReadUnit (single or pair) at a time instead of processing
	// it directly. A pair therefore never straddles a chunk boundary.
	auto produce = [&](ReadUnit &u) -> bool {
	    if( !u.a1 ) u.a1 = bam_init1();
	    if( havePending ){
		std::swap(u.a1, pending);   // take the carried-over record
		havePending = false;
	    }else{
		if( readerEOF ) return false;
		if( sam_read1(in, hdr, u.a1) < 0 ){ readerEOF = true; return false; }
	    }

	    while(1) {
		if( !u.a2 ) u.a2 = bam_init1();
		if( sam_read1(in, hdr, u.a2) < 0 ) {
		    // EOF: u.a1 is a leftover single
		    readerEOF = true;
		    u.paired  = false;
		    return true;
		}
		// If it's paired, both should have the same index, and we
		// save some work.  Since the reads are probably not
		// ordered, check the names first
		if( (u.a1->core.flag & BAM_FPAIRED) &&
		    strcmp(bam_get_qname(u.a1), bam_get_qname(u.a2))==0 ) {
		    u.paired = true;            // u.a2 already holds the mate
		    return true;
		} else {
		    // no match, treat u.a1 as a single and carry u.a2 over
		    u.paired    = false;
		    std::swap(pending, u.a2);   // stash the look-ahead record
		    havePending = true;
		    return true;
		}
	    }
	};

	// work(): the heavy read-group assignment, thread-safe, in place.
	auto work = [&](ReadUnit &u){
	    if( u.paired )
		processPairedEndReads(u.a1,u.a2,printError,unknownSeq,wrongSeq,conflictSeq,failBAM);
	    else
		processSingleEndReads(u.a1,printError,unknownSeq,wrongSeq,conflictSeq,failBAM);
	};

	// consume(): save the processed records in input order.
	auto consume = [&](ReadUnit &u){
	    if( sam_write1(out, hdr, u.a1) < 0 ){ cerr<<"Error writing BAM record"<<endl; fatalExit(); }
	    if( u.paired )
		if( sam_write1(out, hdr, u.a2) < 0 ){ cerr<<"Error writing BAM record"<<endl; fatalExit(); }
	};

	runOrderedPipeline<ReadUnit>(numberOfThreads, chunkSizeUnits, produce, work, consume);

	bam_destroy1(pending);
	sam_hdr_destroy(hdr);
	if( sam_close(in)  < 0 ){ cerr << "Error closing input BAM file  "<<bamFile << endl; }
	if( sam_close(out) < 0 ){ cerr << "Error closing output BAM file "<<outfile<< endl; return 1; }






    }//bam
    //cleaning up

    //Print summary of RG assignment
    if(printSummary){
	map<string,tallyForRG>::iterator it;   
	unsigned int totalRG=0;	
	unsigned int totalAssignRG=0;	

	vector< pair<string,tallyForRG> > toprintVec;
	for ( it=namesMap.begin() ; it != namesMap.end(); it++ ){
	    toprintVec.push_back(  make_pair( it->first , it->second ) );
	    totalRG+=it->second.assigned+it->second.unknown+it->second.conflict+it->second.wrong;
	}

	sort (toprintVec.begin(),   toprintVec.end(),   compareNameTally() ); 
	ofstream fileSummary;
	fileSummary.open(filenameSummary.c_str());

	if (fileSummary.is_open()){

	    fileSummary << "RG\ttotal\ttotal%\tassigned\tassigned%\tunknown\tunknown%\tconflict\tconflict%\twrong\twrong%"<<endl;
	    fileSummary<<dashes<<endl;
	    for(unsigned int i=0;i<toprintVec.size();i++){		
		unsigned int totalForRQ=toprintVec[i].second.assigned+toprintVec[i].second.unknown+toprintVec[i].second.conflict+toprintVec[i].second.wrong;

		fileSummary << toprintVec[i].first << "\t" << totalForRQ << "\t"
			    << 100.0*double(totalForRQ)/double(totalRG) << "%\t" ;

		fileSummary  << toprintVec[i].second.assigned << "\t"
			     << 100.0*double(toprintVec[i].second.assigned)/double(totalForRQ) << "%\t" ;

		fileSummary  << toprintVec[i].second.unknown << "\t"
			     << 100.0*double(toprintVec[i].second.unknown)/double(totalForRQ) << "%\t" ;

		fileSummary <<  toprintVec[i].second.conflict << "\t"
			    << 100.0*double(toprintVec[i].second.conflict)/double(totalForRQ) << "%\t" ;

		fileSummary <<  toprintVec[i].second.wrong << "\t"
			    << 100.0*double(toprintVec[i].second.wrong)/double(totalForRQ) << "%\n" ;
		
		if(toprintVec[i].first != "unknown" )
		    totalAssignRG+=toprintVec[i].second.assigned;
	    }

	    fileSummary<<dashes<<endl;
	    fileSummary<<"ASSIGNED:\t"<< totalAssignRG<<"\t"<<100.0*double(totalAssignRG)/double(totalRG)<<"%"<<endl;
	    fileSummary<<"PROBLEMS:\t"<< (totalRG-totalAssignRG)<<"\t"<<100.0*double(totalRG-totalAssignRG)/double(totalRG)<<"%"<<endl;

	    fileSummary<<"TOTAL:\t"<<totalRG<<"\t100.0%"<<endl;
	}else{
	    cerr << "Unable to print to file "<<filenameSummary<<endl;
	}
	fileSummary.close();
    }




    //Print over-represented sequences in conflict,unknown,wrong
    if(printError){
	vector< pair<string,int> > conflictToPrint( conflictSeq.begin(), conflictSeq.end() ) ;
	vector< pair<string,int> > unknownToPrint(  unknownSeq.begin(),  unknownSeq.end() ) ;
	vector< pair<string,int> > wrongToPrint(    wrongSeq.begin(),    wrongSeq.end() ) ;
     	
	sort (conflictToPrint.begin(),   conflictToPrint.end(),   compareNameRG() ); 
	sort (unknownToPrint.begin(),    unknownToPrint.end(),    compareNameRG() ); 
	sort (wrongToPrint.begin(),      wrongToPrint.end(),      compareNameRG() ); 

	ofstream fileError;
	fileError.open(filenameError.c_str());
	if (fileError.is_open())
	    {
		fileError<<      dashes<<endl<<"Conflict:"<<endl<<dashes<<endl;
		printUnfoundToFile(&conflictToPrint,fileError, printErrorJSONabsent);

		fileError<<endl<<dashes<<endl<<"Unknown:" <<endl<<dashes<<endl;
		printUnfoundToFile(&unknownToPrint,fileError, printErrorJSONabsent);

		fileError<<endl<<dashes<<endl<<"Wrong:"   <<endl<<dashes<<endl;
		printUnfoundToFile(&wrongToPrint,fileError, printErrorJSONabsent);
	    }else{
	    cerr << "Unable to print to file "<<filenameError<<endl;
	}
	fileError.close();
    }

    deallocate();
    if(printError){
	delete trieKnownString;
    }

    if(rgqualFlag)
	rgqualOS.close();

    if(ratioValuesFlag)
	ratioValuesOS.close();
    
    return 0;
}

