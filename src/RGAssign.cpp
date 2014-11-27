#include "RGAssign.h"
// vim:ts=8

// #define DEBUG
// #define DEBUG2
//#define DEBUG3
// #define DEBUGMISPAIR

PrefixTree<int> * indTrie1;
PrefixTree<int> * indTrie2;
indexData values;               

ofstream * ratioValues;
ofstream * rgqual;
bool flag_ratioValues=false;
bool flag_rgqual=false;
bool warnedForDouble=false;


// array containing the log10 of likelihoods corresponding to quality
// scores (these are negative numbers)
double likeMatch[64];
double likeMismatch[64];

bool shiftByOne=false;

void setFileForRatio(ofstream * streamFile){
    flag_ratioValues=true;
    ratioValues=streamFile;
}

void setFileForRGQual(ofstream * streamFile){
    flag_rgqual=true;
    rgqual=streamFile;
}


// Returns log10( pow(10,x)+pow(10,y) ), but does so without causing
// overflow or loss of precision.
// double oplus( double x, double y )
// {
//     return x > y 
//         ? x + log1p( pow( 10, y-x ) ) / log(10)
//         : y + log1p( pow( 10, x-y ) ) / log(10) ;
// }


// compare by likelihood (second component).  since likelihoods are
// negative, the meaning is flipped
struct comparePair 
{
    bool operator() (pair<int,double> i,pair<int,double> j) { 
        return (i.second>j.second); 
    }
} ;

string toUpperCase(string toCheck){
    string toReturn=toCheck;
    transform(toReturn.begin(), toReturn.end(), toReturn.begin(), ::toupper);
    return toReturn;
}

bool isValidDNA(string tocheck){
    for(unsigned int  i=0;i<tocheck.size();i++){
	if(tocheck[i] != 'A' &&
	   tocheck[i] != 'C' &&
	   tocheck[i] != 'G' &&
	   tocheck[i] != 'T' )
	    return false;	  
    }
    return true;
}

void checkRGname(string tocheck){
    if(tocheck == "unknown") {  cerr<<"Error: Cannot use reserved name \"unknown\" as read group name"<<endl;  exit(1);}
    if(tocheck == "conflict"){  cerr<<"Error: Cannot use reserved name \"conflict\" as read group name"<<endl; exit(1);}
    if(tocheck == "wrong")   {  cerr<<"Error: Cannot use reserved name \"wrong\" as read group name"<<endl;    exit(1);}
    //fine
}

indexData intern_readIndex(string filename){
    string line;
    ifstream myFile;

    indexData toReturn;
    toReturn.mlindex1=0;
    toReturn.mlindex2=0;

    bool isFirstLine  =true;


    //initialize the values for the likelihood of matches or mismatches 
    for(int i=0;i<64;i++){
	if(i == 0)
	    likeMatch[i]    = -3.0; // this is vrong, hope it's never accessed
	else
	    likeMatch[i]    = log1p( -pow(10.0,i/-10.0) )/log(10);
	
	likeMismatch[i]     = i/-10.0;	
#ifdef DEBUG2
	cout<<"qual = "<<i<<endl;
	cout<<likeMatch[i]<<endl;
	cout<<likeMismatch[i]<<endl;
#endif
    }

    //reading the files



    // myFile.open(filename.c_str(), ios::in);
    // if (myFile.is_open()){

    vector<string> allLinesIndex = allTokens(filename,'\n');

    //while ( getline (myFile,line)){
    for(unsigned int i=0;i<allLinesIndex.size();i++){
	line = allLinesIndex[i];
	if(line.empty())
	    continue;
	line+=' ';
	// cerr<<"line #"<<line<<"#"<<toReturn.isDoubleIndex<<endl;

	if(isFirstLine){
	    if(line[0] == '#'){
		unsigned int i=0;
		int numberOfFields=0;
		bool inWS=true;
		while(i<line.length()){			
		    if( isspace(line[i])){			    
			inWS=true;
		    }else{
			if(inWS){
			    numberOfFields++;
			}
			inWS=false;			    
		    }
		    i++;
		}
		    
		if(numberOfFields==2){ 
		    toReturn.isDoubleIndex=false; 
		}else{
		    if(numberOfFields==3 || numberOfFields==5){
			toReturn.isDoubleIndex=true; 
		    }else{
			cerr << "Must have 2, 3 or 5 fields"<<endl;
			exit(1);
		    }			
		}

	    }else{
		cerr << "First line must begin with #"<<endl;
		exit(1);
	    }
	    isFirstLine=false;
	}else{
	    int i=0;
	    int fieldIndex=0;
	    bool inWS=false;
	    int lastOneNW=0;
	    string foundName;

	    while(i<int(line.length())){		
		    
		if( isspace(line[i]) && i==0){
		    cerr<<line<<endl;
		    cerr << "First character cannot be a space"<<endl;
		    exit(1);
		}
	    
		if( isspace(line[i]) ){			    
		    if(!inWS){ //found a field

			//first field, first index
			if(fieldIndex==0){
			    toReturn.indices1.push_back(toUpperCase(line.substr(lastOneNW,i-lastOneNW)));

			    if(toReturn.mlindex1 < (i-lastOneNW)){
				toReturn.mlindex1 =(i-lastOneNW);
			    }

			}else{
			    //second field, either name of single ind or second index
			    if(fieldIndex==1){
				if(toReturn.isDoubleIndex){
				    toReturn.indices2.push_back(toUpperCase(line.substr(lastOneNW,i-lastOneNW)));
				    if(toReturn.mlindex2 < (i-lastOneNW)){
					toReturn.mlindex2 =(i-lastOneNW);
				    }
				}else{
				    foundName=line.substr(lastOneNW,i-lastOneNW);
				    //duplicated names ?					
				    if(toReturn.namesMap.find(  foundName  ) !=  toReturn.namesMap.end()){
					cerr<<"Warning: The sequence name is duplicated "<<foundName<<endl;
					//exit(1);
				    }else{
					toReturn.namesMap[ foundName ] = ""; 
				    }

				    toReturn.names.push_back( foundName );

				}
			    }else if(fieldIndex==2){
				//sequence name when two indices
				if(toReturn.isDoubleIndex){
				    //duplicated names
				    foundName=line.substr(lastOneNW,i-lastOneNW);
					
				    if(toReturn.namesMap.find(  foundName  ) !=  toReturn.namesMap.end()){
					cerr<<"Warning: The sequence name is duplicated "<<foundName<<endl;
					//exit(1);
				    }else{
					toReturn.namesMap[ foundName ] = ""; 
				    }

				    toReturn.names.push_back( foundName );
				}else{
				    //it's a comment for single index
				    toReturn.namesMap[ foundName ] +=  line.substr(lastOneNW,i-lastOneNW);
				    // cerr<<"Single index file cannot have 3 fields"<<endl;
				    // exit(1);
				}
			    }else{
				//it's a comment again
				toReturn.namesMap[ foundName ] +=  line.substr(lastOneNW,i-lastOneNW);
			    }

			}
			fieldIndex++;
		    }
		    inWS=true;		    
			
		}else{
		    if(inWS)
			lastOneNW=i;
		    inWS=false;			    
		}
		i++;		
	    } //ending while(i<line.length()){		
	}  // ending else firstline

    }  // ending while myFile.good() ){





    //checking for size
    // cout<<toReturn.indices1.size()<<endl;
    // cout<<toReturn.indices2.size()<<endl;
    // cout<<toReturn.names.size()<<endl;
    if(toReturn.isDoubleIndex)
	if((toReturn.indices1.size() != toReturn.indices2.size()) ){
	    cerr << "Size of the fields inconsistent "<<filename<<endl;
	    exit(1);
	}


    if(toReturn.indices1.size() != toReturn.names.size() ){
	cerr << "Size of the fields inconsistent "<<filename<<endl;
	exit(1);
    }


    //checking for valid dna    
    for(unsigned int i=0;i<toReturn.indices1.size();i++){
	if(!isValidDNA(toReturn.indices1[i])){
	    cerr << "Index " << toReturn.indices1[i] <<" is not a valid DNA sequence"<<endl;
	    exit(1);
	}
	if(toReturn.isDoubleIndex)
	    if(!isValidDNA(toReturn.indices2[i])){
		cerr << "Index " << toReturn.indices2[i] <<" is not a valid DNA sequence"<<endl;
		exit(1);
	    }
    }
    return toReturn;
}

void deallocate() {
    // cout<<"deallocate"<<endl;
    delete indTrie1;
    delete indTrie2;
}

map<string,string>  readIndexFile(string filename,int mismatchesTrie,bool _shiftByOne){
    values = intern_readIndex(filename);
    shiftByOne=_shiftByOne;
    
    indTrie1 = new PrefixTree<int>  ();    
    if(values.isDoubleIndex){
	indTrie2 = new PrefixTree<int>  ();
    }
    // PrefixTree<int> pairsOfIndex;
    for(int i=0;i<int(values.names.size());i++){

	indTrie1->insertIntoTree( values.indices1[i].c_str() , i);
	if(values.isDoubleIndex){
	    indTrie2->insertIntoTree( values.indices2[i].c_str() , i);
	    //pairsOfIndex->insertIntoTree( values.indices1[i].c_str()+values.indices2[i].c_str() , i);	    
	}	
    }
    
    
    
    //detect conflicts ?
    cerr<<"Conflicts for index1:"<<endl;
    for(int i=0;i<int(values.names.size());i++){
	vector< int > matchesind ;
	indTrie1->searchMismatch(values.indices1[i].c_str(),&matchesind,mismatchesTrie);
	string toPrint="";

	for(unsigned int j=0;j<matchesind.size();j++){
	    if(matchesind[j] != i)
		toPrint+=values.names[ matchesind[j]  ]+" ";
	}

	if(toPrint.length() > 0){
	    if(values.names[i] != toPrint)
		cerr<<values.indices1[i]<<" from "<<values.names[i]<<" causes a conflict with "<<toPrint<<""<<endl;
	}
    }

    if(values.isDoubleIndex){
	cerr<<"Conflicts for index2:"<<endl;
	for(int i=0;i<int(values.names.size());i++){

	    vector< int > matchesind;
	    indTrie2->searchMismatch(values.indices2[i].c_str(),&matchesind,mismatchesTrie);

	    string toPrint="";

	    for(unsigned int j=0;j<matchesind.size();j++){
		if(matchesind[j] != i)
		    toPrint+=values.names[ matchesind[j]  ]+" ";
	    }

	    if(toPrint.length() > 0){
		cerr<<values.indices1[i]<<" from "<<values.names[i]<<" causes a conflict with "<<toPrint<<endl;
	    }
	}

	cerr<<"Conflicts for pairs:"<<endl;
	for(int i=0;i<int(values.names.size());i++){
	    vector< int > matchesind1;
	    vector< int > matchesind2;

	    indTrie1->searchMismatch(values.indices1[i].c_str(),&matchesind1,mismatchesTrie);	    
	    indTrie2->searchMismatch(values.indices2[i].c_str(),&matchesind2,mismatchesTrie);

	    string toPrint="";

	    for(unsigned int j=0;j<matchesind1.size();j++){
		//found in second as well
		if( find(matchesind2.begin(),
			 matchesind2.end(),
			 matchesind1[j]) != matchesind2.end() 
		    &&
		    matchesind1[j] != i)
		    toPrint+=values.names[ matchesind1[j]  ]+" ";
	    }

	    if(toPrint.length() > 0){
		cerr<<values.indices1[i]+"#"+values.indices2[i]<<" from "<<values.names[i]<<" causes a conflict with "<<toPrint<<endl;
	    }
	}
	
    }

    return values.namesMap ;
}

// likelihood for an index read given the correct index sequence
// this is effectively the sum of qualities of mismatching bases; a
// negative number is our world.
inline double computeLike(const string & indexRef,const string & indexRead,const vector<int> * quals){
#ifdef DEBUG2
    cerr<<"computeLike() "<<indexRef<<"\t"<<indexRead<<endl;
#endif

    double toReturn=0.0;
    for(unsigned int i=0;i<min(indexRef.length(),indexRead.length());i++){

	
	if( indexRef[i] == indexRead[i] ){    
	    toReturn+=likeMatch[    (*quals)[i] ]; 
#ifdef DEBUG2
	    cerr<<"i="<<i<<"\tmatch="<<indexRef[i]<<"\t"<<indexRead[i]<<"\t"<<(*quals)[i]<<"\t"<<likeMatch[    (*quals)[i] ]<<endl;
#endif

	}else{
	    toReturn+=likeMismatch[ (*quals)[i] ]; 
#ifdef DEBUG2
	    cerr<<"i="<<i<<"\tmismatch="<<indexRef[i]<<"\t"<<indexRead[i]<<"\t"<<(*quals)[i]<<"\t"<<likeMismatch[    (*quals)[i] ]<<endl;
#endif

	}

#ifdef DEBUG
	cerr<<"i="<<i<<"\t"<<indexRef[i]<<"\t"<<indexRead[i]<<"\t"<<toReturn<<endl;
#endif

	
    }

    return toReturn;
}



//computes mismatches between index ref and index from the read
inline int computeMM(const string & indexRef,const string & indexRead){
    int toReturn=0;
    for(unsigned int i=0;i<min(indexRef.length(),indexRead.length());i++){
	
	if( indexRef[i] != indexRead[i] )   
	    toReturn++;
	    
    }

    return toReturn;
}

rgAssignment assignReadGroup(string &index1, 
			     string &index1q, 
			     string &index2,
			     string &index2q, 
			     double rgScoreCutoff,
			     double fracConflict,
			     int mismatchesTrie,
			     int qualOffset){
    //BEGIN DEBUG
    // cout<<"DEBUG"<<endl;
    // vector< int > * test1=new vector<int>();
    // indTrie1->searchMismatch("NTGNNNN",test1,2);
    // //    vector< matches<int> > * test2=indTrie2->searchForWordMismatch("CCGGTAC",0);
    // cout<<test1->size()<<endl;
    // //    cout<<test2->size()<<endl;

    // for(unsigned int j=0;j<test1->size();j++){
    // 	// list< int  >::const_iterator iter;		    
    // 	// for (iter  = (*test1)[j].listOfDeflinesOfMatches->begin(); 
    // 	//      iter != (*test1)[j].listOfDeflinesOfMatches->end(); 
    // 	//      iter++){
    // 	//     cout<<"test "<<*iter<<"\t"<<values.indices1[*iter]<<"\t"<<values.names[*iter]<<endl;
    // 	// }
    // 	cout<<"test "<< (*test1)[j] <<"\t"<<values.indices1[  (*test1)[j] ]<<"\t"<<values.names[ (*test1)[j] ]<<endl;
    // }
    // delete test1;
    //exit(1);
    //END DEBUG

    rgAssignment toReturn;

    vector<int> quals1;
    vector<int> quals2;

    //Find matches using the prefix trie

    vector< int > matchesind1, matchesind2;
    indTrie1->searchMismatch(      index1.c_str(),&matchesind1,mismatchesTrie);

    
    if(shiftByOne){       
	indTrie1->searchMismatch(     ("N"+index1.substr(0, index1.size() -1)     ).c_str(),&matchesind1,mismatchesTrie);
	indTrie1->searchMismatch(         (index1.substr(1, index1.size() -1)+"N" ).c_str(),&matchesind1,mismatchesTrie);
    }

    for(unsigned int i=0;i<index1q.length();i++){
	char tempCtocov   = char(index1q[i]);
	int tempIntTopush = (int( tempCtocov )-qualOffset);
	quals1.push_back(     max(tempIntTopush ,2)  )   ; //since qual scores less than 2 do not make sense
    }


    if(!index2.empty()){
	if(!values.isDoubleIndex){
	    if(!warnedForDouble){
		cerr<<"WARNING : THIS RUN IS DOUBLE INDEXED YET YOUR INDEX IS SINGLE INDEXED"<<endl;
		warnedForDouble=true;
	    }
	}else{
	    indTrie2->searchMismatch(  index2.c_str(),&matchesind2,mismatchesTrie);

	    if(shiftByOne){       
		indTrie2->searchMismatch( ("N"+index2.substr(0, index2.size() -1)     ).c_str(),&matchesind2,mismatchesTrie);
		indTrie2->searchMismatch(     (index2.substr(1, index2.size() -1)+"N" ).c_str(),&matchesind2,mismatchesTrie);
	    }
	
	    for(unsigned int i=0;i<index2q.length();i++){
		char tempCtocov   = char(index2q[i]);
		int tempIntTopush =(int(tempCtocov)-qualOffset);
		quals2.push_back( max( tempIntTopush ,2)  )   ;	//since qual scores less than 2 do not make sense
	    }
	}
    }


    set<int> foundIndices; //set of indices found

    double like1=0.0;
    double like2=0.0;

    vector< pair<int,double>  > sortedLikelihoodAll; //likelihood for both indices
    vector< pair<int,double>  > sortedLikelihood1;   //likelihood for index 1
    vector< pair<int,double>  > sortedLikelihood2;   //likelihood for index 1

    //putting all the index # into the set
    for(unsigned int j=0;j<matchesind1.size();j++)
        foundIndices.insert ( matchesind1[j]  );   

    for(unsigned int j=0;j<matchesind2.size();j++)
        foundIndices.insert ( matchesind2[j]  );   


    //for every # in the set, compute likelihood    
    for (set<int>::iterator sit=foundIndices.begin(); sit!=foundIndices.end(); sit++){

	like1      =  computeLike(index1,values.indices1[ *sit ],&quals1);
	
	if(shiftByOne){ //check if shifting by one improves the likelihood
	    like1      =  max(like1,
			      max(
				  computeLike( ("N"+index1.substr(0, index1.size() -1)     ), values.indices1[ *sit ],&quals1),
				  computeLike( (    index1.substr(1, index1.size() -1)+"N" ), values.indices1[ *sit ],&quals1)
				  )			      
			      );
	}

	sortedLikelihood1.push_back(     make_pair (*sit,like1) );

	if(!index2.empty()  && 	values.isDoubleIndex ){

	    like2  =  computeLike(index2,values.indices2[ *sit ],&quals2);

	    if(shiftByOne){ //check if shifting by one improves the likelihood
		like2      =  max(like2,
				  max(
				      computeLike( ("N"+index2.substr(0, index2.size() -1)     ), values.indices2[ *sit ], &quals2),
				      computeLike( (    index2.substr(1, index2.size() -1)+"N" ), values.indices2[ *sit ], &quals2)
				      )
				  );
	    }

	    sortedLikelihood2.push_back( make_pair (*sit,like2) );
	}else{
	    like2=0.0;
	}

	sortedLikelihoodAll.push_back( make_pair (*sit,like1+like2) ) ;      
    }

    //if nothing was found, unknown
    if(foundIndices.empty() ){
	toReturn.predictedGroup.clear();
	return toReturn;
    }
    //At this point, we will return a RG 

    //sorting by likelihood
    sort (sortedLikelihood1.begin(),   sortedLikelihood1.end(),   comparePair()); 
    sort (sortedLikelihood2.begin(),   sortedLikelihood2.end(),   comparePair()); 
    sort (sortedLikelihoodAll.begin(), sortedLikelihoodAll.end(), comparePair()); 

#ifdef DEBUG
    cerr<<endl;
    cerr<<"first:"<<endl;
    for(unsigned int j=0;j<sortedLikelihood1.size();j++){
    	cerr<< values.names[ sortedLikelihood1[j].first ]<<"\t"<< sortedLikelihood1[j].second<<endl;
    }
    cerr<<"second:"<<endl;
    for(unsigned int j=0;j<sortedLikelihood2.size();j++){
    	cerr<< values.names[ sortedLikelihood2[j].first ]<<"\t"<< sortedLikelihood2[j].second<<endl;
    }
    cerr<<"all:"<<endl;
    for(unsigned int j=0;j<sortedLikelihoodAll.size();j++){
    	cerr<< values.names[ sortedLikelihoodAll[j].first ]<<"\t"<< sortedLikelihoodAll[j].second<<"\t"<<pow(10.0,sortedLikelihoodAll[j].second)<<endl;
    }
    cerr<<endl;
    //exit(1);
#endif

    // DETECT WRONGS
    // Look for a wrong index.  Suppose the top indices do not match up,
    // then those give the likelihood for being wrong, which we compare
    // to that of the top correct pair.
    //
    // Suppose they do, then we have to find the highest scoring wrong
    // pair.  That's either the top first index with the runner-up
    // second index, or vice versa.  Could actually be both to some
    // extent, so we add them.

    if( (sortedLikelihood1.size()>1) && 
	(sortedLikelihood2.size()>1)   ){
        if(sortedLikelihood1[0].first != sortedLikelihood2[0].first) { // mismatch, we found the wrong pair
           
	    //find likelihood of p5 (in sortedLikelihood2)  for p7 top hit (sortedLikelihood1[0].first)
	    int  indexP5LikeForP7Top=-1;
	    for(unsigned int j=0;j<sortedLikelihood2.size();j++){
		if(sortedLikelihood2[j].first == sortedLikelihood1[0].first){
		    indexP5LikeForP7Top = int(j);
		    break;
		}
	    }

	    //find likelihood of p7 (in sortedLikelihood1)  for p5 top hit (sortedLikelihood2[0].first)
	    int  indexP7LikeForP5Top=-1;
	    for(unsigned int j=0;j<sortedLikelihood1.size();j++){
		if(sortedLikelihood1[j].first == sortedLikelihood2[0].first){
		    indexP7LikeForP5Top = int(j);
		    break;
		}
	    }

	    if(indexP7LikeForP5Top == -1 ||
	       indexP5LikeForP7Top == -1 ){
		cerr<<"Internal error, unable to trace back the likelihood for one of the index pairs"<<endl;
		exit(1);
	    }
	       
	  
	    
            toReturn.topWrongToTopCorrect =    1.0*(  sortedLikelihood1[                  0].second + sortedLikelihood2[                   0 ].second )
		                             + 0.3 
		                             - oplus( sortedLikelihood1[                  0].second + sortedLikelihood2[ indexP5LikeForP7Top ].second 
						      , 
						      sortedLikelihood1[indexP7LikeForP5Top].second + sortedLikelihood2[                   0 ].second ) ;
	    

#ifdef DEBUGMISPAIR
	    cout<<"conf\t"<<toReturn.topWrongToTopCorrect<<"\t"<<-10*toReturn.topWrongToTopCorrect<<endl;
#endif
        }else { //stem from the same sample
            // we compare one correct pair to two potentially wrong
            // ones; add 0.3 to make it fair
            toReturn.topWrongToTopCorrect = - 1.0* sortedLikelihoodAll[0].second 
                                 	    - 0.3 
		                            + oplus( sortedLikelihood1[0].second + sortedLikelihood2[1].second 
						     , 
						     sortedLikelihood1[1].second + sortedLikelihood2[0].second ) ;

#ifdef DEBUGMISPAIR
	    cout<<"same\t"<<toReturn.topWrongToTopCorrect<<"\t"<<-10*toReturn.topWrongToTopCorrect<<endl;
#endif
        }
    }else{
	toReturn.topWrongToTopCorrect = -1.0*1000000  ; // +infinity for practical purposes
    }

    // DETECT CONFLICTS
    // Checking likelihood of inferior hits; if the ratio is too low,
    // it's a conflict.  Checking the second best only is already a
    // useable approximation, adding all is more appropriate (and a bit
    // more expensive).
    double probRG    = sortedLikelihoodAll[0].second;
    toReturn.predictedGroup = values.names[ sortedLikelihoodAll[0].first ];
    toReturn.logLikelihoodScore = probRG;

    if(sortedLikelihoodAll.size() > 1){
	double probRG2nd    = sortedLikelihoodAll[1].second;
        for( size_t i = 2 ; i != sortedLikelihoodAll.size() ; ++i )
            probRG2nd = oplus( probRG2nd, sortedLikelihoodAll[i].second ) ; //oplus= log10( pow(10,x)+pow(10,y) )

	toReturn.logRatioTopToSecond = probRG2nd - oplus(probRG,probRG2nd) ;
	double temporaryD=-10.0*toReturn.logRatioTopToSecond;
	if(flag_ratioValues ) // && toReturn.logRatioTopToSecond > -5) //to avoid very small values
	    ratioValues->write( (char *)&(temporaryD), sizeof(toReturn.logRatioTopToSecond));
    }else{
	toReturn.logRatioTopToSecond = 1;
    }

    if(flag_rgqual){
	double temporaryD=-10.0*toReturn.logLikelihoodScore;
	rgqual->write( (char *)&(temporaryD), sizeof(toReturn.logLikelihoodScore));
    }


#ifdef DEBUG3    

    toReturn.numberOfMismatches      = computeMM(index1,values.indices1[ sortedLikelihoodAll[0].first ]);
    if(!index2.empty())
	toReturn.numberOfMismatches += computeMM(index2,values.indices2[ sortedLikelihoodAll[0].first ]);
    
    //cerr<<endl;
    // if(toReturn.numberOfMismatches!=0){
    // 	cerr<<toReturn.logLikelihoodScore<<endl;
    // cerr<<toReturn.topWrongToTopCorrect<<endl;
    // cerr<<toReturn.logRatioTopToSecond<<endl;
    //cerr<<"MM"<<toReturn.numberOfMismatches<<endl;
    //}
    //cerr<<toReturn.numberOfMismatches<<"\t"<<toReturn.logLikelihoodScore<<endl;
#endif
    
    return toReturn;
}
