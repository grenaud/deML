#ifndef RGAssign_h
#define RGAssign_h

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <algorithm>
#include <vector>
#include <math.h>
#include <sys/time.h>
//#include <ctype.h>
#include <string>
#include <map>
#include <set>

#include "PrefixTree.h"

using namespace std;

typedef struct {
    map<string,string>  namesMap;

    vector<string>  names;
    vector<string>  indices1;
    vector<string>  indices2;    
    bool isDoubleIndex;
    int mlindex1;
    int mlindex2;   
} indexData;

typedef struct{
    string predictedGroup;    
    bool conflict;
    bool wrong;
    bool unknown;

    double logRatioTopToSecond;            // ~~ conflict
    double logLikelihoodScore;             // ~~ unknown
    double topWrongToTopCorrect;        // ~~ wrong
    int numberOfMismatches;          //total # of mismatches
} rgAssignment;

string toUpperCase(string toCheck);
bool isValidDNA(string tocheck);
void setFileForRatio(ofstream * streamFile);
void setFileForRGQual(ofstream * streamFile);

indexData intern_readIndex(string filename);
map<string,string>  readIndexFile(string filename,int mismatchesTrie,bool _shiftByOne);
rgAssignment assignReadGroup(string  & index1,string & index1q,string & index2,string & index2q,double rgScoreCutoff,double fracConflict,int mismatchesTrie);
void deallocate();

/* static double cutoffRatioLogLike=0.5; */


#endif
