/*
 * PrefixTree
 * Date: Nov-09-2012 
 * Author : Gabriel Renaud gabriel.reno [at sign here] gmail.com
 *
 */
#include <vector>
#include <stdlib.h>

#include "utils.h"

using namespace std;

/* int longestCommonSubstring(string &s1,string s2); */
/* inline bool matchesPrefixes(const char * s1,const char * s2,int maxMM); */

inline bool matchesPrefixes(const char * s1,const char * s2,int maxMM){
    int i      =0;
    int mmFound=0;

    while(true){
	if(s1[i] == '\0')
	    return true;
	if(s2[i] == '\0')
	    return true;
	if(s1[i] !=  s2[i]){
	    mmFound++;
	    if(mmFound>maxMM)
		return false;
	}
	i++;
    }
    
}

inline int bp2Index(char c){
    if(c == 'A'){ return 0;}
    if(c == 'C'){ return 1;}
    if(c == 'G'){ return 2;}
    if(c == 'T'){ return 3;}
    cerr<<"Invalid character "<<c<<" exiting "<<endl;
    exit(1);
}

inline int bpTarget2Index(char c){
    if(c == 'A'){ return 0;}
    if(c == 'C'){ return 1;}
    if(c == 'G'){ return 2;}
    if(c == 'T'){ return 3;}
    if(c == 'N'){ return 4;}

    cerr<<"Invalid character "<<c<<" exiting "<<endl;
    exit(1);
}

//forward declare for friend
template <typename typeOfVar> 
class PrefixTree;


template <typename typeOfVar>
class Node{
 private:

    Node   ** children ;
    bool leaf;
    bool endOfWord;
    string * characters;

    vector< typeOfVar > * myIdentifiers; //vector of identifiers associated with this node
    
 public:
    Node(const char * toinsert,const vector<typeOfVar> * identifier); //for leaf
    Node(); //for non-leaf

    ~Node();
    bool isLeaf();
    void insertIntoNode(const char * toinsert, const vector<typeOfVar> * identifier);
    void createLeaves(  const char * toinsert1,const vector<typeOfVar> * identifier1,
		        const char * toinsert2,const vector<typeOfVar> * identifier2);
    void addFromNode(   vector<typeOfVar > * foundIds) const;
    void addFromSubtree(vector<typeOfVar > * foundIds) const;
    friend class PrefixTree<typeOfVar>;

    friend ostream& operator<<(ostream& str,Node<typeOfVar >  const& nd){
	str<<"Node: char=\""<<*(nd.characters)<<"\" leaf="<<boolStringify(nd.leaf)<<" eow="<<boolStringify(nd.endOfWord) <<" contains=\""<<vectorToString( *(nd.myIdentifiers) )<<"\""<<endl;
	return str;
    }

};


template <typename typeOfVar>
class PrefixTree{
 private:

    Node<typeOfVar>   ** root ;

    bool test;
    void internPrintTree(string path,Node<typeOfVar> * n);
    void internSearchMismatch(const char * target,Node<typeOfVar > * currentNode,vector<typeOfVar > * foundIds,int nMMfound,int nMMTolerated);

 public:
    PrefixTree();
    ~PrefixTree();
    void insertIntoTree(const char * toinsert,typeOfVar identifier);
    void printTree();
    void searchMismatch(const char * target,  vector<typeOfVar > * foundIds,int numberMisMatches);



};



#define __PREFIXTREE_H__
#include "PrefixTree.cpp"
#undef  __PREFIXTREE_H__

