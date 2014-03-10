/*
 * PrefixTree
 * Date: Nov-09-2012 
 * Author : Gabriel Renaud gabriel.reno [at sign here ] gmail.com
 *
 */

#ifdef __PREFIXTREE_H__


//for DNA
#define SIZEALPHA 4
// #define DEBUG

static string  ALPHA = "ACGT";

// int longestCommonSubstring(const string &s1,const string & s2){
//     int i=0;
//     int ind1=int(s1.length()-1);
//     int ind2=int(s2.length()-1);
//     while(i<=min(ind1,ind2) ){
// 	if(s1[i] != s1[i] )
// 	    return i;
// 	i++;
//     }
//     return 0;
// }



template <typename typeOfVar> 
Node<typeOfVar>::Node(){

    children= new Node<typeOfVar> * [SIZEALPHA];
    for(int i =0;i<SIZEALPHA;i++){
	children[i] = NULL;
    }
    myIdentifiers = new vector< typeOfVar > ();

    characters = new string("") ;

    leaf     =false;
    endOfWord=false;
}



template <typename typeOfVar> 
Node<typeOfVar>::Node(const char * toinsert,const vector<typeOfVar > *identifier){
    children= new  Node<typeOfVar>  * [SIZEALPHA];
    for(int i =0;i<SIZEALPHA;i++){
	children[i] = NULL;
    }
    myIdentifiers = new vector< typeOfVar > (*identifier);

    if(toinsert[0] == '\0'){ //empty char
	characters = new string("") ;
	endOfWord=true; //The path represents a string
    }else{
	characters = new string(toinsert) ;
	endOfWord=false; //The path alone does not represent a string
    }
    leaf     =true;

}



template <typename typeOfVar> 
Node<typeOfVar>::~Node(){

    for(int i =0;i<SIZEALPHA;i++){
	if(children[i] != NULL)
	    delete children[i];
    }
    delete [] children;
    delete characters;
    delete myIdentifiers;
}

template <typename typeOfVar> 
bool Node<typeOfVar>::isLeaf(){
    for(int i =0;i<SIZEALPHA;i++)
	if(children[i] != NULL)
	    return false;
    
    return true;
}




/*  This subroutine is designed to create a set of subtree
 *  when trying to add a new record to an existing leaf.
 *  The leaf is no longer a leaf and the information is contains
 *  is copied into (toinsert2,identifier2) and the previous characters
 *  and identifiers have been emptied (by the calling function, in our case insertIntoNode())
 *  The new information from the record to add is stored in (toinsert1,identifier1). 
 *  This subroutine can create long branches in the subtree, one node for each common character.
 *
 *  There are 4 cases:
 *  1) both strings are empty, this case should not happend because the calling function insertIntoNode() should have handled this case
 *  2) One of them is empty, is will be the parent node, the other will be the child
 *  None are empty and :
 *  3) They start with the same character, the current node becomes a parent to non-leaf child nodes
 *  4) They do not start with the same character, the current node becomes a parent to leaf child nodes
 *  
 */
template <typename typeOfVar> 
void Node<typeOfVar>::createLeaves(const char * toinsert1,const vector<typeOfVar> * identifier1,
				   const char * toinsert2,const vector<typeOfVar> * identifier2){

    //case 1
    //both strings are empty, we just need to hold the identifiers
    if(toinsert1[0] == '\0' &&
       toinsert2[0] == '\0' ){
	
	cerr<<"ERROR: Case 1 in createLeaves(), exiting"<<endl;
	exit(1);
	// leaf     =true; //since it was called on a leaf
	// endOfWord=true; //for insert1 and insert2
	
	// for(int id=0;id<int(identifier1->size());id++)
	//     myIdentifiers->push_back( identifier1->at(id) );	    	

	// for(int id=0;id<int(identifier2->size());id++)
	//     myIdentifiers->push_back( identifier2->at(id) );	    	

	// return;
    }


    //case 2
    if(toinsert1[0] == '\0' ){ //insert 2 not empty, will become child but will hold insert1
	int indexBp=bp2Index(toinsert2[0]);
	leaf     =false;
	endOfWord=true; //for insert1
	
	for(int id=0;id<int(identifier1->size());id++)
	    myIdentifiers->push_back( identifier1->at(id) );	    	

	children[ indexBp ] = new Node< typeOfVar >( toinsert2,identifier2);
	return;
    }

    if(toinsert2[0] == '\0' ){ //insert 1 not empty, will become child but will hold insert2
	int indexBp=bp2Index(toinsert1[0]);
	leaf     =false;
	endOfWord=true; //for insert2

	for(int id=0;id<int(identifier2->size());id++)
	    myIdentifiers->push_back( identifier2->at(id) );

	children[ indexBp ] = new Node< typeOfVar >( toinsert1,identifier1);
	return;
    }


    int indexBp1=bp2Index(toinsert1[0]);
    int indexBp2=bp2Index(toinsert2[0]);

    //case 3
    if(indexBp1 == indexBp2){//share a common character
	//we know the child will not be a leaf
	children[ indexBp1 ] = new Node< typeOfVar >();
	children[ indexBp1 ]->createLeaves(toinsert1+1,identifier1,
					  toinsert2+1,identifier2);
	return ;
    //case 4
    }else{//will split with two kids
	leaf      = false;  
	endOfWord = false; //since we do not hold a terminal word
	children[ indexBp1 ] = new Node< typeOfVar >(toinsert1+1,identifier1);
	children[ indexBp2 ] = new Node< typeOfVar >(toinsert2+1,identifier2);
	return ;

    }
    

}


/*  This subroutine inserts the string into the node or
 *  recursively calls it on one of children.
 *
 *  There are the following 5 cases:
 *  1) We have reached the end of the characters in toinsert, just insert its info into the current node
 *  The current node is a leaf and: 
 *     2) The strings are identical, insert its info into the current node
 *     3) The strings are differ (most complex case), then we have 
 *        - to store the current information 
 *        - delete our "characters" and "myIdentifiers"
 *        - Call createLeaves() using the current information and the new one 
 *          This will create the appropriate subtree
 *        - The current node will become a parent
 *  The current node is a not a leaf and: 
  *    4) It has no child with the same character, create it as a leaf
 *     5) It has already a child with the same character as the string, recursively call the procedure on it
 * 
 */
template <typename typeOfVar> 
void Node<typeOfVar>::insertIntoNode(const char * toinsert,const vector<typeOfVar > * identifier){


    //case 1
    //if the string is empty, need to insert it here
    if(toinsert[0] == '\0'){
	endOfWord=true; //represents an end of word 

	for(int id=0;id<int(identifier->size());id++)
	    myIdentifiers->push_back( identifier->at(id) );

    }else{ //char is not \0
	int indexBp=bp2Index(toinsert[0]);
	
	
	//The current node is a leaf, no children
	if(leaf){
#ifdef DEBUG
	    cout<<"Trying to insert into leaf #"<<indexBp<<" with seq = " <<(toinsert)<<endl;
#endif
	    //case 2
	    if( characters->compare(toinsert) == 0){//the new string to insert is identical to ours
#ifdef DEBUG
		cout<<"Duplicating entries with seq = " <<(toinsert)<<endl;
#endif
		//cerr<<"WARNING: Duplicated entries in prefix tree "<<endl;
		for(int id=0;id<int(identifier->size());id++)
		    myIdentifiers->push_back( identifier->at(id) );	    	
		return ;
		//case 3
	    }else{
		//else
		//if it's a leaf with an empty string
		if(characters->empty()){
#ifdef DEBUG
		    cout<<"Inserting to an empty record a leaf #"<<indexBp<<"  " <<(toinsert)<<endl;
#endif
		    children[ indexBp ] = new Node( (toinsert+1),identifier);
		    leaf=false;
		    return ;
		}else{
#ifdef DEBUG
		    cout<<"Inserting to an non-empty record ("<< *characters<<") a leaf #"<<indexBp<<"  " <<(toinsert)<<endl;
#endif
		    //need to descend to create the common branch until we
		    //reach an unequal char
		    //This node will become a non-leaf node
		    leaf     =false;      //not a leaf anymore
		    endOfWord=false;      //maybe not an end of word anymore, createLeaves() will initialize both properly

		    //copying contents
		    string              * strTemp = new string(*characters);
		    vector< typeOfVar > * idTemp  = new vector< typeOfVar > ();		    
		    idTemp->assign(myIdentifiers->begin(),myIdentifiers->end());

		    //removing contents
		    characters->erase();                                                //need to erase the characters
		    myIdentifiers->erase(myIdentifiers->begin(),myIdentifiers->end() ); //need to erase the identifiers

		    this->createLeaves(toinsert,            identifier,    //new data
				       strTemp->c_str(),    idTemp);       //previous data
		    delete(strTemp);
		    delete(idTemp);

		}
	    }
	}else{ //not a leaf
	    //case 4
	    if(children[ indexBp ] == NULL){ //child node is not created, create it, no longer a leaf
#ifdef DEBUG
		cout<<"Creating child# "<<indexBp<<" with seq = " <<(toinsert+1)<<endl;
#endif
		children[ indexBp ] = new Node< typeOfVar >( (toinsert+1),identifier);
		leaf=false;
		return ;
		//case 5
	    }else{//it already exists, recurse on the child
#ifdef DEBUG
		cout<<"Inserting into child# "<<indexBp<<"  with seq = " <<(toinsert+1)<<endl;
#endif
		children[ indexBp ] -> insertIntoNode( (toinsert+1),identifier);
		return ;
	    }
	}

    }
}

template <typename typeOfVar> 
void Node<typeOfVar>::addFromNode(vector<typeOfVar > * foundIds) const{
    foundIds->insert(foundIds->end(),myIdentifiers->begin(),myIdentifiers->end());
}

/*  This subroutine returns all the ids in the subtree
 *
 */
template <typename typeOfVar> 
void Node<typeOfVar>::addFromSubtree(vector<typeOfVar > * foundIds) const{
    addFromNode(foundIds);
    for(int i =0;i<SIZEALPHA;i++){
    	if(children[i] != NULL){
	    children[i]->addFromSubtree(foundIds);
	}
    }
}




template <typename typeOfVar> 
PrefixTree<typeOfVar>::PrefixTree(){

    root = new Node<typeOfVar> * [SIZEALPHA];
    for(int i =0;i<SIZEALPHA;i++){
    	root[i] = NULL;
    }

    test=true;
}


template <typename typeOfVar> 
PrefixTree<typeOfVar>::~PrefixTree(){
    for(int i =0;i<SIZEALPHA;i++){
    	if(root[i] != NULL)
	    delete root[i];
    }
    delete []  root;
}

template <typename typeOfVar> 
void PrefixTree<typeOfVar>::insertIntoTree(const char * toinsert,typeOfVar identifier){
    if(toinsert[0] == '\0' ){ 
	cerr<<"Cannot add an empty string"<<endl;
	exit(1);
    }

    int indexBp=bp2Index(toinsert[0]);
    vector<typeOfVar> * tempVec=new vector<typeOfVar>();
    tempVec->push_back(identifier);
    
    if(root[ indexBp ] == NULL){ //child node is not created, create it, no longer a leaf
#ifdef DEBUG
	cout<<"Creating root node# "<<indexBp<<"  with seq = " <<(toinsert+1)<<endl;
#endif
	root[ indexBp ] = new Node< typeOfVar >( (toinsert+1),tempVec);
    }else{//it already exists, recurse on the child
#ifdef DEBUG
	cout<<"Inserting into root node# "<<indexBp<<"  with seq = " <<(toinsert+1)<<endl;
#endif
	root[ indexBp ] -> insertIntoNode( (toinsert+1),tempVec);
    }
    delete tempVec;
}


/*  This subroutine prints the tree to stdout
 *  
 */
template <typename typeOfVar> 
void PrefixTree<typeOfVar>::printTree(){
    for(int i =0;i<SIZEALPHA;i++){
    	if(root[i] != NULL){
	    internPrintTree(string(1,ALPHA[i]),root[i]);
	}}
}


/*  This subroutine searches the tree with mismatches
 *  
 */
template <typename typeOfVar> 
void PrefixTree<typeOfVar>::internPrintTree(string path,Node<typeOfVar> * n){

    for(int i =0;i<SIZEALPHA;i++){
    	if(n->children[i] != NULL){
	    internPrintTree(path+ALPHA[i],n->children[i]);
	}}

    if(n->endOfWord){
	cout<<path<<*(n->characters)<<"\t"<<vectorToString<typeOfVar>((*n->myIdentifiers))<<endl;
    }
}





template <typename typeOfVar> 
void PrefixTree<typeOfVar>::internSearchMismatch(const char * target,Node<typeOfVar > * currentNode,vector<typeOfVar > * foundIds,int nMMfound,int nMMTolerated){
#ifdef DEBUG
    cout<<"internSearchMismatch target="<<target<< " node "<<*currentNode<<endl;
#endif
    if(nMMfound>nMMTolerated)
	return;

    //match and also matches everything in subtrees
    if(target[0] == '\0' ){ 
	currentNode->addFromSubtree(foundIds);
	return;
    }

#ifdef DEBUG
    cout<<"Found ids1 "<<vectorToString(*foundIds)<<endl;
#endif

    //signals end of character
    if(currentNode->endOfWord && !currentNode->leaf){ //if it's a leaf, it will get added later
	currentNode->addFromNode(foundIds);
    }
#ifdef DEBUG
    cout<<"Found ids2 "<<vectorToString(*foundIds)<<endl;
#endif    
    //
    if(currentNode->leaf){
#ifdef DEBUG
	cout<<"LEAF  "<<(*currentNode) <<endl;
#endif
	if( matchesPrefixes( (currentNode->characters)->c_str(),target,(nMMTolerated-nMMfound)) ){
	    currentNode->addFromNode(foundIds);
	}
    }else{
#ifdef DEBUG
	cout<<"NLEAF "<<(*currentNode) <<endl;
#endif
	int indexBp=bpTarget2Index(target[0]);
#ifdef DEBUG
	cout<<"id "<<indexBp<<endl;
#endif
	for(int i =0;i<SIZEALPHA;i++){
#ifdef DEBUG
	    cout<<"i "<<i<<endl;
#endif
	    if(indexBp == 4 ) {//tolerate 'N' but discount one mismatch
		if(currentNode->children[i] != NULL)
		    internSearchMismatch(         target+1, currentNode->children[i], foundIds, nMMfound+1 ,nMMTolerated);
	    }else
		if(currentNode->children[i] != NULL){
		    if(indexBp  == i){
			internSearchMismatch( target+1, currentNode->children[i],     foundIds, nMMfound+0 ,nMMTolerated);
		    }else{
			internSearchMismatch( target+1, currentNode->children[i],     foundIds, nMMfound+1 ,nMMTolerated);
		    }
		}
	}
    
    }

}

template <typename typeOfVar> 
void PrefixTree<typeOfVar>::searchMismatch(const char * target,vector<typeOfVar > * foundIds,int numberMisMatches){
    string test (target);
#ifdef DEBUG
    cout<<"test "<<test<<endl;
    cout<<"searchMismatch "<<(target)<<endl;
#endif    

    if(target[0] == '\0' ){ 
	cerr<<"Cannot search for an empty string"<<endl;
	exit(1);
    }

    if(numberMisMatches<0){
	cerr<<"Cannot search with negative mismatches"<<endl;
	exit(1);
    }


    int indexBp=bpTarget2Index(target[0]);

    for(int i =0;i<SIZEALPHA;i++){
	if(indexBp == 4 ) {//tolerate 'N' but discount one mismatch
	    if(numberMisMatches>0)
		if(root[i] != NULL){
#ifdef DEBUG
		    cout<<"from root N  "<<("ACGT"[i])<<"\t"<< *(root[i])<<endl;
#endif    

		    internSearchMismatch(      target+1, root[i], foundIds,  1, numberMisMatches);
		}
	}else
	    if(root[i] != NULL){
		if(indexBp  == i){
#ifdef DEBUG
		    cout<<"from root MA "<<("ACGT"[i])<<"\t"<< *(root[i])<<endl;
#endif    

		    internSearchMismatch(      target+1, root[i], foundIds,  0, numberMisMatches);
		}else{
		    if(numberMisMatches>0){
#ifdef DEBUG
			cout<<"from root MM  "<<("ACGT"[i])<<"\t"<< *(root[i])<<endl;
#endif    

			internSearchMismatch(  target+1, root[i], foundIds,  1, numberMisMatches);
		    }
		}
	    }
    }
    
}



#endif
