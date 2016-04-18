#ifndef INF1_H
#define INF1_H
#include <vector>
#include "math.h"
#include "state.h"
#include "rng.h"
#include "node.h"
#include "infextract.h"

namespace fines
{

/**
    @brief Inference algorithm 1: non-mcmc
*/
class Inf1
{
public:
    Inf1(Data *d,State *s,Data *d2=NULL,int datainference=INFDATA_COUNTS,bool v=false,int tmax=1500,int treescale=0);///<Constructor if everything already known
    Inf1(Data *d,State *s,string newick,bool v=false);///<Constructor if tree already known
    Inf1(Data *d,int initpop=-1,Data *d2=NULL,int datainference=INFDATA_COUNTS,int modeltype=1,bool v=false,int tmax=1500);///<Constructor from data only.  Starts with one population or N populations

    bool doBestMerge(int mmax=100,bool stopattop=true,int treemodification=0);///<Calculates and acts on the best merge to do
    bool doBestSplit(bool stopattop=true);///<Calculates and acts on the best split to do
    void mergeHillClimb(std::ostream * fileout,bool stopattop,int treemodification);///<Merges until the best population is reached
    void splitHillClimb(bool stopattop=true);///<Splits until the best population is reached
    void doHyperParameterMaximisation(int numsteps);///<Does greedy hyperparameter moves to maximise the hyperparameters
    
    void linkNodes(int i, int j,double d,bool usepopnodes,bool intree,bool absage=false);///<Links nodes i and j with distance d (using internal popnode representation or not)
    void createSkeletonTree();///<Creates the tree nodes needed on creation

    inline State* getState(){return(state);}///<Gets the state
    void statusOut(long c1, long numints);///< Prints the status bar
    void exportXmlHead(std::ostream * fileout,std::string fs,std::string type,double c,long x);///<Xml start
    void exportXmlTail(std::ostream * fileout);///<Xml end
    void diagonaliseOnM(std::vector<std::vector<double> > *mat,bool symm);///<diagonalise on the given matrix
    void setcertainty(std::vector<std::vector<double> > *mat);///< Evaluates the certainty of each node
    //void diagonaliseStatesOnM(State *in, std::vector<std::vector<double> > *mat,bool symm);///<diagonalise in over the states on M // NOT WORKING!
    bool inOrder (int i,int j);///< returns true if i appears before j in the same population, false if after, error otherwise
    bool inFullOrder (int i,int j); ///< returns true if state a appears before b, false if b before a, error otherwise
    void reorderState(State *in=NULL);///< reorder a state to match the tree
    inline void exportXmlComment(std::ostream * fileout,std::string comment){
	    *fileout<< "<comment>"<<comment<<"</comment>"<<std::endl;
    };///<Print a comment tag
    void printTree(std::ostream * fileout,bool xml=true);///<Prints the newick tree
    ~Inf1();
    int whichPopNode(int i);
    inline std::vector<std::vector<double>* > swapM(std::vector<std::vector<double>* > m, int a,int b){
	std::vector< std::vector<double>* > tmp;
	for(unsigned int c1=0;c1<m.size();c1++) tmp.push_back(m[a]);
	unsigned int ts=tmp.size();
	for(unsigned int c1=0;c1<m.size();c1++) tmp.push_back(m[b]);
	for(unsigned int c1=0;c1<m.size();c1++) m[a]=tmp[c1+ts];
	for(unsigned int c1=0;c1<m.size();c1++) m[b]=tmp[c1];
	return(m);
    }///<Reorders a matrix according to a and b swapping in the tree
    inline std::vector<Node*> getNodes(){return(nodes);}///<gets avector of pointers to the nodes
    inline Node * getRoot(){return(root);}
    inline vector<int> ignoreIgnorables(vector<int> plist){
      int counter=0;
      unsigned int oldsize=plist.size();
      while(counter<(int) plist.size()){
	  if(data->getIgnore(plist[counter])) plist.erase(plist.begin()+counter);
	  else counter++;
      }
      if(plist.size()!=oldsize && oldsize>1) cerr<<"WARNING!!! AA plist="<<plist.size()<<" But ignoring population!"<<endl;
      return(plist);
    }
    void reorderSuper(State *in);// places super individuals at the start of the population list
    inline bool testIgnorePop(int p){
      bool found=false;
      vector<int> plist=state->getIndInPop(p);
      for(unsigned int c1=0;c1<plist.size();c1++){
	if(data->getIgnore(plist[c1])) found=true;//return(true);
      }
      if(found && plist.size()>1) cerr<<"WARNING!!! VV plist="<<plist.size()<<" But ignoring population!"<<endl;
      return(found);
    }///< Returns true if any individuals in a population are super
protected:
    Data *data;
    State *state;
    bool verbose;
    int test_max;///<Maximum number of tests to perform on each split/merge

    std::vector<int> popnodes;///<Nodes for the populations
    Node *root;///<Root node of the genealogy
    std::vector<Node*> nodes;///<Vector of all the nodes in the genealogy
    int assignednodes;///<Index for the last asssigned node
    bool warned;///< a flag to say if we have issued the "too many combinations" warning
    int treescale;///<how we form the tree; 0:	Discard all ordering and likelihood information, 1: keep ordering, 2 keep everything
    double initheight,currentheight;///<If we are keeping track of heights, this is the reference height
};

} // end namespace fines
#endif
