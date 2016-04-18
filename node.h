#ifndef NODE_H
#define NODE_H
#include <string>
#include <vector>
#include <iostream>
#include "data.h"

namespace fines
{

/**
    @brief A node of the tree
*/
class Node
{
protected:
    Node * left;///<Left daughter node
    Node * right;///<Right daughter node
    Node * father;///<Father node
    double age;///<Age of the node
    double tempdist;///<distance to parent; only used for initial construction
    double weight,weight2;///< Certainty of a node (2 representations)
    bool intree;///< Whether in the tree, or within a population
    int id;///<Index of the node
    double lefting;///<Distance from the far left on a (0,1) scale
    double totheleft;///<Number of nodes to the left of this node
    double numchildren;///<Number of children of this node
    int popnodeid;///<Node is a population node
    Node * popnode;///<Population node of this node
//construction functions    
    int getStartPt(string * newick);///gets the correct colon
    int getEndPt(string * newick,char value=':', unsigned int i0=0);///gets the correct bracket
    double getAge(string * newick,int found,int found2);/// gets the age from the remaining string
public:
    Node();
    Node(Data *d,std::string * newick,Node * father,std::vector<Node*> * nodes);
    ~Node();
    
    std::string newick(Data *d, int p=6) const;///<Convert the node to a Newick string (precision options)
    std::vector<int> tipsUnder();///< A list of the ids under this node
    std::vector<int> internalsUnder();///< A list of the internal ids under this node
    void tellPopNodes(Node * newpopnode);///< Updates whether each node is a popululation node, or within a population
    
    inline bool isPopNode(std::vector<int> indsinpop){
      std::vector<int> tips=tipsUnder();
      if(tips.size()!=indsinpop.size())return(false);
      for(unsigned int c1=0;c1<tips.size();c1++)if(!isinlist(tips[c1],&indsinpop))return(false);
      for(unsigned int c1=0;c1<indsinpop.size();c1++)if(!isinlist(indsinpop[c1],&tips))return(false);
      return(true);
    };///< returns true if all and only those ids listed are under this node
    inline bool isinlist(int i,std::vector<int> *list){
      for(unsigned int c1=0;c1<list->size();c1++){if(i==list->at(c1))return(true);
      }
      return(false);
    };
    
    double getDist(std::vector<std::vector<double>* > m,std::vector<int> tips);///<Gets the distance weighted by m of the tips to the diagonal
    std::vector<std::vector<double>* > swapOrder(std::vector<std::vector<double>* > m,std::vector<int> lu,std::vector<int> lr);///< Swaps the left and right of a node based on m
    double diagonalise(std::vector<std::vector<double>* > m);///< Diagonalises the tree according to the distance from the diagonal weighted by m
    void setcertainty(std::vector<std::vector<double> > *mat);///< Assigns the certainty of a split under the mcmc sample
    inline void assignCertainty(double c,bool isfirst=true){
	if(isfirst) weight=c;
	else weight2=c;
    }///<Sets the certainty externally
    inline int getId() const
    {
        return id;
    }///<Returns the index of the node
    inline void setId(int i)
    {
        id=i;
    }///<Sets the index of the node
    inline double getDist() const
    {
        if (father==NULL)
            return 0;
        return father->age-age;
    }///<Returns the distance of the node to its father
    inline double getTempDist() const
    {
	return(tempdist);
    }///<Returns the distance of the node to its father
    inline double getAge() const
    {
        return age;
    }///<Returns the age of the node
    inline void setAge(double a)
    {
        age=a;
    }///<Sets the age of the node
    inline Node * getLeft() const
    {
        return left;
    }///<Returns the left daughter node
    inline Node * getRight() const
    {
        return right;
    }///<Returns the right daughter node
    inline Node * getFather() const
    {
        return father;
    }///<Returns the father node
    inline void setFather(Node * f)
    {
        father=f;
    }///<Sets the father of the Node
    inline void setLeft(Node * l)
    {
        left=l;
    }///<Sets the left child of the Node
    inline void setRight(Node * r)
    {
        right=r;
    }///<Sets the right child of the Node
    inline void changeAge(double a)
    {
        age=age+a;
    }///<Changes the age of the node by adding an amount a
    inline double getWeight(bool first=true){
      if(first) return(weight);
      return(weight2);
    }
    inline void setInTree(){
	intree=true;
    }///<Sets that this node is a real tree node
    inline bool inTree(){
	return(intree);
    }///<Is this node a real tree node?
    double calcChildrenBelow(double popscale=1.0);
    inline double getChildrenBelow(){return(numchildren);};
    double calcToTheLeft(double totheleft=0,double popscale=1.0);
    inline double getToTheLeft(){return(totheleft);};
    void calcLefting(double N,bool popnodesonly);
    inline void setLefting(double l){
      lefting=l;
    }
    inline double getLefting(){return(lefting);};
    inline int getPopNodeId(){
      return(popnodeid);
    }
    inline bool isPopNode(){
      if(popnodeid==id) return(true);
      return(false);
    };///< Is this node the root for a population?
    inline bool isUnderPopNode(){
      if(popnodeid>0) return(true);
      return(false);
    };///<Is this node within a population?
};

} // end namespace fines
#endif
