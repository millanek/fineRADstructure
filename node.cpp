#include "node.h"
#include <sstream>
#include <math.h>

using namespace std;
namespace fines
{

Node::Node()
{
    left=NULL;
    right=NULL;
    father=NULL;
    tempdist=0.0;
    age=0.0;
    weight=-1.0;
    weight2=-1.0;
    numchildren=1;
    totheleft=0;
    lefting=-1;
    intree=false;
    popnodeid=-1;
    popnode=NULL;
    id=0;
}

Node::Node(Data *d,string * newick,Node * father,vector<Node*> * nodes)
{
    this->father=father;
    int found=getStartPt(newick);
    int found2=getEndPt(newick);    
    int foundbracket=found2;
    for(int i=found2;i>found;i--) {
      if(newick->at(i)==')') {foundbracket=i; break;}  // if not a "):" sequence, we have an edge length
    }
    tempdist=0.0;
    age=0.0;
    weight=-1.0;
    weight2=-1.0;
    numchildren=-1;
    totheleft=0;
    lefting=-1;
    popnodeid=-1;
    intree=true;
    popnode=NULL;
    // set age
    if(foundbracket<found2-1) {
      weight=1.0-getAge(newick,foundbracket,found2); // getage is actually just a string->char conversion!
      tempdist=getAge(newick,found2,newick->length());
    }else{
      tempdist=getAge(newick,found2,newick->length());
    }
    if (found2==0)found2=newick->length(); // now read till the end of the string
    if (found==0)
    {
        //Leaf
	age=0.0;
        left=NULL;
        right=NULL;
        istringstream input(*newick);
	string name =newick->substr(found,found2);
//        input >> name;
	id=d->getIndex(name);//nodes->size();
        nodes->at(id)=this;
    }
    else
    {
        //Internal node
        id=-1;
        string leftStr =newick->substr(1,found-1);
        string rightStr=newick->substr(found+1,found2-2-found);
        left=new Node(d,&leftStr,this,nodes);
        right=new Node(d,&rightStr,this,nodes);
	age = left->getAge() + left->getTempDist();
        nodes->push_back(this);
    }
}

int Node::getStartPt(string * newick){
   //Find central comma if it exists, or the start of the string if not
    int depth=0;
    int found=0;
    for (unsigned int i=0;i<newick->length();i++)
    {
        if (newick->at(i)=='(')
        {
            depth++;
            continue;
        }
        if (newick->at(i)==')')
        {
            depth--;
            continue;
        }
        if (newick->at(i)==',' && depth==1)
        {
            found=i;
            break;
        }
    }
    return(found);
}

int Node::getEndPt(string * newick,char value, unsigned int i0){
    //Find last value = ':' by default, returns 0 if not found
    int found2=0;
    int depth=0;
    for (unsigned int i=i0;i<newick->length();i++)
    {
        if (newick->at(i)=='(')
            depth++;
        if (newick->at(i)==value)
            found2=i;
        if (newick->at(i)==')')
        {
            depth--;
            if (depth==0 && value!=')')
                found2=i;
        }
    }
    
    return(found2);
}

double Node::getAge(string * newick,int found,int found2)
{
   double rval=0.0;
    if (found2>0)
    {
//	string tstr=newick->substr(found2+1,newick->length()-found2);
	string tstr=newick->substr(found+1,found2-found-1);
        istringstream input(tstr);
	input.precision(10);
        input >> rval;
	cout.precision(10);
    }
    return(rval);
}

string Node::newick(Data * d, int p) const
{
    ostringstream idstream;
    ostringstream diststream;
    ostringstream weightstream;
    diststream.setf(ios::scientific, ios::floatfield);
    diststream.precision(p);
    if(id<d->getDim()) idstream << d->getnames(id);
    diststream<<getDist();
    weightstream.setf(ios::fixed);
    weightstream.precision(3);
    if(intree && weight>=0) {
      // only non-population edges are "in the tree"; don't print edge weights for them if they are not
      weightstream<<1-weight;
      if(weight2>=0) weightstream<<"|"<<1-weight2;
    }
    if((left==NULL && right!=NULL)||(left!=NULL && right==NULL) ) {cerr<<"Unmatched node A"<<endl;throw(std::string("Unmatched node!"));}
    if (left==NULL)
        return idstream.str()+":"+diststream.str();
    else{
		string leftstring=left->newick(d,p);
		string rightstring=right->newick(d,p);
		string ret="("+leftstring+","+rightstring+")"+weightstream.str()+idstream.str()+":"+diststream.str();
//		cout<<ret<<endl;
        return ret;
	}
}

vector<int> Node::tipsUnder()
{
  if((left==NULL && right!=NULL)||(left!=NULL && right==NULL) ) {cerr<<"Unmatched node B"<<endl;throw(std::string("Unmatched node!"));}
    if (left==NULL) return(vector<int>(1,id));
    vector<int> lu=left->tipsUnder();
    vector<int> ru=right->tipsUnder();
    for(unsigned int c1=0;c1<ru.size();c1++) lu.push_back(ru[c1]);

    return(lu);
}

vector<int> Node::internalsUnder()
{
  if((left==NULL && right!=NULL)||(left!=NULL && right==NULL) ) {cerr<<"Unmatched node C"<<endl;throw(std::string("Unmatched node!"));}
    if (left==NULL) return(vector<int>(0,0));
    vector<int> ret;
    ret.push_back(id);
    vector<int> lu=left->internalsUnder();
    vector<int> ru=right->internalsUnder();
    for(unsigned int c1=0;c1<lu.size();c1++) ret.push_back(lu[c1]);
    for(unsigned int c1=0;c1<ru.size();c1++) ret.push_back(ru[c1]);
    return(ret);
}


std::vector<std::vector<double>* > Node::swapOrder(std::vector<std::vector<double>* > m,std::vector<int> lu,std::vector<int> lr)
{
	Node * tnode;
	std::vector< std::vector<double>* > tmp;
	for(unsigned int c1=0;c1<lu.size();c1++) tmp.push_back(m[lu[c1]]);
	for(unsigned int c1=0;c1<lr.size();c1++) tmp.push_back(m[lr[c1]]);
	for(unsigned int c1=0;c1<lr.size();c1++) m[lr[c1]]=tmp[c1+lu.size()];
	for(unsigned int c1=0;c1<lu.size();c1++) m[lu[c1]]=tmp[c1];
	tnode=left;
	left=right;
	right=tnode;
	return(m);
}

double Node::getDist(std::vector<std::vector<double>* > m,vector<int> tips)
{
	double d=0.0;
	for(unsigned int c1=0;c1<tips.size();c1++){
		for(unsigned int c2=c1+1;c2<tips.size();c2++){
			d+=m[tips[c1]]->at(tips[c2]) * 4.0*(tips[c1]-tips[c2])*(tips[c1]-tips[c2]);
		}
	}
	return(d);
}	

double Node::diagonalise(std::vector<std::vector<double>* > m){
   if((left==NULL && right!=NULL)||(left!=NULL && right==NULL) ) {cerr<<"Unmatched node D"<<endl;throw(std::string("Unmatched node!"));}
	//cout<<"Node "<<id<<":"<<endl;
    if (left==NULL) // tip node
	return(0.0);
    else{
//        double dl=left->diagonalise(m);
//	double dr=right->diagonalise(m);
	vector<int> lu=left->tipsUnder();
    	vector<int> ru=right->tipsUnder();
	vector<int> tips=lu;
	for(unsigned int c1=0;c1<ru.size();c1++) tips.push_back(ru[c1]);
	double startd=getDist(m,tips);
//cout<<"dl="<<dl<<" and dr="<<dr<<" and dt="<<dl+dr<<" and dist="<<startd<<endl;
	m=swapOrder(m,lu,ru);
	double otherd=getDist(m,tips);
	if(startd<otherd) {
		m=swapOrder(m,ru,lu);
		return(startd);
	}else return(otherd);
    }
}

void Node::setcertainty(std::vector<std::vector<double> > *mat)
{
	vector<int> tips=tipsUnder();
//	cout<<"Id="<<id<<" tipsize="<<tips.size()<<endl;
	double count=0;
	double denom=0;
	for(unsigned int c1=0;c1<tips.size();c1++) {
	  for(unsigned int c2=0;c2<tips.size();c2++) if(c1!=c2) count +=mat->at(tips[c1])[tips[c2]];
	  for(unsigned int c2=0;c2<mat->at(tips[c1]).size();c2++) if(tips[c1]!=(int)c2)denom+=mat->at(tips[c1])[c2];
	}
	//if(tips.size()>1) weight=count / (tips.size() * (tips.size()-1.0)/2.0);
	if(tips.size()>1 && denom>0) weight=count/denom;
	if(left!=NULL) left->setcertainty(mat);
	if(right!=NULL) right->setcertainty(mat);
}

double Node::calcChildrenBelow(double popscale)
{
  
  double numchildrenret=0;
  if(left!=NULL) numchildrenret+=left->calcChildrenBelow(popscale);
  if(right!=NULL) numchildrenret+=right->calcChildrenBelow(popscale);
  if(left==NULL &&right==NULL) numchildrenret=1.0;
  
  if(popscale==1.0) numchildren=numchildrenret;
  else numchildren=pow(numchildrenret,popscale);
  return(numchildrenret);
}

/*
double Node::calcToTheLeft(double parenttotheleft,double popscale)
{  
  totheleft=parenttotheleft;
  if(left==NULL &&right==NULL) totheleft+=1;
  else if(left!=NULL) totheleft=left->calcToTheLeft(totheleft);
  double ret=totheleft;
  if(right!=NULL) ret=right->calcToTheLeft(totheleft);
  return(ret);
}*/


double Node::calcToTheLeft(double parenttotheleft,double popscale)
{  
  double underleft=0,underright=0;
  if(left==NULL &&right==NULL) {
    underleft=1;
  }else if(left!=NULL) {
    underleft=left->calcToTheLeft(parenttotheleft,popscale)-parenttotheleft;
  }
  //double ret=parenttotheleft + underleft;
  if(right!=NULL) underright=right->calcToTheLeft(parenttotheleft + underleft,popscale)-(parenttotheleft + underleft);
  double ret=parenttotheleft + underleft + underright;
  if(isPopNode()){
//     totheleft=parenttotheleft + (underleft+underright)/2.0;
    totheleft=parenttotheleft + numchildren/2.0;
    ret = parenttotheleft + numchildren;
//    cout<<"POPNODE "<<id<<" and I'm at "<<totheleft<<" and have "<<underleft<<" left and "<<numchildren<<" beneath"<<" and parenttotheleft="<<parenttotheleft<<endl;    //if(popscale!=1.0) totheleft=parenttotheleft + pow(underleft+underright,popscale)/2.0;
    //else totheleft=parenttotheleft + (underleft+underright)/2.0;
  }else{
    totheleft=parenttotheleft + underleft;
//    if(isUnderPopNode()) cout<<"INTERNAL Node "<<id<<" has totheleft="<<totheleft<<endl;
//    else cout<<"EXTERNAL Node "<<id<<" has totheleft="<<totheleft<<endl;
  }
  return(ret);
//  if(right!=NULL) underright=right->calcToTheLeft(totheleft,popscale) - totheleft;
//  return(totheleft + underright);
}

void  Node::calcLefting(double N,bool popnodesonly)
{
  double offset=0.0;
  lefting=(totheleft+offset)/N;
  if(isPopNode()){
    lefting=(totheleft+offset)/N;
//    cout<<"POPNODE... "<<id<<" has lefting="<<lefting<<" from N="<<N<<endl;
  }else if(isUnderPopNode()){
    if(popnode==NULL) throw(string("Node error: inconsistency with the popnodes?"));
    if(popnodesonly) lefting=popnode->getLefting();
    else {
      lefting=(totheleft-0.5)/N;
    }
  }else if(left!=NULL && right!=NULL) {
   // lefting=(totheleft  - left->getChildrenBelow() +getChildrenBelow()/2.0)/N;
    lefting=(totheleft+offset)/N;
 //   cout<<"EXTERNAL Node "<<id<<" has lefting="<<lefting<<" from N="<<N<<endl;
    
  }
  //  else if(left==NULL  || right==NULL) throw(string("Illegal nodeconfiguration?!"));
//  lefting=0.5*(left->getToTheLeft()+right->getToTheLeft())/N;*/
  if(left!=NULL)left->calcLefting(N,popnodesonly);
  if(right!=NULL)right->calcLefting(N,popnodesonly);
}

Node::~Node()
{
}

void Node::tellPopNodes(Node * newpopnode)
{
  popnode=newpopnode;
  popnodeid=newpopnode->getId();
  if(left!=NULL)left->tellPopNodes(newpopnode);
  if(right!=NULL)right->tellPopNodes(newpopnode);
  
}

} // end namespace fines

