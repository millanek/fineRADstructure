#include "inf1.h"
#include "infmcmc.h"
#include <iomanip>
#include <iostream>

#define MIN_LINK 0.00001
#define BASE_LINK 1
using namespace std;
namespace fines
{

Inf1::Inf1(Data *d,State *s,Data *d2,int datainference,bool v,int tmax,int treescale)
{
	warned=false;
	test_max=tmax;
	data=d;
	state=s;
	verbose=v;
//	reorderSuper(state);
	this->treescale=treescale;
	currentheight=0;
	nodes=vector<Node*>(d->getDim()  + d->getDim() - d->numIgnore()-1);
  for (unsigned int i=0;i<nodes.size();i++) {
    nodes[i]=new Node();
    nodes[i]->setId(i);
  }
  root=nodes[d->getDim()  + d->getDim() - d->numIgnore()-2];
  createSkeletonTree();
    if(verbose) state->setprint(&cout);
	//if(verbose) {cout<<"Printing Data"<<endl;data->print(&cout);}
}

Inf1::Inf1(Data *d,State *s,string newick,bool v)
{
  	warned=false;
	test_max=0;
	data=d;
	state=s;
	verbose=v;
    int n=d->getDim();

    nodes=vector<Node*>(n);
    for(unsigned int c1=0;c1<nodes.size();c1++){nodes[c1]=NULL;}
    root = new Node(d,&newick,NULL,&nodes);
    //cout<<"root="<<n+n-2<<endl;
   // cout<<"but root="<<nodes.size()<<endl;
    for(unsigned int c1=0;c1<nodes.size();c1++){
      if(nodes[c1]==NULL){
	  nodes[c1]=new Node();
      }
      nodes[c1]->setId(c1);
    }
    // Nodes n to 2n-1 are created from the root
//    for(long c1=n;c1<n;c1++){
  //  }
    ///NOTE: The 0..n-1 nodes are NOT INDEXED BY THEIR ID. Sort HERE if this is needed in the future.

    if(verbose) state->setprint(&cout);


}

Inf1::Inf1(Data *d,int initpop,Data *d2,int datainference,int modeltype,bool v,int tmax)
{
	warned=false;
	test_max=tmax;
	verbose=v;
	data=d;
	if(initpop<0 || initpop>data->getDim()) initpop=data->getDim();
	state=new State(data,initpop,vector<double>(1.0),1.0,BETAMOD_CONST,2.0,NULL,-1,modeltype);
	state->setVerbose(verbose);
  int n=d->getDim() - d->numIgnore();
  nodes=vector<Node*>(n+n-1);
  for (int i=0;i<n+n-1;i++) {
    nodes[i]=new Node();
    nodes[i]->setId(i);
  }
  root=nodes[n+n-2];
  createSkeletonTree();

    if(verbose) state->setprint(&cout);
	//if(verbose) {cout<<"Printing Data"<<endl;data->print(&cout);}
}

void Inf1::createSkeletonTree()
{
  popnodes=vector<int>(state->getP());
  assignednodes=data->getDim();

  for(int i=0;i<state->getP();i++) {
	vector<int> plist=ignoreIgnorables(state->getIndInPop(i));
	for(unsigned int j=1;j<plist.size();j++){
	  linkNodes(plist[0],plist[j],MIN_LINK,false,false);
	}
	if(plist.size()>0){
	  Node *fi=nodes[plist[0]];
	  while(fi->getFather()!=NULL) {fi=fi->getFather();}//fi=nodes[fi->getId()]->getFather(); }
	  popnodes[i]=fi->getId();
	}else popnodes[i]=-1;
  }
}

void Inf1::statusOut(long c1, long numints)
{
	long maxdisp=myMin(numints,50);
	int gap=((numints)/maxdisp);
	int numhash=1;
	if ((c1)%gap==0){
	  //if(numints<maxdisp) while(numhash<c1*(1.0-1.0/gap)) numints++;// this should account for doing less than 50 total prints.  But it doesn't!
		ios::fmtflags f = cout.flags();
//if (maxdisp*c1/(numints)>0) 
		cout<<"\b\b\b\b\b";
		for(int i=0;i<numhash;i++) cout<<"#";
		cout.width(4);
		cout.precision(0);
		cout.fill(' ');
		cout<<right<<fixed<<100l*(double)c1/(numints)<<"%"<<flush;
		cout.flags(f); // reset flags
	}
	if (c1+1==numints) cout<<"\b\b\b\b\b# 100%"<<endl<<flush;

}

void Inf1::printTree(ostream * fileout,bool xml)
{
	if(xml) *fileout<<"<Tree>";
	*fileout<<root->newick(data)<<";";
	if(xml) *fileout<<"</Tree>";
	*fileout<<endl;
}

void Inf1::setcertainty(std::vector<std::vector<double> > *mat)
{
	if((int)mat->size()!=data->getDim()) {cerr<<"setcertainty: wrong matrix size"<<endl;;throw(string("setcertainty: wrong matrix size!"));}
	root->setcertainty(mat);
}

/*
// NOT WORKING
void Inf1::diagonaliseStatesOnM(State *in, std::vector<std::vector<double> > *mat,bool symm)
{
	std::vector<std::vector<double>* > mptr;
	for(unsigned int c1=0;c1<mat->size();c1++) {
		mptr.push_back(& mat->at(c1));
	}

	vector<int> allvec=in->allIndInOrder();
	for(int a=0;a<in->getP();a++) {
	  vector<int> pvec=in->getIndInPop(a);
	  for(int c1=1;c1<pvec.size();c1++) {
	    for(int c2=c1-1;c2>=0;c2--) {
		double d1=root->getDist(mptr,allvec,true);
		mptr=swapM(mptr,pvec[c1],pvec[c2]);
		double d2=root->getDist(mptr,allvec,true);
		cout<<"A d1-s2="<<d1/d2<<endl;
		if(d2+0.0001<d1) in->reorderIndiv(a,c1,c2);
		else mptr=swapM(mptr,pvec[c1],pvec[c2]);
	    }
	  }
	}
}*/

void Inf1::diagonaliseOnM(std::vector<std::vector<double> > *mat,bool symm)
{
	std::vector<std::vector<double>* > mptr;
	for(unsigned int c1=0;c1<mat->size();c1++) {
		mptr.push_back(& mat->at(c1));
	}
	if(!symm){
		for(unsigned int c1=0;c1<mat->size();c1++) {
		for(unsigned int c2=0;c2<mat->size();c2++) {
			mat->at(c1)[c2]=(mat->at(c1)[c2] + mat->at(c2)[c1])/2.0;
		}
		}
	}
	root->diagonalise(mptr);
	//std::vector<int> tu=root->tipsUnder();
	//state->setIndInPop(root->tipsUnder(),0);
//	reorderState();
}

bool Inf1::inOrder (int i,int j) 
{ // i and j are individuals
	int pop=state->getPop(i);
	if(verbose) cout<<"Comparing "<<i<<" and "<<j<<" in pop "<<pop<<endl;
	if(state->getPop(j)!=pop) {cerr<<"inOrder:Comparing individuals in different populations!"<<endl;throw(std::string("inOrder:Comparing individuals in different populations!"));}
	vector<int> plist=state->getIndInPop(pop);
	for(unsigned int c1=0;c1<plist.size();c1++) {
		if(plist[c1]==i) return(true);
		else if (plist[c1]==j) return(false);
	}
	cerr<<"inOrder:states not found!"<<endl;
	throw(std::string("inOrder: states not found!"));
/*	int i1=-1,i2=-1,itest=0;
	while(i1<0 || i2<0) {
		if(plist[itest]==i)i1=itest;
		if(plist[itest]==j)i2=itest;
		itest++;
		if(itest>=(int)plist.size()) break;
	}
	if(i1<0 || i2<0) throw(std::string("inOrder:compared individuals not found!"));
	return (i1<i2);*/
}

bool Inf1::inFullOrder (int i,int j) 
{ // i and j are indivs
	// find the state that occurs first
	vector<int> allindiv=state->allIndInOrder();
	for(unsigned int c1=0;c1<allindiv.size();c1++) {
		if(allindiv[c1]==i) return(true);
		else if (allindiv[c1]==j) return(false);
	}
	cerr<<"statesInOrder: states not found!"<<endl;
	throw(std::string("statesInOrder: states not found!"));
}

void Inf1::reorderSuper(State *in)
{// puts all the super individuals at the start
  int nfound=0;
  for(int c1=0;c1<in->getP();c1++) {
	  if(testIgnorePop(c1)) {
	    if(nfound<c1) {
		  if(verbose) cout<<"Reordering "<<c1<<" ("<<data->getnames(in->getIndInPop(c1)[0])<<") and "<<nfound<<" ("<<data->getnames(in->getIndInPop(nfound)[0])<<")"<<endl;
		  in->reorderPop(c1,nfound);
		  nfound++;
	    }
	  }else if(verbose)cout<<"Number"<<c1<<" ("<<data->getnames(in->getIndInPop(c1)[0])<<") ok!"<<endl;
	  if(verbose) {
		for(int i=0;i<=c1;i++) cout<<data->getnames(in->getIndInPop(i)[0])<<",";
	  	cout<<endl;
	  }
  }
}

void Inf1::reorderState(State *in)
{
	if(in==NULL) in=state;
	// reorder the individuals within a state
/*	for(int a=0;a<in->getP();a++){
		if(in->getIndInPop(a).size()!=in->getPsize(a)) throw("Error: indinp disagrees wih psize!");
//cout<<"a="<<a<<" of "<<in->getP()<<endl;
	  for(int c1=0;c1<in->getPsize(a)-1;c1++) {
//cout<<"c1="<<c1<<" of "<<in->getPsize(a)-1<<endl;
	    int npos=c1-1;
	    while(npos>=0) {
//		cout<<"npos="<<npos<<" of "<<flush;
		if(!inOrder(in->getIndInPop(a)[npos-1],in->getIndInPop(a)[c1])) {
			npos++;break;
		}else npos--;
	    }
cout<<"Stopped with npos="<<npos<<endl;
	    if(npos<c1 && npos>=0) {
		in->reorderIndiv(a,c1,npos);
		cout<<"Reordering from "<<c1<<" to "<<npos<<endl;
		}
	  }
	}
cout<<"Done with indivs"<<endl;*/
	// reorder the states using insertion sort (not fast, but doesn't matter)
	for(int c1=1;c1<in->getP();c1++) {
	  int npos=c1-1;
	  if(!inFullOrder(in->getIndInPop(npos)[0],in->getIndInPop(c1)[0])) {
	    while(npos>=0) {
		  if(verbose) cout<<"Comparing "<<c1<<" ("<<data->getnames(in->getIndInPop(c1)[0])<<") and "<<npos<<" ("<<data->getnames(in->getIndInPop(npos)[0])<<")"<<endl;
		  if(!inFullOrder(in->getIndInPop(c1)[0],in->getIndInPop(npos)[0])) {
			  npos++;break;
		  }else npos--;
	    }
	    if(npos<c1 && npos>=0) {
		  if(verbose) cout<<"Reordering "<<c1<<" ("<<data->getnames(in->getIndInPop(c1)[0])<<") and "<<npos<<" ("<<data->getnames(in->getIndInPop(npos)[0])<<")"<<endl;
		  in->reorderPop(c1,npos);
	    }
	  }else if(verbose)cout<<"Number"<<c1<<" ("<<data->getnames(in->getIndInPop(c1)[0])<<") ok!"<<endl;
	  if(verbose) {
		for(int i=0;i<=c1;i++) cout<<data->getnames(in->getIndInPop(i)[0])<<",";
	  	cout<<endl;
	  }
	}
}

void Inf1::mergeHillClimb(ostream * fileout,  bool stopattop,int treemodification)
{
//  cout<<"MERGEHILLCLIMB MOD="<<treemodification<<" stopattop="<<stopattop<<endl;
	double pprob=state->posteriorProb();
	if(fileout!=NULL) *fileout<<"<Posterior>"<<pprob<<"</Posterior>"<<endl;
	bool stillgoing=true;
	if(state->getP()==1) stillgoing=false;
	long c1=0;
	long numints=0;
	if(!stopattop) numints=state->getP()-1-data->numIgnore();
	initheight=-INFINITY;
	if(treemodification==2) doHyperParameterMaximisation(0);
	while(stillgoing) {
		if(!stopattop) statusOut(c1,numints);
		c1++;
 		stillgoing=doBestMerge(test_max,stopattop,treemodification);
		if(verbose && stillgoing)cout<<"Merge attempt successful."<<endl;
		else if(verbose && !stillgoing)cout<<"Merging finished."<<endl;
		if(state->getP()==data->numIgnore()+1) stillgoing=false;
		if(verbose) {
		cout<<"Accepted"<<endl;
		state->setprint(&cout);
		cout<<"Posterior 1 ="<<state->posteriorProb()<<endl;
		}
		if(treemodification>0 && stillgoing) {state->setSuperIndivRule(true);}
		if(treemodification==2) doHyperParameterMaximisation(100);
	}
	if(fileout!=NULL) printTree(fileout);
	if(treemodification>0) state->setSuperIndivRule(false);// reset things
}

void Inf1::splitHillClimb(bool stopattop)
{
	bool stillgoing=true;
	if(state->getP()==data->getDim()-data->numIgnore()) stillgoing=false;
	while(stillgoing) {
		stillgoing=doBestSplit(stopattop);
		if(state->getP()==data->getDim()) stillgoing=false;
		if(verbose && stillgoing)cout<<"Split attempt successful."<<endl;
		else if(verbose && !stillgoing)cout<<"Splitting finished."<<endl;
		if(verbose) {
		  cout<<"Accepted"<<endl;
		  state->setprint(&cout);
		  cout<<"Posterior="<<state->posteriorProb()<<endl;
		}
	}
}

void Inf1::doHyperParameterMaximisation(int numsteps){
    InfMCMC tmcmc(data,state,NULL,INFDATA_COUNTS,0,false);
    for(int c1=0;c1<numsteps;c1++) {
      tmcmc.moveHyper(true);
    }
	tmcmc.exportXmlIter(&cout,0); // we do not have a natural height available
    state=new State(tmcmc.getState());
}

bool Inf1::doBestMerge(int mmax,bool stopattop,int treemodification)
{
	if(treemodification>0)state->setSuperIndivRule(true);
	State * olds=new State(state);
	State * bests=new State(state);
	State * origstate=new State(state);
	vector <int> li,lj;
	double prevpost=state->posteriorProb();

// construct a list of merges to consider
	if(mmax<0 || (state->getP()-data->numIgnore())*(state->getP()-data->numIgnore()-1)/2 < mmax) {
	  if(verbose) cout<<"Comparing all merges"<<endl;
	  //for(int i=state->getP()-1;i>=data->numIgnore();i--) for(int j=i-1;j>=data->numIgnore();j--) {li.push_back(i);lj.push_back(j);}
	  for(int i=state->getP()-1;i>=0;i--) for(int j=i-1;j>=0;j--) {li.push_back(i);lj.push_back(j);}
	}else {
	  if(!warned) {warned=true;
		cout<<"WARNING!  NOT TESTING ALL "<<state->getP()*(state->getP()-1)/2<<" COMBINATIONS! (max "<<mmax<<")"<<endl;
	  }
	  for(int i=0;i<mmax;i++) {
		int a=0,b=0;
		while(a==b) {
			a=RandomInteger(data->numIgnore(),state->getP()-1);
			b=RandomInteger(data->numIgnore(),state->getP()-1);
			if(a>b) {int c=a;a=b;b=c;}
		}
		li.push_back(a);
		lj.push_back(b);
	  }
	}
// test the merges to find the best
	int curbest=-1;
	double bestpost;
	if(stopattop) bestpost=prevpost;
	else bestpost=-10e30;
	if(!(initheight>-INFINITY)){
	    initheight=state->posteriorProb();
	}
	for(unsigned int i=0;i<li.size();i++) {
		if(testIgnorePop(li[i])||testIgnorePop(lj[i])) {
		  continue;
		  //cerr<<"Error:doBestMerge Trying to merge a superpop!"<<endl;
		  //throw(string("Error in doBestMerge: Merging a superpopulation!"));
		}
		string s1=data->getnames(state->getIndInPop(li[i])[0]);
		string s2=data->getnames(state->getIndInPop(lj[i])[0]);
		state->merge(li[i],lj[i]);
		double curpost=state->posteriorProb();

		if(curpost > bestpost) 
		{// this is our current best
			if(verbose) {cout.precision(9);cout<<"Found "<<li[i]<<"("<<s1<<")"<<","<<lj[i]<<"("<<s2<<")"<<": " <<curpost<<" vs "<<bestpost<<endl;}
			curbest=i;
			bestpost=curpost;
			delete(bests);
			bests=new State(state);
		}//else if(verbose) {cout.precision(9);cout<<"REJECT "<<li[i]<<"("<<s1<<")"<<","<<lj[i]<<"("<<s2<<")"<<": " <<curpost<<" vs "<<bestpost<<endl;}


		if(state!=NULL) delete(state);
		state=new State(olds);
	}
	if(treemodification>0)state->setSuperIndivRule(false);
	if(curbest>=0) {
		delete(state);
		state=bests;
/*		olds->setSuperIndivRule(false);
		state->setSuperIndivRule(false);
		prevpost=olds->posteriorProb();
		bestpost=state->posteriorProb();*/
		delete(olds);
		if(verbose) cout<<"Merging "<<li[curbest]<<" and "<<lj[curbest]<<" with posterior "<<bestpost<<"-"<<prevpost<<"="<<bestpost-prevpost<<endl;
		if(stopattop) linkNodes(li[curbest],lj[curbest],MIN_LINK,true,true);
 		else if(treescale==1) {
		  currentheight=currentheight+BASE_LINK;
		  linkNodes(li[curbest],lj[curbest],currentheight,true,true,true);// force absolute age
		}else if(treescale==2){
		  double relprob=origstate->posteriorProb();// where "0" would be for us now with the flattened data
		  double effectivestep=relprob-bestpost;
		  //if(bestpost-initheight+BASE_LINK>currentheight) currentheight=bestpost-initheight+BASE_LINK;// upside down
		  if(effectivestep>BASE_LINK) currentheight+=effectivestep;
		  //if(initheight-bestpost+BASE_LINK>currentheight) currentheight=-bestpost+initheight+BASE_LINK;
		  else currentheight+=BASE_LINK;	
		 // cout<<setprecision(10)<<"K="<<state->getP()-1<<" with posterior "<<bestpost<<" currentheight "<<bestpost-initheight<<" ="<<currentheight<<" IH:"<<initheight <<endl;
		  linkNodes(li[curbest],lj[curbest],currentheight,true,true,true);
		}else{ //treescale==0 and default:
		  linkNodes(li[curbest],lj[curbest],BASE_LINK,true,true);
		}
		return(true);
	}
	delete(olds);
	return (false);	
}

void Inf1::linkNodes(int i, int j,double d,bool usepopnodes,bool intree,bool absage)
{
  Node *fi ,*fj;
if(i>j){int k=j;j=i;i=k;}
  if(usepopnodes) {
    if(popnodes[i]<0 || popnodes[j]<0){ cerr<<"BAD NODE"<<endl;throw(string("bad node in linknodes"));};
    fi=nodes[popnodes[i]], fj=nodes[popnodes[j]];
    //cout<<"Popi[0]="<<
  }else {fi=nodes[i], fj=nodes[j];}
    //cout<<"BBB fi="<<fi->getId()<<" fj="<<fj->getId()<<" fathernode="<<assignednodes<<endl;
  while(fi->getFather()!=NULL) fi=fi->getFather();
  while(fj->getFather()!=NULL) fj=fj->getFather();
  if(fi==fj){// this is an error at present
    cerr<<"Error: tried to merge a population with itself!"<<endl;
    throw(string("Error: tried to merge a population with itself!"));
  }
  if(assignednodes<=(int)nodes.size()) {
    //cout<<"TTT fi="<<fi->getId()<<" fj="<<fj->getId()<<" fathernode="<<assignednodes<<endl;
    fi->setFather(nodes[assignednodes]);
    fj->setFather(nodes[assignednodes]);
    nodes[assignednodes]->setLeft(fi);
    nodes[assignednodes]->setRight(fj);
    if(!absage){
      if(fi->getAge()>fj->getAge()) d+=fi->getAge();
      else d+=fj->getAge();
    }
    nodes[assignednodes]->setAge(d);
    if(intree) {
	nodes[assignednodes]->setInTree();
	nodes[assignednodes]->getLeft()->setInTree();
	nodes[assignednodes]->getRight()->setInTree();
/*	if(nodes[assignednodes]->getLeft()!=NULL){
		nodes[assignednodes]->getLeft()->getLeft()->setInTree();
		nodes[assignednodes]->getLeft()->getRight()->setInTree();
	}
	if(nodes[assignednodes]->getRight()!=NULL){
		nodes[assignednodes]->getRight()->getLeft()->setInTree();
		nodes[assignednodes]->getRight()->getRight()->setInTree();
	}*/
    }
    if(usepopnodes) {popnodes[i]=assignednodes;
    popnodes.erase(popnodes.begin()+j);}
    assignednodes++;
  }else cerr<<"ERROR! Out of nodes to allocate!"<<endl;
}

bool Inf1::doBestSplit(bool stopattop)
{
	State * olds=new State(state);
	State * bests=new State(state);

// test the splits to find the best
	int curbest=-1;
	double bestpost;
	if(stopattop) bestpost=state->posteriorProb();
	else bestpost=-10e30;
	for(int i=0;i<state->getP();i++) {
	if(verbose) cout<<"Comparing i="<<i<<" of "<<state->getP()<<endl;
	  if(state->getPlength(i)>=2) {
		state->splitSAMSgreedy(i,test_max);
		double curpost=state->posteriorProb();
		if(curpost > bestpost) 
		{// this is our current best
			if(verbose) cout<<"Old "<<curbest<<" with prob "<<bestpost<<". "<<curbest<<" is the new best with posterior="<<bestpost<<endl;
			curbest=i;
			bestpost=curpost;
			delete(bests);
			bests=new State(state);
		}else if(verbose)  cout<<"Old "<<curbest<<" with prob "<<bestpost<<". "<<i<<" is inferior with posterior="<<curpost<<endl;
		if(state!=NULL) delete(state);
		state=new State(olds);
	  }
	}
	delete(olds);
	if(curbest>=0) {
		delete(state);
		state=new State(bests);
		return(true);
	}else return (false);	
}

void Inf1::exportXmlHead(std::ostream * fileout,std::string fs,std::string type,double c,long x)
{
	*fileout<<"<?xml version = '1.0' encoding = 'UTF-8'?>"<<endl<<"<outputFile>"<<endl;
  *fileout<< "<header>"<<endl;
  *fileout<< "<runtype>"<<type<<"</runtype>"<<endl;
  *fileout<<setprecision(12);
  if(c>=0) *fileout<< "<inflation>"<<c<< "</inflation>"<<endl;
  if(x>=0) *fileout<< "<burnin>"<<x<< "</burnin>"<<endl;
  *fileout<<"<datafilename>"<<data->getFileName()<<"</datafilename>"<<endl;
  *fileout<<"<mcmcfilename>"<<fs<<"</mcmcfilename>"<<endl;
  *fileout<< "<copymodel>"<<state->getModelType()<< "</copymodel>"<<endl;    
  *fileout<<"</header>"<<endl;	
}

void Inf1::exportXmlTail(std::ostream * fileout)
{
	*fileout<<"</outputFile>"<<endl;
}


int Inf1::whichPopNode(int i)
{
  return(popnodes[i]);
  
  cout<<"looking for population "<<i<<endl;
  for(unsigned int c1=0;c1<nodes.size();c1++){
    cout<<"Node "<<c1<<endl;
    if(nodes[c1]->getId()==popnodes[i])return(c1);
  }
  cerr<<"WARNING!: Unavailable population node requested!"<<endl;
  return(-1);
}

Inf1::~Inf1()
{
	for(unsigned int c1=0;c1<nodes.size();c1++) if(nodes[c1]!=NULL) delete(nodes[c1]);
	if(state!=NULL) delete(state);
}


} // end namespace fines
