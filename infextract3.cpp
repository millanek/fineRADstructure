#include "infextract3.h"
#include <iomanip>
#include <iostream>

#define MIN_LINK 0.000001
using namespace std;
namespace fines
{

InfExtract3::InfExtract3(Data *d,FsXml *fs,std::vector<Node*> nodes,bool v)
{
	weights=vector<double>(nodes.size(),0.0);
	weights2=vector<double>(nodes.size(),0.0);
	counts=0;

	data=d;
	verbose=v;
	streampos fpos;
	if(verbose) cout<<"Creating average coincidence..."<<endl;
	fpos=fs->gotoLineContaining("<Iteration>");
	int counts=0;
	while(!fs->eof() && fpos>=0) {
		if(fpos>0) state=new State(data,fs,vector<double>(1.0),1.0,BETAMOD_CONST,2.0);
		vector<double> tmp=getSplits(state,nodes,true);//****
		vector<double> tmp2=getSplits(state,nodes,false);//****
		for(unsigned int c1=0;c1<tmp.size();c1++) { weights[c1]+=tmp[c1]; weights2[c1]+=tmp2[c1];}
		counts++;
		delete(state);
		fpos=fs->gotoNextLineContaining("<Iteration>");
	}
	for(unsigned int c1=0;c1<nodes.size();c1++) if(nodes[c1]->getFather()!=NULL){
	  if(nodes[c1]->inTree() || nodes[c1]->getFather()->inTree()){
		nodes[c1]->assignCertainty(weights[c1]/(double)counts,true);
	//	nodes[c1]->assignCertainty(weights2[c1]/(double)counts,false);
	// NOTE THe weights2 is not functioning correctly.
	  }
	}
	if(verbose) cout<<"Finished creating mean distance matrix."<<endl;
}


double InfExtract3::getSingleSplit(State *state, vector<int> tips)
{
	unsigned int tipsfound=0;
	vector<int> numfound(state->getP(),0);
	for(int c1=0;c1<state->getP();c1++) {
		vector<int> v1=state->getIndInPop(c1);
		bool stateisin=false;
		for(unsigned int c2=0;c2<v1.size();c2++) {
		  bool thisfound=false;
		  for(unsigned int c3=0;c3<tips.size();c3++) {
			if(v1[c2]==tips[c3]) {
			  if(c2==0) stateisin=true;
			  else if(!stateisin) numfound[c1]++;//return(0);// this tip is in but the state is not
			  tipsfound++;
			  c3=tips.size();
			  thisfound=true;
			}
		  }
/*		  if(!thisfound){// we looked through all the tips
			if(stateisin && strict)return(0);// this state is in but this tip is not
		  }*/
		}
	}
	double maxagree=0;
	for(unsigned int c1=0;c1<numfound.size();c1++) {
	  double agreelevel=2.0*numfound[c1]/(tips.size()+state->getPsize(c1));
	  if(maxagree<agreelevel) maxagree=agreelevel;
	}
	return(maxagree);
}

double InfExtract3::getSingleSplitStrict(State *state, vector<int> tips)
{
	unsigned int tipsfound=0;
	for(int c1=0;c1<state->getP();c1++) {
		vector<int> v1=state->getIndInPop(c1);
		bool stateisin=false;
		for(unsigned int c2=0;c2<v1.size();c2++) {
		  bool thisfound=false;
		  for(unsigned int c3=0;c3<tips.size();c3++) {
			if(v1[c2]==tips[c3]) {
			  if(c2==0) stateisin=true;
			  else if(!stateisin) return(0);// this tip is in but the state is not
			  tipsfound++;
			  c3=tips.size();
			  thisfound=true;
			}
		  }
		  if(!thisfound){// we looked through all the tips
			if(stateisin)return(0);// this state is in but this tip is not
		  }
		}
		if(tipsfound==tips.size()) return(1);
	}
	return(1);
}

vector<double> InfExtract3::getSplits(State *state, std::vector<Node*> nodes,bool strict)
{
	vector<double> ret=vector<double>(nodes.size(),0.0);
	for(int c1=0;c1<(int)nodes.size();c1++) {
		if(nodes[c1]->inTree()&& strict) ret[c1]+=getSingleSplitStrict(state,nodes[c1]->tipsUnder());
		if(nodes[c1]->inTree()&& !strict) ret[c1]+=getSingleSplit(state,nodes[c1]->tipsUnder());
	}
	return(ret);
}

InfExtract3::~InfExtract3()
{
}



} // end namespace fines
