#include "infextract.h"
#include <iomanip>
#include <iostream>

#define MIN_LINK 0.000001
using namespace std;
namespace fines
{

InfExtract::InfExtract(Data *d,FsXml *fs,bool v)
{
	data=d;
	verbose=v;
	streampos fpos;
	int nind=d->getDim();
	for(int c1=0;c1<nind;c1++) {
	  order.push_back(c1);
	  meanX.push_back(vector<double>());
	  for(int c2=0;c2<nind;c2++) {
		meanX[c1].push_back(0.0);
	  }
	}
	if(verbose) cout<<"Creating average coincidence..."<<endl;
	fpos=fs->gotoLineContaining("<Iteration>");
	int counts=0;
	while(!fs->eof() && fpos>=0) {
		counts++;
		if(fpos>0) state=new State(data,fs,vector<double>(1.0),1.0,BETAMOD_CONST,2.0);
		state->addCoincidence(&meanX);
		delete(state);
		fpos=fs->gotoNextLineContaining("<Iteration>");
	}
	for(int c1=0;c1<nind;c1++) {
	  for(int c2=0;c2<nind;c2++) {
		meanX[c1][c2]/=counts;
	  }
	}
	if(counts==0) throw(string("No iterations found in file!"));
	if(verbose) cout<<"Finished creating mean distance matrix."<<endl;
}

void InfExtract::makeMinSquaresState(double pen,State * startstate)
{
	//if(state!=NULL && startstate != state && startstate!=NULL) {cout<<"Deleting State"<<endl;delete(state);}
	//else 
	if(startstate!=NULL) {state=new State(startstate);}
	if(state==NULL) throw("No working state available!");
	double curdistsq=state->getDistanceSq(&meanX),distchange=-1;
	if(verbose) cout<<"Making minimum distance state..."<<endl;
	while(distchange<0) {
		distchange=state->minDistanceMove(&meanX,pen);
		if(distchange<=0) curdistsq=curdistsq+distchange;
	}
	if(verbose) cout<<"Finished creating minimum distance state."<<endl;
}


void InfExtract::printMeanX(std::ostream * out)
{
	*out<<"coincidence,";
	for(unsigned int i=0;i<meanX.size()-1;i++) *out<<data->getnames(order[i])<<", ";
	*out<<data->getnames(order[meanX.size()-1])<<endl;
	for(unsigned int i=0;i<meanX.size();i++) {
	  *out<<data->getnames(order[i])<<", ";
	  for(unsigned int j=0;j<meanX[i].size()-1;j++) *out<<meanX[i][j]<<",";
	  *out<<meanX[i][meanX[i].size()-1]<<endl;
	}

}

void InfExtract::reorder(std::vector<int> allvec)
{
	std::vector<std::vector<double> > oldmeanX=meanX;
	if(allvec.size()< meanX.size()){
	cout<<"allvec.size="<<allvec.size()<<" meanX.size="<<meanX.size()<<endl;
	throw(string("reorder: allvec<meanX!"));
	}
	for(unsigned int c1=0;c1<meanX.size();c1++){
	  if(allvec.size()< meanX[c1].size())throw(string("reorder: allvec<meanX[c1]!"));
	  for(unsigned int c2=0;c2<meanX[0].size();c2++){
		if(allvec[c1]>(int)meanX.size() || allvec[c1]<0 || allvec[c2]>(int)meanX[c1].size() || allvec[c2]<0) throw(string("reorder: Invalid allvec!"));
		meanX[c1][c2]=oldmeanX[allvec[c1]][allvec[c2]];
	  }
	}
	order=allvec;
}

InfExtract::~InfExtract()
{
}


} // end namespace fines
