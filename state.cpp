
#include <algorithm>
#include "state.h"
#include "rng.h"
#include "fsxml.h"
#include <vector>

#include <iostream>
#include <iomanip>

using namespace std;
namespace fines
{
State::State(Data *d,int npop,vector<double> bvec,double a,int betamodel,double corfactor,Data *d2,int datainference,int modeltype)
{
	this->modeltype=modeltype;
	dlength=d2;
	this->datainference=datainference;
	setdataused();
	this->betamodel=betamodel;
	this->corfactor=corfactor;
	hyperprior=bvec;
	data=d;
	vector<string> unames,names;// unique population names
	if(npop==-1) npop=getDim(); //-1 means the number of individuals
	else if(npop==-2){ // -2 means try to infer from labels
		for(long i=0;i<getDim();i++) {
			names.push_back(data->getnames(i));
			names.back()=removeNumbers(names.back());
		}
		unames=names; // get the vector of unique names
 		sort (unames.begin(),unames.end() , less<string>());
		vector<string>::iterator it=unique (unames.begin(),unames.end() );
		unames.resize( it - unames.begin() );
		npop=unames.size();
	}
	if(npop<=0) throw(std::string("Bad initial population!"));
	indinp.resize(npop,vector<int>());
	psize.resize(npop,0);
	for(int i=0;i<getDim();i++) {
		if(getDim()==npop) ind.push_back(i);
		else if(unames.size()>0){
		  for(unsigned long j=0;j<unames.size();j++) if(names[i]==unames[j]){
			ind.push_back(j);
			j=unames.size();
		  }
		}else ind.push_back(RandomInteger(0,npop-1));
		psize[ind[i]]+=data->nindiv(i);
		indinp[ind[i]].push_back(i);
	}
	for(unsigned int i=0;i<psize.size();i++) if(psize[i]==0) removePop(i);
	createParams(a);
}

State::State(Data *d,FsXml *infile,vector<double> bvec,double a,int betamodel,double corfactor,bool readp,Data *d2,int datainference,int modeltype)
{
	this->modeltype=modeltype;
	dlength=d2;
	this->datainference=datainference;
	this->betamodel=betamodel;
	this->corfactor=corfactor;
	setdataused();
	hyperprior=bvec;
	data=d;
	alpha=a;
	readState(getFromFile(infile,readp));
	createParams(getAlpha());
	//cout<<setprecision (12)<< "Approx log-likelihood="<<getApproxLogLikelihood()<<endl;
}

State::State(Data *d,std::string stringin,vector<double> bvec,double a,int betamodel,bool isfile,double corfactor,Data *d2,int datainference,int modeltype)
{
	dlength=d2;
	this->modeltype=modeltype;	
	this->datainference=datainference;
	this->betamodel=betamodel;
	this->corfactor=corfactor;
	setdataused();
	hyperprior=bvec;
	data=d;
	string statestr;
	size_t found;
	if(isfile){// its a file but not a valid state file
	  ifstream file;
   	  file.open(stringin.data());//Open file
 	  if(!file.good()) throw(std::string("Invalid file!"));
	  string curline;
/*	  while (1){
		getline(file,curline);//Read next line from file
		found=curline.find_first_of('<');
		if(curline.length()>5+found) if(curline.substr(found,5).compare("<Pop>")==0) statestr=curline;
		if (file.eof()) break;
	  }
	  file.close();//Close file
	  if(statestr.length()==0) throw(std::string("Population tag not found!"));*/
	  while (1){
		getline(file,curline);//Read next line from file
		found=curline.find_first_of('(');
		if(found!=string::npos) statestr=curline;
		if (file.eof()) break;
	  }
	  file.close();//Close file
	  if(statestr.length()==0) throw(std::string("Population tag not found!"));
//	  cout<<"Using state "<<statestr.c_str()<<endl;
	}else statestr=stringin;
	readState(statestr);
	createParams(a);
}

State::State(State *in)
{
	copyState(in);
}

State::State(Data *d,std::vector< std::vector<double> > *vecin,bool mergerule, double mergeval,std::vector<double> bvec,double a,int betamodel,double corfactor,Data *d2,int datainference,int modeltype)
{
  	this->modeltype=modeltype;
	dlength=d2;
	this->datainference=datainference;
	this->betamodel=betamodel;
	this->corfactor=corfactor;
	setdataused();
	hyperprior=bvec;
	data=d;
// create a state with all individuals separate
	indinp.resize(getDim(),vector<int>());
	psize.resize(getDim(),0);
	for(int i=0;i<getDim();i++) {
		ind.push_back(i);
		psize[ind[i]]+=data->nindiv(i);
		indinp[ind[i]].push_back(i);
	}
// and apply the merging algorithm to them according to the pairwise coincidences
	createParams(a);
	mergeOnCoincidence(vecin,mergerule,mergeval);
}

std::string State::getFromFile(FsXml *infile,bool readparams)
{
	string ret;
	std::streampos sp=infile->tellg();
	while(1) {
		string res=infile->getLine();
		if(res.find("<Pop>")!=string::npos){
			ret=res;
		}
		if(readparams){
		  if(res.find("<beta>")!=string::npos)setBetaFromString(res);
		  if(res.find("<alpha>")!=string::npos)setAlphaFromString(res);
		  if(res.find("<delta>")!=string::npos)setDeltaFromString(res);
		  if(res.find("<F>")!=string::npos)setFFromString(res);
		}
		if (infile->eof() || res.find("</Iteration>")!=string::npos) break;
	}
	infile->seekg(sp);
	return(ret);
}

void State::mergeOnCoincidence(std::vector< std::vector<double> > *vecin,bool mergerule, double mergeval)
{
	if(mergerule) {
	  for(unsigned int c1=0;c1<vecin->size();c1++) {
		for(unsigned int c2=c1+1;c2<vecin->at(c1).size();c2++) {
		  if(vecin->at(c1)[c2] > mergeval && getPop(c1)!=getPop(c2)) {
			merge(getPop(c1),getPop(c2));}
		}
	  }
	}
}

void State::createParams(double a){
	mergepreva=-1;mergeprevb=-1;mergenewab=-1;
	alpha=a;
	sumbeta=0.0;
	beta=vector<double>();
	if(betamodel==BETAMOD_F || betamodel==BETAMOD_F2 || betamodel==BETAMOD_F2_COPYMAT) {
	  if(hyperprior.size()<4) throw(string("Param creation error: beta has wrong size for model!"));
	  if(hyperprior[0]<=-1) throw(string("Param creation error: fixed beta must be in the (0,1) range!"));
	  if((hyperprior[0]>0 && hyperprior[2]<0) || (hyperprior[1]>0 && hyperprior[3]<0)) throw(string("Param creation error: hyperpriors for beta must all be positive."));
	}

	if(betamodel==BETAMOD_COPYMAT || betamodel==BETAMOD_F2_COPYMAT){// only initialise once
		int offset=0;
		if(betamodel==BETAMOD_F2_COPYMAT)offset=4;
		  copyPrior=vector< vector<double> >(data->getDim(),vector<double>(data->getDim(),0));
		if((int)hyperprior.size()<data->getDim()*data->getDim()){cerr<<"matrix too small:"<<hyperprior.size()<<"<"<<data->getDim()*data->getDim()<<endl;throw(string("Invalid beta matrix"));}
		  for(int c1=0;c1<data->getDim();c1++) {
		    for(int c2=0;c2<data->getDim();c2++) {
			copyPrior[c1][c2]=hyperprior[c1*data->getDim()+c2+offset];
		    }
		  }
	}
	for(unsigned int i=0;i<psize.size();i++) {// one entry per population in these priors
		if(betamodel==BETAMOD_F){
		  if(hyperprior[0]>0) betaF.push_back(sampleBetaP(0,0)); else betaF.push_back(-hyperprior[0]);
		  if(hyperprior[1]>0) delta.push_back(sampleBetaP(1,0)); else delta.push_back(-hyperprior[1]);
		}else if(betamodel==BETAMOD_F2 || betamodel==BETAMOD_F2_COPYMAT) {
		  if(i==0 && delta.size()==0) {
			betaF.push_back(sampleBetaP(0,0));///**** UNDESIREABLE TO DO THIS WHEN GENERATING TREES
		  	delta.push_back(sampleBetaP(1,0));
		  }
		}else{
		  beta.push_back(fabs(hyperprior[0]));
		  sumbeta+=fabs(hyperprior[0]);
		}
	}
	popX=vector<vector<double> >();
	for(unsigned int a=0;a<psize.size();a++) {
		popX.push_back(vector<double>());
		for(unsigned int b=0;b<psize.size();b++) {
			popX.back().push_back(calcSumXab(a,b));
		}
	}
	if(uselengths) {
		popL=vector<vector<double> >();
		poplogLpi=vector<vector<double> >();
		poplogLgamma=vector<vector<double> >();
		for(unsigned int a=0;a<psize.size();a++) {
			popL.push_back(vector<double>());
			poplogLpi.push_back(vector<double>());
			poplogLgamma.push_back(vector<double>());
			for(unsigned int b=0;b<psize.size();b++) {
				popL.back().push_back(calcSumLab(a,b));
				poplogLpi.back().push_back(calcSumlogLpi(a,b));
				poplogLgamma.back().push_back(calcSumlogLgamma(a,b));
			}
		}
		vector<double> hyperlen=getLengthHyperPrior();
		for(int c1=0;c1<NUMHYPERPARAMLENGTH;c1++) if(hyperlen[c1]>0) priorLengths.push_back(sampleBetaP(c1,1));else priorLengths.push_back(fabs(hyperlen[c1]));
	}
	if(usesums) {
		muvec=vector<double>(data->getDim(),0);
		for(unsigned int c1=0;c1<muvec.size();c1++) muvec[c1]=calcMu(c1);
		vector<double> hyperlen2=getTotalChunksHyperPrior();
		for(int c1=0;c1<NUMHYPERPARAMTOTLENGTH;c1++){
		 if(hyperlen2[c1]>0) priorMeanSizes.push_back(sampleBetaP(c1,2)); else priorMeanSizes.push_back(fabs(hyperlen2[c1]));
		}
	}
	verbose=false;
}

void State::readState(std::string statein)
{
	indinp=vector<vector<int> >();
	psize=vector<int>();
	ind=vector<int>(getDim(),-1);
	int popon=-1;
	size_t found=statein.find_first_of("(),"), pos=0;
	while(found!=string::npos){
	  if(statein.at(pos)=='(') {
	    popon++;
	    indinp.push_back(vector<int>());
	    psize.push_back(0);
	    if(found==0) found=statein.find_first_of("(),",pos+1);
	  }
	  if(statein.at(pos)=='(' || statein.at(pos)==',') {
		int indiv=data->getIndex(statein.substr(pos+1,found-pos-1));
		if(indiv<0) {
		  cerr<<"Individual "<<statein.substr(pos+1,found-pos-1).c_str()<<" not recognised!"<<endl;
		  cerr<<"Options are: ";
		  for(int c1=0;c1<getDim();c1++)cerr<<data->getnames(c1)<<" ";
		  cerr<<endl;
		  throw(std::string("State creation error: name not found!"));
		}
		ind[indiv]=popon;
		psize.back()+=data->nindiv(indiv);
		indinp.back().push_back(indiv);
	  }
	  pos=found;
	  found=statein.find_first_of("(),",pos+1);
	}
	for(int i=0;i<getDim();i++) if(ind[i]<0) {
		cerr<<"Error: Individual "<<i<<" named "<<data->getnames(i)<<" not found in state input string."<<endl;
		throw(std::string("State creation error"));
	}
}

void State::copyState(State *in)
{
	data=in->data;
	psize=in->psize;
	ind=in->ind;
	dlength=in->dlength;
	usecounts=in->usecounts;
	uselengths=in->uselengths;
	usesums=in->usesums;

	datainference=in->datainference;
	alpha=in->alpha;
	hyperprior=in->hyperprior;
	beta=in->beta;
	betaF=in->betaF;
	delta=in->delta;
	betamodel=in->betamodel;
	corfactor=in->corfactor;
	sumbeta=in->sumbeta;

	priorLengths=in->priorLengths;
	verbose=in->verbose;
	popX=in->popX;
	popL=in->popL;
	poplogLpi=in->poplogLpi;
	poplogLgamma=in->poplogLgamma;
	copyPrior=in->copyPrior;
	for(unsigned int a=0;a<popX.size();a++) popX[a]=in->popX[a];
	if(uselengths) {
		for(unsigned int a=0;a<popL.size();a++) {
			popL[a]=in->popL[a];
			poplogLpi[a]=in->poplogLpi[a];
			poplogLgamma[a]=in->poplogLgamma[a];
		}
	}
	indinp=in->indinp;
	for(unsigned int a=0;a<indinp.size();a++) indinp[a]=in->indinp[a];
	diagmod=in->diagmod;
	epopsize=in->epopsize;
	muvec=in->muvec;
	priorMeanSizes=in->priorMeanSizes;
	mergepreva=in->mergepreva;mergeprevb=in->mergeprevb;mergenewab=in->mergenewab;
	modeltype=in->modeltype;	
}

void State::removePop(int lose,int keep)
{
  // update super individual properties
	if(keep>=0 && diagmod.size()>0 && epopsize.size()>0) {
		diagmod[keep]+=diagmod[lose];
//		epopsize[keep]=(2.0*psize[lose]*psize[keep])/(psize[lose]+psize[keep]);
		epopsize[keep]=(2.0*psize[lose]*psize[keep])/(psize[lose]+psize[keep]);
		//epopsize[keep]=sumXab(keep,keep) + sumXab(lose,lose);//sumXab(lose,keep) + sumXab(keep,lose)
		epopsize.erase(epopsize.begin() + lose);
		diagmod.erase(diagmod.begin() + lose);
	}
// update the individuals and populations
	if(keep>=0) {
		for(unsigned int a=0;a<indinp[lose].size();a++) indinp[keep].push_back(indinp[lose][a]);
		psize[keep]+=psize[lose];
	}
	indinp.erase(indinp.begin() + lose);
	for(unsigned int i=0;i<ind.size();i++) {
		if(ind[i]==lose && keep>=0) ind[i]=keep;
		else if(ind[i]==lose && keep<0) throw(std::string("No replacement population supplied!"));
		else if(ind[i]>lose) ind[i]--;
	}
// update the population copy counts
	if(keep>=0) {
		for(unsigned int a=0;a<psize.size();a++) {
			popX[a][keep]+=popX[a][lose];
			popX[keep][a]+=popX[lose][a];
		}
		popX[keep][keep]+=popX[lose][lose];
		if(uselengths){
			for(unsigned int a=0;a<psize.size();a++) {
				popL[a][keep]+=popL[a][lose];
				popL[keep][a]+=popL[lose][a];
				poplogLpi[a][keep]+=poplogLpi[a][lose];
				poplogLpi[keep][a]+=poplogLpi[lose][a];
				poplogLgamma[a][keep]+=poplogLgamma[a][lose];
				poplogLgamma[keep][a]+=poplogLgamma[lose][a];
			}
			popL[keep][keep]+=popL[lose][lose];
			poplogLpi[keep][keep]+=poplogLpi[lose][lose];
			poplogLgamma[keep][keep]+=poplogLgamma[lose][lose];
		}
	}
	for(unsigned int a=0;a<psize.size();a++) {popX[a].erase(popX[a].begin() + lose);}
	popX.erase(popX.begin() + lose);
	if(uselengths){
		for(unsigned int a=0;a<psize.size();a++) {
			popL[a].erase(popL[a].begin() + lose);
			poplogLpi[a].erase(poplogLpi[a].begin() + lose);
			poplogLgamma[a].erase(poplogLgamma[a].begin() + lose);
		}
		popL.erase(popL.begin() + lose);
		poplogLpi.erase(poplogLpi.begin() + lose);
		poplogLgamma.erase(poplogLgamma.begin() + lose);
	}

// update the other vectors
	psize.erase(psize.begin() + lose);
	if(betamodel==BETAMOD_F){
	  betaF.erase(betaF.begin()+lose);
	  delta.erase(delta.begin()+lose);
	}else if(betamodel==BETAMOD_F2_COPYMAT){
	  // nothing to do
	}else if(betamodel!=BETAMOD_F2 && betamodel!=BETAMOD_F2_COPYMAT){
	  sumbeta-=beta[lose];
	  beta.erase(beta.begin()+lose);
	}

}

void State::moveInd(int i,int popto){
	int popfrom=getPop(i);

  	popX[popto][popto]+=X(i,i);// self copying gets counted twice
	popX[popfrom][popfrom]+=X(i,i);
	popX[popto][popfrom]-=X(i,i);// self copying gets counted twice
	popX[popfrom][popto]-=X(i,i);
 	for(unsigned int a=0;a<psize.size();a++) {
		double sxa=sumX(i,a),say=sumY(a,i);
		popX[a][popto]+=say;
		popX[popto][a]+=sxa;
		popX[a][popfrom]-=say;
		popX[popfrom][a]-=sxa;
		if(uselengths) {
			double sxa=sumLx(i,a),say=sumLy(a,i);
			popL[a][popto]+=say;
			popL[popto][a]+=sxa;
			popL[a][popfrom]-=say;
			popL[popfrom][a]-=sxa;
			sxa=sumlogLpix(i,a);say=sumlogLpiy(a,i);
			poplogLpi[a][popto]+=say;
			poplogLpi[popto][a]+=sxa;
			poplogLpi[a][popfrom]-=say;
			poplogLpi[popfrom][a]-=sxa;
			sxa=sumlogLgammax(i,a);say=sumlogLgammay(a,i);
			poplogLgamma[a][popto]+=say;
			poplogLgamma[popto][a]+=sxa;
			poplogLgamma[a][popfrom]-=say;
			poplogLgamma[popfrom][a]-=sxa;
		}
	}
// move the individual between the two population lists
	int found=-1;
	for(unsigned int a=0;a<indinp[popfrom].size();a++){
		if(indinp[popfrom][a]==i) {found=a;break;}
	}
	if(found<0) {
		cerr<<"ERROR: population missing its individual!"<<endl;
		setprint(&cout);
		for(unsigned int a=0;a<indinp.size();a++) {
			cout<<"Pop "<<a<<": ";
			for(unsigned int b=0;b<indinp[a].size();b++) cout<<indinp[a][b]<<",";
			cout<<endl;
		}
		throw(std::string("Individual not found in its population!"));
	}
	indinp[popto].push_back(i);
	indinp[popfrom].erase(indinp[popfrom].begin()+found);
// update the population sizes and individuals record of its pop
	psize[popto]+=data->nindiv(i);
	psize[popfrom]-=data->nindiv(i);
	ind[i]=popto;
}

void State::setIntendedSplit(int i, int j){
	mergepreva=i;
	mergeprevb=j;
	mergenewab=mergepreva; // this population is kept
};
    
void State::merge(int a,int b)
{
  if(a==b) {cerr<<"Error in state:merge - merging identical populations!"<<endl;throw(string("Merging identical populations!"));}
	int keep=min(a,b);
	int lose=max(a,b);
	vector<int> tmp=getIndInPop(keep);
	if(tmp.size()>0) mergepreva=tmp[0];
	tmp=getIndInPop(lose);
	if(tmp.size()>0) mergeprevb=tmp[0];
	mergenewab=mergepreva; // this population is kept

/*      cout<<"Merging Pop "<<ind[mergepreva]<<":";
       vector<int> indina=getIndInPop(ind[mergepreva]);
      for(int i =0;i<indina.size();i++) cout<<indina[i]<<",";
      cout<<endl;
	cout<<"With Pop "<<ind[mergeprevb]<<":";      
       vector<int> indinb=getIndInPop(ind[mergeprevb]);
      for(int i =0;i<indinb.size();i++) cout<<indinb[i]<<",";
      cout<<endl;*/
	removePop(lose,keep);
/*	cout<<"To make Pop "<<ind[mergenewab]<<":";      
       vector<int> indinab=getIndInPop(ind[mergenewab]);
      for(int i =0;i<indinab.size();i++) cout<<indinab[i]<<",";
      cout<<endl;	*/
}

int State::addEmptyPop()
{
// add beta
	if(betamodel==BETAMOD_F){
	  betaF.push_back(sampleBetaP(0,0)); //***
	  delta.push_back(sampleBetaP(1,0)); //***
	}else if(betamodel!=BETAMOD_F2){
	  double avep=0.0;
	  for(unsigned int i=0;i<beta.size();i++) avep+=beta[i]/(double)beta.size();
	  beta.push_back(avep);
	  sumbeta+=avep;
	}
// add to popX
	for(unsigned int a=0;a<psize.size();a++) {
		popX[a].push_back(0);
	}
	popX.push_back(vector<double>(psize.size()+1,0));
	if(uselengths){
		for(unsigned int a=0;a<psize.size();a++) {
			popL[a].push_back(0);
			poplogLpi[a].push_back(0);
			poplogLgamma[a].push_back(0);
		}
		popL.push_back(vector<double>(psize.size()+1,0));
		poplogLpi.push_back(vector<double>(psize.size()+1,0));
		poplogLgamma.push_back(vector<double>(psize.size()+1,0));
	}
// add to popsize
	indinp.push_back(vector<int>());
	psize.push_back(0);
	return(psize.size()-1);// returns the index of the new pop
}

double State::probOfSplitSAMS(int i, int j, State * relstate)
{
//  verbose=true;
	bool valid=true;
	int a=ind[i],b=ind[j];// get the initial populations
	if(verbose) cout<<"Merging pops of individuals "<<i<<","<<j<<" in populations "<<a<<","<<b<<"; testing the split probability."<<endl;
	double logpofsplit=0.0;
	// move i and j to new populations
	moveInd(i,addEmptyPop());
	moveInd(j,addEmptyPop());
// apply the PCA automerge
	//if(pca!=NULL) valid=testPcaAuto(a,b,i,j,pca->getDist(i,j),pca);
	//if(logpofsplit<0) return(logpofsplit); // -INFINITY, i.e. impossible to get the split this way
	// but can't return yet since we need a completed merge to get a valid state
	// later will implement forcing this merge...
	
// List the other individuals
	vector<int> olda=getIndInPop(a);
	vector<int> oldb=getIndInPop(b);
	vector<int> opop=olda;
	opop.insert (opop.end(),oldb.begin(),oldb.end());
	vector<int> frompop;
	for(unsigned int c1=0;c1<olda.size();c1++) frompop.push_back(0);
	for(unsigned int c1=0;c1<oldb.size();c1++) frompop.push_back(1);
// permute the list
	rpermute(&opop,&frompop);
// do the merge
	merge(a,b);// merge the populations
	if(getPlength(a)==0) removePop(a);
	if(getPlength(b)==0) removePop(b);
	if(opop.size()==0) return(0.0);// there was only two individuals
	int opindex=ind[opop[0]]; //original population index
	a=ind[i],b=ind[j];
// Do the split
	for(unsigned int k=0;k<opop.size();k++) {
  // get the probabilities of the current state
		  double oldlptoi=posteriorSetProbPartial(ind[i],relstate);
		  double oldlptoj=posteriorSetProbPartial(ind[j],relstate);
		  double oldlp0=posteriorSetProbPartial(opindex,relstate);
  // move individual and get new probabilities
		  moveInd(opop[k],ind[i]);
		  double lptoi=posteriorSetProbPartial(ind[i],relstate)+posteriorSetProbPartial(opindex,relstate)-oldlptoi-oldlp0;
		  double psizetoi=(double)psize[ind[i]];
		  //if(verbose) {cout<<"trying to move "<<data->getnames(opop[k])<<" to pop "<<ind[i]<<": " <<psizetoi<<" + log "<<lptoi <<endl;}
  // move individual to the other population and get probabilities
		  moveInd(opop[k],ind[j]);
		  double lptoj=posteriorSetProbPartial(ind[j],relstate)+posteriorSetProbPartial(opindex,relstate)-oldlptoj-oldlp0;
		  double psizetoj=(double)psize[ind[j]];
		  //if(verbose) {cout<<"trying to move "<<data->getnames(opop[k])<<" to pop "<<ind[j]<<": " <<psizetoj<<" + log "<<lptoj <<endl;}
  // evaluate the probability of each move.
  // NOTE: r < Aa/(Aa+Bb) => 1/r > 1 + Bb/Aa = 1 + (b/a)*exp(log(B)-log(A))
  /////////////////////////////****************
		double lpi=-log(1.0 +(psizetoj/psizetoi)*exp(lptoj-lptoi));
  		double lpj=-log(1.0 +(psizetoi/psizetoj)*exp(lptoi-lptoj));
//		  lpi=-log(0.5);
//		  lpj=-log(0.5);
		if(verbose) cout<<"trying to move "<<data->getnames(opop[k])<<" to pop "<<ind[i]<<": " <<psizetoi<<" + "<<lptoi <<endl;

		if(frompop[k]==0) {
			if(verbose) cout<<"Moved to "<<ind[i]<<endl;
			logpofsplit+=lpi;
			moveInd(opop[k],ind[i]);
		}else {
			if(verbose) cout<<"Moved to "<<ind[j]<<endl;
			logpofsplit+=lpj;
			moveInd(opop[k],ind[j]);
		}
	}
	removePop(opindex);
	if(valid) return(logpofsplit);
	return(-DBL_MAX);
}

double State::splitSAMS(int i, int j,bool greedy, State * relstate)
{

	if(i==j) {cerr<<"WARNING: Tried to split a single individual!"<<endl;return(0);}
	if(ind[i]!=ind[j]) {cerr<<"WARNING: Tried to split individuals in different populations!"<<endl;return(0);}
	int opindex=ind[i]; //original population index
	// move i and j to new populations
	vector<int> opop=getIndInPop(opindex);
	if(verbose) {cout<<"Splitting pop "<<opindex<<" containing:";
	  for(unsigned int c1=0;c1<opop.size();c1++) cout<<data->getnames(opop[c1])<<",";
	  cout<<endl; setprint(&cout);}
// move and create new populations
	moveInd(i,addEmptyPop());
	moveInd(j,addEmptyPop());
// permute the rest of the individuals
	//if(pca!=NULL) movePcaAuto(opindex,i,j,pca->getDist(i,j),pca);
	opop=getIndInPop(opindex);
	if(opop.size()>0) rpermute(&opop);
	if(verbose) cout<<"SEEDING WITH "<<data->getnames(i)<<" and "<<data->getnames(j)<<endl;
	double logpofsplit=0.0;
	for(unsigned int k=0;k<opop.size();k++) {
		double lpi=0;
		double lpj=0;
 // obtain the probabilities of the initial populations
		  double oldlptoi=posteriorSetProbPartial(ind[i],relstate);
		  double oldlptoj=posteriorSetProbPartial(ind[j],relstate);
		  double oldlp0=posteriorSetProbPartial(opindex,relstate);
  // move and get new probabilities
		  moveInd(opop[k],ind[i]);
		  double lptoi=posteriorSetProbPartial(ind[i],relstate)-oldlptoi+posteriorSetProbPartial(opindex,relstate)-oldlp0;
		  double psizetoi=(double)psize[ind[i]];
		  if(verbose) cout<<"trying to move "<<data->getnames(opop[k])<<" to pop "<<ind[i]<<": " <<psizetoi<<" + log "<<lptoi <<endl;
  // move to other population and get new probabilities
		  moveInd(opop[k],ind[j]);
		  double lptoj=posteriorSetProbPartial(ind[j],relstate)-oldlptoj+posteriorSetProbPartial(opindex,relstate)-oldlp0;
		  double psizetoj=(double)psize[ind[j]];
		  //if(verbose) cout<<"trying to move "<<data->getnames(opop[k])<<" to pop "<<ind[j]<<": " <<psizetoj<<" + "<<lptoj <<endl;
  // calculate the probabilities
  // NOTE: r < Aa/(Aa+Bb) => 1/r > 1 + Bb/Aa = 1 + (b/a)*exp(log(B)-log(A))
  /////////////////////////////****************
  		lpi=-log(1.0 +(psizetoj/psizetoi)*exp(lptoj-lptoi));
  		lpj=-log(1.0 +(psizetoi/psizetoj)*exp(lptoi-lptoj));
//		  lpi=-log(0.5);
//		  lpj=-log(0.5);

		if(verbose) cout<<"Prob of pop "<<ind[i]<<"="<<lpi<<" and of pop "<<ind[j]<<"="<<lpj<<endl;
		if(greedy) {
			if(lpi>log(0.5000001) || (lpi>log(0.499999) &&rnd()<0.5)) {
				moveInd(opop[k],ind[i]);
				if(verbose) cout<<"Moved to "<<ind[i]<<endl;
			}else {
				moveInd(opop[k],ind[j]);
				if(verbose) cout<<"Moved to "<<ind[j]<<endl;
			}
		}else if(log(rnd()) < lpi) {
				if(verbose) cout<<"Moved to "<<ind[i]<<endl;
			logpofsplit+=lpi;
			moveInd(opop[k],ind[i]);
		}else {
				if(verbose) cout<<"Moved to "<<ind[j]<<endl;
			logpofsplit+=lpj;
			moveInd(opop[k],ind[j]);
		}
	};
	removePop(opindex);
	if(ind[i]<ind[j]){
	  mergepreva=i;mergeprevb=j;
	}else{
	  mergeprevb=i;mergepreva=j;
	}
	return(logpofsplit);
}

void State::splitSAMS(int a, State * relstate)
{
	vector<int> pop=getIndInPop(a);
	if(pop.size()==1) {cerr<<"WARNING: splitting population of size 1!"<<endl;return;}
	int i=RandomInteger(0,pop.size()-1);
	int j=i;
	while(j==i) j=RandomInteger(0,pop.size()-1);
	splitSAMS(pop[i],pop[j],false,relstate);
}


void State::splitSAMSgreedy(int a,int cmax, State * relstate)
{
	vector<int> pop=getIndInPop(a);
	if(pop.size()==1) {cerr<<"WARNING: splitting population of size 1!"<<endl;return;}
	if(cmax<0) throw(std::string("Need cmax to be set!"));
	State oldstate(this);
	State beststate(this);
	double splitbest=-10e30,cursplit;
	bool rand=true;
	if((int) (pop.size()*(pop.size()-1))<2*cmax) {cmax=pop.size();rand=false;}
	else cmax=(int)ceil(sqrt(cmax));
	int i,j;
	if(verbose) print(&cout,false);
// get a list of all the possible split seeds, or use random seeds
	for(int c1=0;c1<cmax;c1++) {
 	  for(int c2=c1+1;c2<cmax;c2++) {
		if(rand) {
		  i=RandomInteger(0,pop.size()-1);
		  j=i;
		  while(j==i) j=RandomInteger(0,pop.size()-1);
		}else{
		  i=c1;
		  j=c2;
		}
// do the splits
		splitSAMS(pop[i],pop[j],true,relstate);
	if(verbose) cout<<"Testing state with probability "<<posteriorProb(relstate)<<endl;
	if(verbose) print(&cout,false);
// test whether to keep the split
		cursplit=posteriorProb(relstate);
		if(cursplit<splitbest) {
		// do nothing! reject this split
		}else {// accept the split
		  splitbest=cursplit;
		  beststate.copyState(this);
		}
    	  	copyState(&oldstate);
	  }
	}
 	copyState(&beststate);
}

vector <int> State::getIndInPop(int pop)
{
	if(pop<0 || pop>=(int)indinp.size()) throw(string("Error in getindinp: invalid index!"));
//cout<<"pop="<<pop<<" indinsize="<<indinp.size()<<endl;
//cout<<"indinp[pop]="<<indinp[pop].size()<<endl;
	return(indinp[pop]);
}

void State::setIndInPop(std::vector<int> v,int i)
{
	if(v.size()!=indinp[i].size()) throw(string("Setting population to invalid value!"));
	vector<int> tmp=indinp[i];
	for(unsigned int c1=0;c1<v.size();c1++) {
	  bool found=false;
	  for(unsigned int c2=0;c2<tmp.size();c2++) {
		if(tmp[c2]==v[c1]) found=true;
	  }
	  if(found==false) throw(string("Invalid individual set in setindinpop!"));
	  indinp[i][indinp[i].size()-c1-1]=v[c1];
	}
}

double State::sumX(int i, int pop,bool cf)
{
	double s=0;
	vector<int> jlist=getIndInPop(pop);
	for(unsigned int j=0;j<jlist.size();j++) s+=X(i,jlist[j],cf);
	return(s);
}

double State::sumY(int pop, int i,bool cf)
{
	double s=0;
	vector<int> jlist=getIndInPop(pop);
	for(unsigned int j=0;j<jlist.size();j++) s+=X(jlist[j],i,cf);
	return(s);
}

double State::calcSumXab(int a, int b,bool cf)
{
	double s=0;
	vector<int> ilist=getIndInPop(a);
	for(unsigned int i=0;i<ilist.size();i++) s+=sumX(ilist[i],b,cf);
	return(s);
}

double State::sumXab(int a, int b,bool cf)
{
//	cout<<"a="<<a<<" b="<<b<<" stored="<<popX[a][b]<<" calc="<<calcSumXab(a,b)<<" psize=["<<getPsize(a)<<","<<getPsize(b)<<"] plength=["<<getPlength(a)<<","<<getPlength(b)<<"]"<<endl;
//	return(calcSumXab(a,b,cf));
//	if(popX[a][b]<0) {if(popX[a][b]>-0.00001) {popX[a][b]=0;}else{cout<<"popx of "<<a<<" "<<b<< "is "<<popX[a][b]<<endl;exit(0);}}
	if (cf) return(popX[a][b]);
	return(popX[a][b]*corfactor);
}

double State::sumLx(int i, int pop)
{
	double s=0;
	vector<int> jlist=getIndInPop(pop);
	for(unsigned int j=0;j<jlist.size();j++) s+=L(i,jlist[j]);
	return(s);
}

double State::sumLy(int pop, int i)
{
	double s=0;
	vector<int> jlist=getIndInPop(pop);
	for(unsigned int j=0;j<jlist.size();j++) s+=L(jlist[j],i);
	return(s);
}

double State::sumlogLpix(int i, int pop)
{
	double s=0;
	vector<int> jlist=getIndInPop(pop);
	for(unsigned int j=0;j<jlist.size();j++) s+=logLpi(i,jlist[j]);
	return(s);
}

double State::sumlogLpiy(int pop, int i)
{
	double s=0;
	vector<int> jlist=getIndInPop(pop);
	for(unsigned int j=0;j<jlist.size();j++) s+=logLpi(jlist[j],i);
	return(s);
}

double State::sumlogLgammax(int i, int pop)
{
	double s=0;
	vector<int> jlist=getIndInPop(pop);
	for(unsigned int j=0;j<jlist.size();j++) if(data->get(i,jlist[j])!=0) s+=mylgamma(data->get(i,jlist[j]));
	return(s);
}

double State::sumlogLgammay(int pop, int i)
{
	double s=0;
	vector<int> jlist=getIndInPop(pop);
	for(unsigned int j=0;j<jlist.size();j++) if(data->get(i,jlist[j])!=0) s+=mylgamma(data->get(jlist[j],i));
	return(s);
}

void State::printX(ostream * out,bool perindiv)
{
	if(perindiv) *out<<"relXmat,";
	else  *out<<"Xmat,";
	double denom;
	for(unsigned int i=0;i<psize.size()-1;i++) *out<<i<<", ";
	*out<<psize.size()-1<<endl;
	for(unsigned int i=0;i<psize.size();i++) {
	  *out<<i<<", ";
	  if(perindiv) {
		double partdenom=sumXa(i);
		for(unsigned int j=0;j<psize.size()-1;j++) {
//			if(i==j)denom=(psize[i])*(psize[i]-1);
//			else denom=psize[i]*psize[j];
			if(i==j)denom=(partdenom)*(psize[i]-1);
			else denom=partdenom*psize[j];
//			if(i==j)denom=(psize[j]-1);
//			else denom=psize[j];
			if(denom==0) *out<<"0.0,";
			else *out<<sumXab(i,j)/denom<<",";
		}
		unsigned int j=psize.size()-1;
//		if(i==j)denom=(psize[j]-1)*psize[i]-1;
//		else denom=psize[i]*psize[j];
		if(i==j)denom=(partdenom)*(psize[i]-1);
		else denom=partdenom*psize[j];
//			if(i==j)denom=(psize[j]-1);
//			else denom=psize[j];
		if(denom==0)*out<<"0.0"<<endl;
		else *out<<sumXab(i,j)/denom<<endl;
	  }else{
		for(unsigned int j=0;j<psize.size()-1;j++) *out<<sumXab(i,j)<<",";
	  	*out<<sumXab(i,psize.size()-1)<<endl;
	  }
	}
//	*out<<"</X>"<<endl;
}


void State::printBeta(ostream * out)
{
	*out<<"beta,";
	for(unsigned int i=0;i<psize.size()-1;i++) *out<<i<<", ";
	*out<<psize.size()-1<<endl;
	for(unsigned int i=0;i<psize.size();i++) {
	  *out<<i<<", ";
	  for(unsigned int j=0;j<psize.size()-1;j++) *out<<getBeta(i,j)<<",";
	  *out<<getBeta(i,psize.size()-1)<<endl;
	}
//	*out<<"</beta>"<<endl;
}

void State::setSuperIndivRule(bool super){
	if(super){// set up diagmod and epopsize
	  diagmod=std::vector<double>(getP(),0);
	  epopsize=std::vector<double>(getP(),0);
	  for(unsigned int a=0;a<psize.size();a++) {
	    vector<int> popa=getIndInPop(a);
	    double maxa=0;
	    for(unsigned int b=0;b<psize.size();b++) {if(b!=a){
		vector<int> popb=getIndInPop(b);
		for(unsigned int c1=0;c1<popa.size();c1++)for(unsigned int c2=0;c2<popb.size();c2++) if(data->get(popa[c1],popb[c2])/data->nindiv(popa[c1])/data->nindiv(popb[c2]) > maxa) maxa=data->get(popa[c1],popb[c2])/data->nindiv(popa[c1])/data->nindiv(popb[c2]);
	    }}
	    for(unsigned int c1=0;c1<popa.size();c1++){for(unsigned int c2=0;c2<popa.size();c2++){ if(data->get(popa[c1],popa[c2])/data->nindiv(popa[c1])/data->nindiv(popa[c2])> maxa){
		//double diff=data->get(popa[c1],popa[c2]) - maxa;
		//cout<<setprecision(10)<<"Flattening "<<popa[c1]<<","<<popa[c2]<<" from "<<data->get(popa[c1],popa[c2])/data->nindiv(popa[c1])/data->nindiv(popa[c2])<<" to "<<maxa<<endl;
		data->set(popa[c1],popa[c2],maxa*data->nindiv(popa[c1])*data->nindiv(popa[c2]));
		popX[a][a]=calcSumXab(a,a);
	    }}}
	  }
	}else {diagmod=std::vector<double>(0);epopsize=std::vector<double>(0);}
}   

void State::print(ostream * out,bool chat)
{
	if(chat) *out<<"Printing State"<<endl<<"IND:, ";
	for(unsigned int i=0;i<ind.size()-1;i++) *out<<ind[i]<<", ";
	*out<<ind[ind.size()-1]<<endl;
	if(chat) *out<<"POP:"<<endl;
	for(unsigned int i=0;i<psize.size();i++) {
	  *out<<"pop "<<i<<" size="<<psize[i]<<": ";
	  vector<int> tp=getIndInPop(i);
	  for(unsigned int j=0;j<tp.size()-1;j++) *out<<data->getnames(tp[j])<<",";
	  *out<<data->getnames(tp[tp.size()-1])<<endl;
	}

	if(chat) {
	*out<<"XMAT (to)"<<endl;
	for(int i=0;i<getDim();i++) {
	  *out<<data->getnames(i)<<", ";
	  for(unsigned int j=0;j<psize.size()-1;j++) {
	    *out<<sumX(i,j)<<", ";
	  }
	  *out<<sumX(i,psize.size()-1)<<endl;
	}
	*out<<"YMAT (from)"<<endl;
	for(int i=0;i<getDim();i++) {
	  *out<<data->getnames(i)<<", ";
	  for(unsigned int j=0;j<psize.size()-1;j++) {
	    *out<<sumY(j,i)<<", ";
	  }
	  *out<<sumY(psize.size()-1,i)<<endl;
	}
	}

}

void State::setprint(ostream * out)
{
	*out<<"<Pop>"<<flush;
	for(int popon=0;popon<getP();popon++) {
	  vector<int> pop=getIndInPop(popon);
	  if(pop.size()>0){
	  *out<<"("<<flush;
	  for(unsigned int i=0;i<pop.size()-1;i++) {
		*out<<data->getnames(pop[i])<<","<<flush;
	  }
		*out<<data->getnames(pop[pop.size()-1])<<")"<<flush;
	  }
	}
	*out<<"</Pop>"<<endl;
}

double State::posteriorProb(State * relstate)
{
	double logp=psize.size() * log(alpha);
	for(unsigned int a=0;a<psize.size();a++) {
		logp+=mylgamma(getPsize(a));
		if(usecounts && (modeltype==MODELTYPE_NORMALISED || modeltype==MODELTYPE_NORMALISEDMERGEONLY)) logp+=posteriorSetProb(a,relstate);
		else if(usecounts && modeltype==MODELTYPE_INDIVIDUAL) logp+=posteriorIndSetProb(a);
		else if(usecounts)logp+=posteriorSetProb(a); // Actually we do a lot of stuff in infmcmc that we don't need to if we use this model
		//if(usecounts) logp+=posteriorSetProb(a);
		if(uselengths) logp+=posteriorLengthSetProb(a);
		if(usesums) logp+=posteriorNumChunksSetProb(a);
	}
	return(logp);
}

double State::posteriorSetProb(int a)
{
      std::vector<int> tmp=getIndInPop(a);

        double logp=0.0;
        double sXa=sumXa(a);
// Normalising factor
        if(psize[a]==0) return(0);// catch when merges remove population
        double sbeta=getSumBeta(a);
        logp+=mylgamma(sbeta) - mylgamma(sXa+sbeta);
	if(checkPopForIgnoredSuper(a)){
	  return(-INFINITY);
	}
// product over other populations
        for(unsigned int b=0;b<psize.size();b++) {
          double sumab=sumXab(a,b);
          double betaab=getBeta(a,b);
          int epop=getPsize(b);
          if(a==(int)b) epop--;
          if(epop>0) {
                logp+=mylgamma(betaab+sumab) - mylgamma(betaab) - sumab*log(epop);
          }
        }
        return(logp);
}

	

double State::posteriorSetProb(int a,State * relstate)
{
  if(relstate==NULL) return(posteriorSetProb(a));
  if(relstate->mergepreva<0 || relstate->mergeprevb<0 || relstate->mergenewab<0) return(posteriorSetProb(a));
  //||getP()<=relstate->getP() was also in above if
  if(psize[a]==0) return(0);// catch when merges remove population

	double logp=0.0;
	//double sXa=sumXa(a);
	double sXa=0;
//	for(int b=0;b<psize.size();b++)sXa+=sumXab(a,b);
	double smod=0;
//	bool userelstate=false;
	//double sbeta=getSumBeta(a);
	double sbeta=0;
    /*if(relstate!=NULL) if(relstate->mergepreva>=0) if(a==ind[relstate->mergepreva] &&getP()>relstate->getP()){
      cout<<"PSP: a="<<a<<" newa="<<relstate->ind[relstate->mergenewab]<<" This K="<<getP()<<" compared to "<<relstate->getP()<<" mnewab="<<mergenewab<<endl;
      vector<int> indinab=relstate->getIndInPop(relstate->ind[relstate->mergenewab]);
      vector<int> indina=getIndInPop(a);
      for(int i =0;i<indina.size();i++) cout<<indina[i]<<",";
      cout<<endl;
      for(int i =0;i<indinab.size();i++)cout<<indinab[i]<<",";
      cout<<endl;
    }*/
	  
// product over other populations

	for(int b=0;b<(int)psize.size();b++) {
	  int epop=getPsize(b);
	  if(a==(int)b) epop--;
	  double epop2=getPsize(a);
	  if(a==b)epop2-=1.0;	  
	  if(epop<=0) continue;
	  if(epop2<=0) continue;
	  smod=0;
	  double sumab=sumXab(a,b);
	  double betaab=getBeta(a,b);
	    vector<int> indina=getIndInPop(a);
	    vector<int> indinb=getIndInPop(b);
	    int relpop=relstate->ind[relstate->mergenewab];
	    int relpopa=relstate->ind[indina[0]],relpopb=relstate->ind[indinb[0]];
	    //if(relpopa!=relpopb || relpopa!=relpop) continue;
	    double psizea=0,psizeb=0,psizeab=0;
	    vector<int> indinab=relstate->getIndInPop(relpop);
	      //if(data->nindiv(indina[0])>1 || data->nindiv(indinb[0])>1) break;
	    if(getPsize(a)==(int)indina.size() && getPsize(b)==(int)indinb.size()) {
	      if((a==ind[relstate->mergepreva] || a==ind[relstate->mergeprevb])&&
	       (b==ind[relstate->mergepreva] || b==ind[relstate->mergeprevb])) {
		// Both pops a and b were merged
		relpopa=relpop;
		relpopb=relpop;
		psizea=relstate->getPsize(relpopa)-1;
		psizeb=relstate->getPsize(relpopb)-1;
		psizeab=relstate->getPsize(relpopa)*(relstate->getPsize(relpopb)-1);
	      }else if((a==ind[relstate->mergepreva] || a==ind[relstate->mergeprevb])){
	      // Only pop a was merged
	      	relpopa=relpop;
		psizea=relstate->getPsize(relpopa);
		psizeb=relstate->getPsize(relpopb);
		psizeab=relstate->getPsize(relpopa)*(relstate->getPsize(relpopb));
	      }else if((b==ind[relstate->mergepreva] || b==ind[relstate->mergeprevb])){
	      // Only pop a was merged
	      	relpopb=relpop;
		psizea=relstate->getPsize(relpopa);
		psizeb=relstate->getPsize(relpopb);
		psizeab=relstate->getPsize(relpopa)*(relstate->getPsize(relpopb));
	      }
//	      if((a==ind[relstate->mergepreva] || a==ind[relstate->mergeprevb]) ||
//	          (b==ind[relstate->mergepreva] || b==ind[relstate->mergeprevb])){
		// calculate the sum modifier
		if(psizeab>0)smod+= getPsize(a)*epop* 2.0*relstate->sumXab(relpopa,relpopb)/psizeab;
		for(unsigned int i =0;i<indina.size();i++){
		  if(psizeb>0)smod-=epop*relstate->sumX(indina[i],relpopb)/psizeb;
		}
		for(unsigned int j =0;j<indinb.size();j++){
		  if(psizea>0)smod-=epop2*relstate->sumY(relpopa,indinb[j])/psizea;
		}
		if(sumab+smod<0) {
		    cout<<setprecision(9)<<"a="<<a<<" b="<<b<<" sumab="<<sumab<<" sumabmod="<<sumab+smod<<endl;
		  cerr<<"Error: sumab="<<sumab<<" but smod="<<smod<<endl;
		  cerr<<"psizea="<<psizea<<" psizeb="<<psizeb<<" psizeab="<<psizeab<<" epop="<<epop<<" epop2="<<epop2<<" indina.size()="<<indina.size()<<" indinb.size()="<<indinb.size()<<endl;
		  cerr<<"old pops="<<relpopa<<" "<<relpopa<<" rp="<<relpop<<endl;
		  cerr<<int(a==ind[relstate->mergepreva]) <<int(a==ind[relstate->mergeprevb])<<int(b==ind[relstate->mergepreva]) << int(b==ind[relstate->mergeprevb])<<endl;
		  cerr<<"! in pop "<<a<<" (";
		  for(unsigned int i =0;i<indina.size();i++)cerr<<data->getnames(indina[i])<<","<<flush;
		  cerr<<" size "<<getPsize(a)<<" | "<<flush;
		  for(unsigned int i =0;i<indinb.size();i++)cerr<<data->getnames(indinb[i])<<","<<flush;
		  cerr<<" size "<<getPsize(b)<<")"<<endl;
		  smod=-sumab;
		  throw(string("Error with Smod in PosteriorSetProb!"));
		}
//	      }
	    }// end matching psize
	   // cout<<setprecision(9)<<"a="<<a<<" b="<<b<<" sumab="<<sumab<<" sumabmod="<<sumab+smod<<endl;
	    // arbitrary penalty for merging continents
	    //smod=0;
//	  if( getPsize(a)!=(int)indina.size() && (int)indina.size()>1) logp-= 10000000000000000000.0;// efectively forbid merge events for continents
//	  if( getPsize(b)!=(int)indinb.size() && (int)indinb.size()>1) logp-= 10000000000000000000.0;
	  bool calcpost=false;
	  if(modeltype==MODELTYPE_NORMALISED) calcpost=true;
	  if(relstate->ind[indina[0]]==relpop && relstate->ind[indinb[0]]==relpop)calcpost=true;
	  if(calcpost) {
	    logp+=mylgamma(betaab+sumab+smod) - mylgamma(betaab) - (sumab+smod)*log(epop);
	//    cout<<setprecision(10)<<"a="<<a<<" b="<<b<<" logp="<<logp<<" sXa="<<sumab<<" Mod="<<sumab+smod<<endl;
	//    cout<<setprecision(10)<<"a2="<<ind[relstate->mergepreva]<<"b2="<<ind[relstate->mergeprevb]<<endl;
	    sXa+=sumab+smod;
	    sbeta+=betaab;
	  }
	}
	if(sbeta>0) logp+=mylgamma(sbeta) - mylgamma(sXa+sbeta);
	if(sbeta<=00 && sXa>0) cerr<<"sb="<<sumbeta<<" sxa="<<sXa<<endl;
	//return(logp + posteriorIndSetProb(a));
	//return(0.5*(logp + posteriorIndSetProb(a)));
	return(logp);
}


double State::posteriorIndSetProb(int a)
{//***
	if(psize[a]==0) return(0);// catch when merges remove population
	std::vector<int> tmp=getIndInPop(a);
	double logp=0.0;
	
// Normalising factor
	//cout<<getSumBeta(a)<<" + "<<sXa+getSumBeta(a)<<endl;
	for(unsigned int i=0;i<tmp.size();i++){
	  //cout<<"i="<<i<<endl;
	  double sXi=sumXall(tmp[i]);
	  double sbeta=getSumBeta(a)/psize[a];//1.0;//getSumBeta(a);
	  logp+=mylgamma(sbeta) - mylgamma(sXi+sbeta);
  // product over other populations
	  for(unsigned int b=0;b<psize.size();b++) {
	    //cout<<"b="<<b<<endl;
	    double sumib=sumX(tmp[i],b);
	    double betaib=getBeta(a,b)/psize[a];//1.0/psize.size();
	    int epop=getPsize(b);
	    if(a==(int)b) epop--;
	    if(epop>0) {
		  logp+=mylgamma(betaib+sumib) - mylgamma(betaib) - sumib*log(epop);
	    }

	  }
	}
	//cout<<"Popsize "<<tmp.size()<<" log posterior="<<logp<<endl;
	return(logp);
}


double State::posteriorSetProbPartial(int a,State * relstate)
{
	double logp=0;
	if(usecounts) logp+=posteriorSetProb(a,relstate);
	if(uselengths) logp+=posteriorLengthSetProb(a);
	if(usesums) logp+=posteriorNumChunksSetProb(a);
	return(logp);
///< This is all ignored; it was an approximation that proved not worth it.
/*	double logp=0.0;
	double sXa=sumXa(a);
// Normalising factor
	logp+=mylgamma(getSumBeta(a));
// product over other populations
	for(unsigned int b=0;b<psize.size();b++) {
	  double sXb=sumXa(b);
	  logp-= - mylgamma(sXa+sXb+getSumBeta(a))/getP();
	  double sumab=sumXab(a,b);
	  double sumba=sumXab(b,a);
	  int epop=getPsize(b);
	  int epopa=getPsize(a);
	  if(a==(int)b) {epop--; epopa--;}
	  if(epop+epopa>0) {logp+=mylgamma(getBeta(a,b)+sumab+sumba) - mylgamma(getBeta(a,b)) - (sumab+sumba)*log(epop+epopa);
	  }
	}
	return(logp);*/
}

double State::posteriorLengthSetProb(int a)
{
	double ret=0;
	for(unsigned int b=0;b<psize.size();b++) {ret+=posteriorLengthSetProb(a,b);}
	return(ret);
}

double State::posteriorLengthSetProb(int a,int b)
{
/*	vector<int> pa=getIndInPop(a);
	vector<int> pb=getIndInPop(b);
	cout<<"pop in a :";
	for(unsigned int c1=0;c1<pa.size();c1++) cout<<pa[c1]<<",";
	cout<<"   pop in b :";
	for(unsigned int c1=0;c1<pb.size();c1++) cout<<pb[c1]<<",";
	cout<<"    sumXab="<<calcSumXab(a,b,false)<<" sumLab="<<sumLab(a,b)<<endl;
*/
/*	if(fabs(sumXab(a,b,false)-calcSumXab(a,b,false))>0.00001) {cout<<sumXab(a,b,false)<<" - "<<calcSumXab(a,b,false)<<endl;exit(0);}
	if(fabs(sumLab(a,b)-calcSumLab(a,b))>0.00001) {cout<<sumLab(a,b)<<" - "<<calcSumLab(a,b)<<endl;exit(0);}
	if(fabs(sumlogLpi(a,b)-calcSumlogLpi(a,b))>0.00001) {cout<<sumlogLpi(a,b)<<" - "<<calcSumlogLpi(a,b)<<endl;exit(0);}
	if(fabs(sumlogLgamma(a,b)-calcSumlogLgamma(a,b))>0.00001) {cout<<sumlogLgamma(a,b)<<" - "<<calcSumlogLgamma(a,b)<<endl;exit(0);}
*/
	double ret=mylgamma(getLalpha(a,b) +  sumXab(a,b,false)) - mylgamma(getLalpha(a,b)) - sumlogLgamma(a,b);
	ret += getLalpha(a,b)*log(getLbeta(a,b)) + sumlogLpi(a,b) - (getLalpha(a,b) +  sumXab(a,b,false))*log(getLbeta(a,b) + sumLab(a,b));
	//	cout<<"a="<<a<<" b="<<b<<" postprob="<<ret<<endl;
//	cout<<"LG: "<<mylgamma(getLalpha(a,b) + calcSumXab(a,b,false)) <<" - "<<mylgamma(getLalpha(a,b))<<" - "<<sumlogLgamma(a,b)<<endl;
//	cout<<getLalpha(a,b)*log(getLbeta(a,b))<<" + "<<sumlogLpi(a,b) <<" - "<<(getLalpha(a,b) +  calcSumXab(a,b,false)-1.0)*log(getLbeta(a,b) + sumLab(a,b))<<endl;


	return(ret);
}

double State::posteriorNumChunksSetProb(int a)
{
	double effsamplesize=getpriorNumChunksParam(2);
	//effsamplesize=1.0;
	vector<int> pa=getIndInPop(a);
	double summu=0.0,sumlgamma=0.0, sumlpi=0.0,sumdat=0.0;
	for(unsigned int c1=0;c1<pa.size();c1++){
		summu+=mu(pa[c1]);
//		cout<<"mlg("<<effsamplesize*data->rowsum(pa[c1])<<")"<<endl;
		sumlgamma+=mylgamma(effsamplesize*data->rowsum(pa[c1]));
		sumlpi+=(effsamplesize*data->rowsum(pa[c1])-1) * log(mu(pa[c1]));
		sumdat+=effsamplesize*data->rowsum(pa[c1]);
	}
	//for(int i=0;i<pa.size();i++) cout<<mu(pa[i])<<", ";
	//cout<<"mean="<<summu/pa.size()<<endl;
	//double fac1=(pa.size()-1.0)/pa.size();
	//double fac2=1.0 - fac1;
//	cout<<"mlg2("<<sumdat + getpriorNumChunksParam(0)<<"),"<<getpriorNumChunksParam(0)<<endl;
	double ret=mylgamma(sumdat + getpriorNumChunksParam(0)) - mylgamma(getpriorNumChunksParam(0)) - sumlgamma + getpriorNumChunksParam(0)*log(getpriorNumChunksParam(1)) + sumlpi - (getpriorNumChunksParam(0) + sumdat)*log(getpriorNumChunksParam(1) + summu);
	//double modprior = mylgamma(sumdat*fac1 + getpriorNumChunksParam(0)) + mylgamma(sumdat*fac2 + getpriorNumChunksParam(0)) - mylgamma(sumdat + getpriorNumChunksParam(0));
	//modprior +=(getpriorNumChunksParam(0) + sumdat)*log(getpriorNumChunksParam(1) + summu) - (getpriorNumChunksParam(0) + fac1*sumdat)*log(getpriorNumChunksParam(1) + fac1*summu) - (getpriorNumChunksParam(0) + fac2*sumdat)*log(getpriorNumChunksParam(1) + fac2*summu);
	//ret-=modprior;
	return(ret);
}

double State::getApproxLogLikelihood()
{
  double logp=0.0;
  for(unsigned int a=0;a<psize.size();a++) {
	logp+=getApproxLogLikelihood(a);
  }
  return(logp);
}

double State::getApproxLogLikelihood(int a)
{
      std::vector<int> tmp=getIndInPop(a);
      std::vector<double> bestP(psize.size(),0);
      double sXa=sumXa(a);
      double logp=0.0;
      for(unsigned int b=0;b<psize.size();b++) {
	bestP[b] = (double)sumXab(a,b)/sXa;
	if(a==(int)b && psize[b]>1) logp+=sumXab(a,b) * log(bestP[b]/(psize[b]-1));
	else if(a!=(int)b) logp+=sumXab(a,b) * log(bestP[b]/psize[b]);
	cout<<"Approx LL for pop "<<a<< ","<<b<<" = "<<sumXab(a,b) * log(bestP[b]/psize[b])<<" (P_{ab}="<<bestP[b]<<")"<<endl;
      
      }
      cout<<"Approx LL for pop "<<a<<" = "<<logp<<endl;
      return(logp);
	
	// sum_b p_{ab} = 1
	// E(p_{ab}) = x_{ab} / \sum_b x_{ab}
	// likelihood = (p_{ab}/n_b)^{x_ab})
	// log-likelihood  = x_{ab} log(p_{ab}/n_b)
	
}



double State::admixtureLogLikelihood(std::vector<std::vector<double> > * Q,std::vector<std::vector<double> > * P, std::vector<double> *qColSums)
{
	double logp=0.0;
	for(int a=0;a<(int)Q->size();a++) {
		logp+=admixtureLogLikelihoodIndiv(a,Q,P,qColSums);
	}
	return(logp);
}

double State::admixtureLogLikelihoodIndiv(int a, std::vector<std::vector<double> > * Q,std::vector<std::vector<double> > * P, std::vector<double> *qColSums)
{
	double logp=0.0;
	for(int b=0;b<(int)Q->size();b++) {
		logp+=admixtureSetLogLikelihood(a,b,Q,P,qColSums);
	}
	return(logp);
}

double State::admixtureSetLogLikelihood(int a,int b,std::vector<std::vector<double> > * Q,std::vector<std::vector<double> > * P, std::vector<double> *qColSums)
{
	double prodsum=0.0;

	for(unsigned int c1=0;c1<psize.size();c1++) {
	  for(unsigned int c2=0;c2<psize.size();c2++) {
		if(qColSums->at(c2)-Q->at(a)[c2]!=0) prodsum+=Q->at(a)[c1] * P->at(c1)[c2] * Q->at(b)[c2]/(qColSums->at(c2)-Q->at(a)[c2]);
		// *** Is this normalised correctly and using the correct rows/columns?
	  }
	}
//	return((sumXab(a,b) + getBeta(a,b))*log(prodsum));
	if(prodsum<=0) return(0);
	return((X(a,b))*log(prodsum));
}

double State::admixtureApproxSetLogLikelihood(int set,std::vector<std::vector<double> > * Q,std::vector<std::vector<double> > * P)
{// Approximate ll for the admixture model treating copies from as fixed.
	vector<int> pset=getIndInPop(set);
	double prodsum=0.0;
	for(unsigned int c1=0;c1<pset.size();c1++) {
		prodsum+=admixtureApproxIndLogLikelihood(pset[c1],Q,P);
	}
	return(prodsum);
}

double State::admixtureApproxIndLogLikelihood(int ind,std::vector<std::vector<double> > * Q,std::vector<std::vector<double> > * P)
{// Approximate ll for the admixture model treating copies from as fixed.
	double prodsum=0.0;

	for(unsigned int c1=0;c1<psize.size();c1++) {
		prodsum+=admixtureApproxIndSetLogLikelihood(ind,c1,Q,P);
	}
	return(prodsum);
}

double State::admixtureApproxIndSetLogLikelihood(int ind,int pop,std::vector<std::vector<double> > * Q,std::vector<std::vector<double> > * P)
{
	double prodsum=0.0;

	for(unsigned int c1=0;c1<psize.size();c1++) {
//		if(getPsize(pop)-Q->at(ind)[pop] > 0) prodsum+=Q->at(ind)[c1] * P->at(c1)[pop]/(getPsize(pop)-Q->at(ind)[pop]);
		if(getPsize(pop)-Q->at(ind)[pop] > 0) prodsum+=Q->at(ind)[c1] * P->at(c1)[pop]/(getPsize(pop));
	}
	if(prodsum==0) {return(-10000000);}
	return(sumX(ind,pop)*log(prodsum));
//	return(sumY(pop,ind)*log(prodsum));
}


double State::mylgamma(double z)
{
/* LGAMMA function

   double_value = lgamma(<double_value > 0.>)

   returns the natural log of the gamma function

Uses Lanczos-type approximation to ln(gamma) for z > 0.
Reference:
 Lanczos, C. 'A precision approximation of the gamma
    function', J. SIAM Numer. Anal., B, 1, 86-96, 1964.

Original was in FORTRAN
Accuracy: About 14 significant digits except for small regions
          in the vicinity of 1 and 2.
Programmer: Alan Miller
          CSIRO Division of Mathematics & Statistics
Latest revision - 17 April 1988

Translated and modified into C by Peter Beerli 1997
Tested against Mathematica's Log[Gamma[x]]
*/

  double a[9] = { 0.9999999999995183, 676.5203681218835,
    -1259.139216722289, 771.3234287757674, -176.6150291498386,
    12.50734324009056, -0.1385710331296526, 9.934937113930748e-6,
    1.659470187408462e-7 };
  double lnsqrt2pi = 0.9189385332046727;
  double result;
  long j;
  double tmp;
  if (z <= 0.)
    {
	cerr<<"Lgamma error: z="<<z<<endl;
      throw(std::string("lgamma function failed with wrong input"));
    }
  result = 0.;
  tmp = z + 7.;
  for (j = 9; j >= 2; --j)
    {
      result += a[j - 1] / tmp;
      tmp -= 1.;
    }
  result += a[0];
  result = log (result) + lnsqrt2pi - (z + 6.5) + (z - 0.5) * log (z + 6.5);
  return result;
}  /* lgamma */

/*---------------------------------*/
double State::LDirichletProb(vector<double> prior,vector<double> post)
/*returns the log probability of a vector "post",
  given a Dirichlet process with prior "prior". */
{
  double sumprior = 0.0;
  double logsum;

  for (unsigned int i=0; i<prior.size(); i++){
//cout<<"prior["<<i<<"]="<<prior[i]<<endl;
    sumprior += prior[i];
  }
  logsum = mylgamma(sumprior);

  for (unsigned int i=0; i<prior.size(); i++){
//cout<<"post["<<i<<"]="<<post[i]<<endl;
	logsum += (prior[i]-1.0)*log(post[i]) - mylgamma(prior[i]);
  }
  return logsum;
}

double State::LGammaDistProb(double alpha,double beta, double y)
/*returns the log probability of a gamma-distributed random
variable "y", with parameters alpha and beta, where the mean
is alpha*beta, and the variance is alpha*beta*beta*/
{

  double logsum;

  logsum = -1*mylgamma(alpha) - alpha*log(beta) + (alpha-1)*log(y) - y/beta;

  return logsum;

}

void State::rpermute(vector<int>* in)
{
// Fisher-Knuth-Yates algorithm
	unsigned int swapto;
	int tmp;
	for(unsigned int max=in->size()-1;max>0;max--) {
	  swapto=RandomInteger(0,max);
	  if(swapto!=max) {
	    tmp=in->at(max);
	    in->at(max)=in->at(swapto);
	    in->at(swapto)=tmp;
	  }
	}
}

void State::rpermute(vector<int>* in1,vector<int>* in2)
{
	if(in1->size()!=in2->size()) throw(std::string("Error in rpermute: Cannot idenitcally permute vectors of different length"));
	else if(in1->size()==0) return;
// Fisher-Knuth-Yates algorithm
	unsigned int swapto;
	int tmp1,tmp2;
	for(unsigned int max=in1->size()-1;max>0;max--) {
	  swapto=RandomInteger(0,max);
	  if(swapto!=max) {
	    tmp1=in1->at(max);
	    in1->at(max)=in1->at(swapto);
	    in1->at(swapto)=tmp1;
	    tmp2=in2->at(max);
	    in2->at(max)=in2->at(swapto);
	    in2->at(swapto)=tmp2;
	  }
	}
}

double State::sampleBetaP(int param, int modelpart)
{
	if(hyperprior.size()<3) throw("Hyperprior too small for Beta model!");
	if(modelpart==0){
		if(param==0) {// its F for the F model i.e. in (0,1)
			if(hyperprior[0]<0) return(-hyperprior[0]);
			double r=-1;
			while(r<0 || r>1) r=RGamma(hyperprior[0],1.0/hyperprior[2]);
			return(r);
		}else if(param==1) {
			if(hyperprior[1]<0) return(-hyperprior[1]);
			return(RGamma(hyperprior[1],1.0/hyperprior[3]));
		}
	}else if(modelpart==1){
		vector<double> tvec=getLengthHyperPrior();
		return(RGamma(tvec[param],1.0/tvec[param+NUMHYPERPARAMLENGTH]));
	}else if(modelpart==2){
		vector<double> tvec=getTotalChunksHyperPrior();
		return(RGamma(tvec[param],1.0/tvec[param+NUMHYPERPARAMTOTLENGTH]));
	}else throw(string("Unimplemented model requested in samplebetaP"));
	return(-1);// shouldn't get here
}

State::~State()
{
}

void State::addCoincidence(vector< vector<double> > *vecin)
{
	for(int c1=0;c1<getP();c1++) {
	  std::vector<int> ilist=getIndInPop(c1);
	  for(unsigned int c2=0;c2<ilist.size();c2++) {
	  for(unsigned int c3=c2;c3<ilist.size();c3++) {
		vecin->at(ilist[c2])[ilist[c3]]++;
		if(c3!=c2) vecin->at(ilist[c3])[ilist[c2]]++;
	  }}
	}
}

double State::getDistanceSq(vector< vector<double> > *meanX)
{
	double d=0;
	for(unsigned int c1=0;c1<meanX->size();c1++) {
	  for(unsigned int c2=c1+1;c2<meanX->at(c1).size();c2++) {
		d += (double)(getPop(c2)==getPop(c1)) - meanX->at(c1)[c2];
	  }
	}
	return(d);
}

double State::popDist(vector< vector<double> > *meanX,int pop,int indiv,double penalty)
{
	double td=0.0;
	std::vector<int> ilist=getIndInPop(pop);
	double mx;
	for(unsigned int c1=0;c1<ilist.size();c1++) {if(indiv!=ilist[c1]) {
	  mx=pow(meanX->at(indiv)[ilist[c1]],penalty);
	  td+= mx*mx - (1.0 - mx)*(1.0 - mx);
	}}
	return(td);
}

double State::minDistanceMove(vector< vector<double> > *meanX,double penalty)
{
	double distchange=0.0,td,tdtemp;
	int bestpopto=-1,bestind=-1;
	int opop;
	for(unsigned int c1=0;c1<meanX->size();c1++) {
	  opop=getPop(c1);
	  td=popDist(meanX,opop,c1,penalty);
	  for(int c2=0;c2<getP();c2++) {// test moving to each population
		if(c2!=opop) {
		  tdtemp=-popDist(meanX,c2,c1,penalty);
		  if(td+tdtemp<distchange){bestind=c1;bestpopto=c2;distchange=td+tdtemp;}
		}
	  }
	  if(td<distchange){bestind=c1;bestpopto=-1;distchange=td;}// test moving to new population
	}
	if(bestind>=0) {// move to best population
	  opop=getPop(bestind);
	  if(bestpopto<0){
		moveInd(bestind,addEmptyPop());
	  }else{
		moveInd(bestind,bestpopto);
	  }
	  if(getPlength(opop)==0) removePop(opop);
	}
	return(distchange);
//	for(int c1=0;c1<getP();c1++) {
}

std::vector<double> State::readOutputVectorD(std::string str){
		std::vector<double> res;
	  	size_t f1=str.find('>'),f2=str.find('<',f1);
		str=str.substr(f1+1,f2-f1-1);
	f2=str.find(',');f1=0;
	res.push_back(atof(str.substr(f1,f2-f1).c_str()));
	while(f2!=std::string::npos){
	f1=f2;f2=str.find(',',f1+1);
	res.push_back(atof(str.substr(f1+1,f2-f1-1).c_str()));
	}
	return(res);
}

void State::reorderPop(int from, int to)
{// insert "to" into "from" in all the vectors where these are accounted for
    if(from>=getP()||from<0) throw(string("reorderPop:from invalid!"));
    if(to>=getP()||to<0) throw(string("reorderPop:to invalid!"));
    int efffrom=from;
    if(from>to)efffrom++;// account for insertion
	//cout<<"rop: From="<<from<<" to="<<to<<" efffrom="<<efffrom<<" (psize="<<psize.size()<<")"<<endl;
    psize.insert(psize.begin()+to,psize[from]);
    psize.erase(psize.begin()+efffrom);
    if(beta.size()==psize.size()){
    	beta.insert(beta.begin()+to,beta[from]);
    	beta.erase(beta.begin()+efffrom);
    }
    if(betaF.size()==psize.size()){
    	betaF.insert(betaF.begin()+to,betaF[from]);
    	betaF.erase(betaF.begin()+efffrom);
    }
    if(delta.size()==psize.size()){
    	delta.insert(delta.begin()+to,delta[from]);
    	delta.erase(delta.begin()+efffrom);
    }
    popX.insert(popX.begin()+to,popX[from]);
    popX.erase(popX.begin()+efffrom);
    if(uselengths) {
    	popL.insert(popL.begin()+to,popL[from]);
    	popL.erase(popL.begin()+efffrom);
	poplogLpi.insert(poplogLpi.begin()+to,poplogLpi[from]);
    	poplogLpi.erase(poplogLpi.begin()+efffrom);
	poplogLgamma.insert(poplogLgamma.begin()+to,poplogLgamma[from]);
    	poplogLgamma.erase(poplogLgamma.begin()+efffrom);
    }
    indinp.insert(indinp.begin()+to,indinp[from]);
    indinp.erase(indinp.begin()+efffrom);
}

void State::reorderIndiv(int pop,int from, int to)
{// insert individual "from" at location "to" in population pop (before element at to)
    if(pop>=getP()||pop<0) throw(string("reorderIndiv:pop invalid!"));
    if(from>=(int)popX[pop].size()||from<0) {cerr<<"From="<<from<<"; valid is(0,"<<popX[pop].size()<<")"<<endl;throw(string("reorderIndiv:from invalid!"));}
    if(to>=(int)popX[pop].size()||to<0) {cerr<<"To="<<to<<"; valid is(0,"<<popX[pop].size()<<")"<<endl;throw(string("reorderIndiv:to invalid!"));}
    int efffrom=from;
    if(from>to)efffrom++;// account for insertion
//	cout<<"From="<<from<<" to="<<to<<" efffrom="<<efffrom<<" (of "<<getIndInPop(pop).size()<<")"<<endl;
//	cout<<"From="<<data->getnames(getIndInPop(pop)[from])<<" to="<<data->getnames(getIndInPop(pop)[to])<<endl;
    popX[pop].insert(popX[pop].begin()+to,popX[pop][from]);
    popX[pop].erase(popX[pop].begin()+efffrom);
    if(uselengths) {
    popL[pop].insert(popL[pop].begin()+to,popL[pop][from]);
    popL[pop].erase(popL[pop].begin()+efffrom);
    poplogLpi[pop].insert(poplogLpi[pop].begin()+to,poplogLpi[pop][from]);
    poplogLpi[pop].erase(poplogLpi[pop].begin()+efffrom);
    poplogLgamma[pop].insert(poplogLgamma[pop].begin()+to,poplogLgamma[pop][from]);
    poplogLgamma[pop].erase(poplogLgamma[pop].begin()+efffrom);
    }
    indinp[pop].insert(indinp[pop].begin()+to,indinp[pop][from]);
    indinp[pop].erase(indinp[pop].begin()+efffrom);
}

/*int State::movePcaAutoInd(int testi,int i, int j,PcaData *pca)
{
    if(pca==NULL) return(0);
    double targd=pca->getDist(i,j);
    double di=pca->getDist(testi,i);
    double dj=pca->getDist(testi,j);
    cout<<"d="<<d1<<","<<targd<<","<<d2<<":"<<getCosAngle(d1,targd,d2)<<endl;
    cout<<"d="<<d2<<","<<targd<<","<<d1<<":"<<getCosAngle(d2,targd,d1)<<endl;
    if(d1>targd && d1<d2) {nc++; moveInd(opop[c1],ind[i]);
    }else if(d2<targd && d2<d1) {nc++; moveInd(opop[c1],ind[j]);}
//    if(isObtuseAngle(d1,targd,d2)) {nc++; moveInd(opop[c1],ind[i]);
//    }else if(isObtuseAngle(d2,targd,d1)) {nc++; moveInd(opop[c1],ind[j]);}
    if(isObtuseAngle(di,targd,dj)) { //moveInd(testi,ind[i]);
      return(1);
    }else if(isObtuseAngle(dj,targd,di)) { //moveInd(testi,ind[j]);
      return(2);
    }
    return(0);
} */

bool State::isObtuseAngle(double d1, double d2, double d3){
  //return(false);
  return((d1*d1 +d2*d2 - d3*d3)/2/d1/d2 < 0.070 && d1<d2);
  //cout<<d3*d3 -(d1*d1 +d2*d2)<<":::"<<(d3*d3 > d1*d1 +d2*d2)<<endl;
//  return(d3*d3 > d1*d1 +d2*d2);
//  return(d1<d3 && d1<d2);
}
  
  /*
void State::movePcaAuto(int opindex,int i, int j,double targd,PcaData *pca)
{
  int nc=0;
  vector<int> opop=getIndInPop(opindex);  
  for(unsigned int c1=0;c1<opop.size();c1++){
    if(movePcaAutoInd(opop[c1],i,j,targd,pca)>0) nc ++;
  }
  cout<<"movePcaAuto nc="<<nc<<" of "<<opop.size()<<endl;
}

bool State::testPcaAuto(int a, int b,int i, int j,double targd,PcaData *pca)
{
  int nc=0, nf=0;
  vector<int> opopa=getIndInPop(a);  
  vector<int> opopb=getIndInPop(b);  
  bool valid=true;
  for(unsigned int c1=0;c1<opopa.size();c1++){
    double d1=pca->getDist(opopa[c1],i);
    double d2=pca->getDist(opopa[c1],j);
    if(isObtuseAngle(d1,targd,d2)) {nc++; moveInd(opopa[c1],ind[i]);
    }else if(isObtuseAngle(d2,targd,d1)) {nc++; moveInd(opopa[c1],ind[j]);valid=false;nf++;}
  }
  for(unsigned int c1=0;c1<opopb.size();c1++){
    double d1=pca->getDist(opopb[c1],i);
    double d2=pca->getDist(opopb[c1],j);
    if(isObtuseAngle(d1,targd,d2)) {nc++; moveInd(opopb[c1],ind[i]);valid=false;nf++;
    }else if(isObtuseAngle(d2,targd,d1)) {nc++; moveInd(opopb[c1],ind[j]);}    
  }
  cout<<"testPcaAuto nc="<<nc<<" of "<<opopa.size()+opopb.size()<<endl;
  if(!valid) cout<<"Invalid merge: nf="<<nf<<" of "<<opopa.size()+opopb.size()<<endl;
  return(valid);
}
*/

bool State::checkPopForIgnoredSuper(int a) {
  vector<int> pa=getIndInPop(a);
  for(unsigned int c1=0;c1<pa.size();c1++) {
    if(data->getIgnore(pa[c1]) && pa.size()>1) return(true);
  }
  return(false);
}
  
} // end namespace weakarg

