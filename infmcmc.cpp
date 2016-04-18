#include "infmcmc.h"
#include "math.h"

using namespace std;
namespace fines
{

InfMCMC::InfMCMC(Data *d,State *s,Data *d2,int datainference,double pcaprob, bool v)
{
	data=d;
	dlength=d2;
	this->datainference=datainference;
	setdataused();
	state=s;
	verbose=v;
	logstateprob=state->posteriorProb();
	resetCounters();
	P_SAMS=MCMCPROB/MCMCNUMVARS;
	P_MERGESPLIT=P_SAMS+MCMCPROB/MCMCNUMVARS;
	P_IND=P_MERGESPLIT+MCMCPROB/MCMCNUMVARS;
	//P_HYPER =1-sum(these)
}

  void InfMCMC::fixK()
  {
	P_SAMS=0;
	P_MERGESPLIT=P_SAMS+MCMCPROB/(MCMCNUMVARS-1.0);
	P_IND=P_MERGESPLIT+MCMCPROB/(MCMCNUMVARS-1.0);
    
  }


bool InfMCMC::moveMergeAndSplit(bool greedy)
{
	if(state->getP()==1) return(0);// Cannot merge only one population
	double logpofsplit1=0.0,logpofsplit2=0.0;
	double logposteriorold=logstateprob,logposteriornew;
	int i=RandomInteger(0,data->getDim()-1),j=i;
	// force a merge of two different populations
	while(i==j || state->getPop(i)==state->getPop(j)) j=RandomInteger(0,data->getDim()-1);
	if(state->getPlength(state->getPop(i))==1 && state->getPlength(state->getPop(j))==1){
		if(verbose) {cout<<"Rejected: Single size populations"<<endl;
			state->setprint(&cout);}
		return(0);// will always give identical populations
	}
	if(verbose){cout<<"MergeAndSplit: INITIAL STATE:"<<endl;
		state->setprint(&cout);}
 	logpofsplit1= state->probOfSplitSAMS(i,j);
	State *oldstate=new State(state);
	state->merge(state->getPop(i),state->getPop(j));
	logpofsplit2 = state->splitSAMS(i,j,false);
	logposteriornew=state->posteriorProb();

	  if(verbose){cout<<"PROPOSED STATE:"<<endl;
		state->setprint(&cout);
	  	cout<<"PROPOSED MERGE&SPLIT LOG PROB="<<logpofsplit1<<"-"<<logpofsplit2<<endl;
	  }

	if((!greedy && log(rnd())<logposteriornew-logposteriorold + logpofsplit2-logpofsplit1 && fabs(logposteriorold-logposteriornew)>1e-07) || (greedy && logposteriornew>logposteriorold)) {
// insist we only accept if the states are different
		if(verbose) {cout<<"Accepted:"<<"LL="<<logposteriornew<<" versus oldLL="<<logposteriorold<<endl;
			state->setprint(&cout);}
		logstateprob=logposteriornew;
		delete(oldstate);
		return(true);
	}else{
		if(verbose) {cout<<"Rejected:"<<"LL="<<logposteriornew<<" versus oldLL="<<logposteriorold<<endl;
			state->setprint(&cout);}
		delete(state);
		state=oldstate;
		return(false);
	}
}

bool InfMCMC::moveSAMS(bool greedy)
{
	if(data->getDim()==1) return(0);// cannot merge-split on one individual
	double logpofsplit=0.0;
	double logposteriorold=logstateprob,logposteriornew;
// Choose two random individuals
	int i=RandomInteger(0,data->getDim()-1),j=i;
	while(i==j) j=RandomInteger(0,data->getDim()-1);
// Split or merge as appropiate
	if(verbose){cout<<"BEFORE MOVE"<<endl;state->setprint(&cout);}
 	if(state->getPop(i)==state->getPop(j)) {// try a split
	  logpofsplit = state->splitSAMS(i,j,false);
	  logposteriornew=state->posteriorProb();
	  if(verbose){cout<<"PROPOSED STATE:"<<endl;
		state->setprint(&cout);
	  	cout<<"PROPOSED SPLIT WITH LOG PROB="<<logpofsplit<<endl;
	  }
	  if((!greedy && log(rnd())<logposteriornew-logposteriorold - logpofsplit) || (greedy && logposteriornew>logposteriorold)) {// accept split
		if(verbose) {cout<<"Accepted:"<<"LL="<<logposteriornew<<" versus oldLL="<<logposteriorold <<"PofSplit(-)"<<logpofsplit <<endl;
			state->setprint(&cout);}
		logstateprob=logposteriornew;
		return(true);
	  }else {
		state->merge(state->getPop(i),state->getPop(j));
		if(verbose) {cout<<"Rejected:"<<"LL="<<logposteriornew<<" versus oldLL="<<logposteriorold <<"PofSplit(-)"<<logpofsplit <<endl;
			state->setprint(&cout);}
		return(false);
	  }
	}else { // try a merge
	  logpofsplit= state->probOfSplitSAMS(i,j);
	  State *oldstate=new State(state);
	  state->merge(state->getPop(i),state->getPop(j));
	  logposteriornew=state->posteriorProb();
	  if(verbose){cout<<"PROPOSED STATE:"<<endl;
		state->setprint(&cout);
	  	cout<<"PROPOSED MERGE WITH split LOG PROB="<<logpofsplit<<endl;
	  }
	  if(log(rnd())<logposteriornew-logposteriorold + logpofsplit || (greedy && logposteriornew>logposteriorold)) {// accept merge
		logstateprob=logposteriornew;
		if(verbose) {cout<<"Accepted:"<<"LL="<<logposteriornew<<" versus oldLL="<<logposteriorold <<"PofSplit(+)"<<logpofsplit <<endl;
			state->setprint(&cout);}
		delete(oldstate);
		return(true);
	  }else {
		delete(state);
		state=oldstate;
		if(verbose) {cout<<"Rejected:"<<"LL="<<logposteriornew<<" versus oldLL="<<logposteriorold <<"PofSplit(+)"<<logpofsplit <<endl;
			state->setprint(&cout);}
		return(false);
	  }
	}
	// should never get here:
	return(false);
}

bool InfMCMC::moveIndiv(bool greedy)
{
	if(state->getP()==1) return(false);// can't move if only one population
// Choose a random individual
	int i=RandomInteger(0,data->getDim()-1);
	if(state->getPlength(state->getPop(i))==1) return(false);//cant undo moving single individual; thats a split/merge move
	  
// choose a random other population
	int oldpop = state->getPop(i),newpop=oldpop;
	while(newpop==oldpop) newpop=RandomInteger(0,state->getP()-1);

	vector <int> tmppop=state->getIndInPop(state->getPop(i));
	int otherindina=-1;
	for(int j=0;j<(int)tmppop.size();j++) if(tmppop[j]!=i) {otherindina=tmppop[j];break;}
	if(otherindina<0) throw(string("Error in moveIndiv: Can't find any other individuals!"));
	State *mergestateold=new State(state);
	mergestateold->merge(mergestateold->getPop(i),newpop);
	
	double logposteriorold=state->posteriorProb(mergestateold);
	//double logposteriorold=state->posteriorProb(NULL);
	if(verbose) {cout<<"MOVEINDIV: BEFORE ";state->setprint(&cout);}
	state->moveInd(i,newpop);
	if(verbose) {cout<<"MOVEINDIV: PROPOSED ";state->setprint(&cout);}
	State *mergestatenew=new State(state);
	mergestatenew->merge(mergestatenew->getPop(otherindina),newpop);
	double logposteriornew=state->posteriorProb(mergestatenew);
	//double logposteriornew=state->posteriorProb(NULL);
	delete(mergestateold);
	delete(mergestatenew);
// Acceptance/rejection step
	if(verbose) {cout<<"MOVEINDIV: log(p) = "<<logposteriornew<<" - "<<logposteriorold<<endl;}
	if( (!greedy && log(rnd())< logposteriornew-logposteriorold) || (greedy && logposteriornew>logposteriorold)) {
		//accept
	if(verbose) cout<<"MOVEINDIV ACCEPTED"<<endl;
		logstateprob=logposteriornew;
		return(true);
	}else {
	if(verbose) cout<<"MOVEINDIV REJECTED"<<endl;
		state->moveInd(i,oldpop);
		return(false);
	}
}

bool InfMCMC::moveF(int a,bool greedy)
{
    double d0=state->getBetaF(a);
    double d1=exp(log(d0) + (rnd()*2.0-1.0)*1);
    while(d1<=0 || d1>=1) {
	if(d1<0) d1=-d1;
	if(d1>=1) d1=2-d1;
	if(d1==0) d1=rnd()*0.001;
	if(d1==1) d1=1.0*rnd() *0.001;
    }
    double p0=state->posteriorProb();
    state->setBetaF(a,d1);
    double p1=state->posteriorProb();
    double postdiff=p1-p0 + state->LGammaDistProb(state->getHyperPrior()[0],state->getHyperPrior()[2],d1) - state->LGammaDistProb(state->getHyperPrior()[0],state->getHyperPrior()[2],d0);
    if((!greedy && log(rnd())< postdiff) || (greedy &&postdiff>0)) {
	if(verbose) cout<<"MOVEF ACCEPTED: from "<<d0<<" to "<<d1<<endl;
	return(true);
    }else {
	state->setBetaF(a,d0);
	if(verbose) cout<<"MOVEF REJECTED"<<endl;
	return(false);
    }
}

bool InfMCMC::moveDelta(int a,bool greedy)
{
//   state->setDelta(a,state->sampleBetaP(true));
//   return(true);
    double d0=state->getDelta(a);
    double d1=exp(log(d0) + (rnd()*2.0-1.0)*1);
    if(d1<0) d1=-d1;
    double p0=state->posteriorProb();
    state->setDelta(a,d1);
    double p1=state->posteriorProb();
//cout<<"d1="<<d1<<", p1="<<p1<<" + "<<state->LGammaDistProb(state->getHyperPrior()[1],state->getHyperPrior()[3],d1)<<" vs d0="<<d0<<", p0="<<p0<<" + "<<state->LGammaDistProb(state->getHyperPrior()[1],state->getHyperPrior()[3],d0)<<", diff="<<p1-p0<<endl;
    double postdiff=p1-p0 + state->LGammaDistProb(state->getHyperPrior()[1],state->getHyperPrior()[3],d1) - state->LGammaDistProb(state->getHyperPrior()[1],state->getHyperPrior()[3],d0);
    if((!greedy && log(rnd())< postdiff)||(greedy && postdiff>0)) {
	if(verbose) cout<<"MOVEDELTA ACCEPTED: from "<<d0<<" to "<<d1<<endl;
	return(true);
    }else {
	state->setDelta(a,d0);
	if(verbose) cout<<"MOVEDELTA REJECTED"<<endl;
	return(false);
    }
}

bool InfMCMC::moveHyperParamLength(int par,bool greedy)
{
    double d0=state->getpriorLengthsParam(par);
    double d1=exp(log(d0) + (rnd()*2.0-1.0)*1);
    if(d1<0) d1=-d1;
    double p0=state->posteriorProb();
    state->setpriorLengthsParam(par,d1);
    double p1=state->posteriorProb();
    vector<double>hypervec=state->getLengthHyperPrior();
//cout<<"d1="<<d1<<", p1="<<p1<<" + "<<state->LGammaDistProb(state->getHyperPrior()[1],state->getHyperPrior()[3],d1)<<" vs d0="<<d0<<", p0="<<p0<<" + "<<state->LGammaDistProb(state->getHyperPrior()[1],state->getHyperPrior()[3],d0)<<", diff="<<p1-p0<<endl;
    double postdiff=p1-p0 + state->LGammaDistProb(hypervec[par],hypervec[par+NUMHYPERPARAMLENGTH],d1) - state->LGammaDistProb(hypervec[par],hypervec[par+NUMHYPERPARAMLENGTH],d0);
    if((!greedy && log(rnd())< postdiff)||(greedy && postdiff>0)) {
	if(verbose) cout<<"MOVEHYPERPARAMLENGTH "<<par<<" ACCEPTED: from "<<d0<<" to "<<d1<<endl;
	return(true);
    }else {
	state->setpriorLengthsParam(par,d0);
	if(verbose) cout<<"MOVEHYPERPARAMLENGTH "<<par<<" REJECTED"<<endl;
	return(false);
    }
}

bool InfMCMC::moveHyperParamTotLength(int par,bool greedy)
{
    double d0=state->getpriorNumChunksParam(par);
    double d1=exp(log(d0) + (rnd()*2.0-1.0)*1);
    if(d1<0) d1=-d1;
    double p0=state->posteriorProb();
    state->setpriorNumChunksParam(par,d1);
    double p1=state->posteriorProb();
    vector<double>hypervec=state->getTotalChunksHyperPrior();
//cout<<"d1="<<d1<<", p1="<<p1<<" + "<<state->LGammaDistProb(state->getHyperPrior()[1],state->getHyperPrior()[3],d1)<<" vs d0="<<d0<<", p0="<<p0<<" + "<<state->LGammaDistProb(state->getHyperPrior()[1],state->getHyperPrior()[3],d0)<<", diff="<<p1-p0<<endl;
    double postdiff=p1-p0 + state->LGammaDistProb(hypervec[par],hypervec[par+NUMHYPERPARAMTOTLENGTH],d1) - state->LGammaDistProb(hypervec[par],hypervec[par+NUMHYPERPARAMTOTLENGTH],d0);
    if((!greedy && log(rnd())< postdiff)||(greedy && postdiff>0)) {
	if(verbose) cout<<"MOVEHYPERPARAMTOTLENGTH "<<par<<" ACCEPTED: from "<<d0<<" to "<<d1<<endl;
	return(true);
    }else {
	state->setpriorNumChunksParam(par,d0);
	if(verbose) cout<<"MOVEHYPERPARAMTOTLENGTH "<<par<<" REJECTED"<<endl;
	return(false);
    }
}

double InfMCMC::moveHyper(bool greedy)
{
	double pacc=0;
	double denom=0;
	vector<double> hyper1=state->getHyperPrior(),hyper2,hyper3;
	if(uselengths)hyper2=state->getLengthHyperPrior();
	if(usesums)hyper3=state->getTotalChunksHyperPrior();
	for(int c1=0;c1< 2*state->getHyperParamLength();c1++) if(hyper1[c1]>0) denom++;
	if(uselengths) {for(unsigned int c1=0;c1< hyper2.size();c1++) if(hyper2[c1]>0) denom++;}
	if(usesums) {for(unsigned int c1=0;c1< hyper3.size();c1++) if(hyper3[c1]>0) denom++;}

	if(state->getBetaModel()==BETAMOD_F2) {
		if(usecounts){
			if(hyper1[0]>0) pacc+=moveDelta(0)/denom; 
			if(hyper1[1]>0) pacc+=moveF(0)/denom;
		}
	}
	else if(state->getBetaModel()==BETAMOD_F) {
		for(int i=0;i<state->getHyperParamLength();i++) {
			if(usecounts){
			if(hyper1[0]>0) pacc+=moveDelta(i)/denom; 
			if(hyper1[1]>0) pacc+=moveF(i)/denom;
			}
		}
	}
	if(uselengths){
		for(unsigned int i=0;i<NUMHYPERPARAMLENGTH;i++)  if(hyper2[i]>0) pacc+=moveHyperParamLength(i)/denom; 
	}
	if(usesums){
		for(unsigned int i=0;i<NUMHYPERPARAMTOTLENGTH;i++) if(hyper3[i]>0) pacc+=moveHyperParamTotLength(i)/denom;
	}
	return(pacc);
}

void InfMCMC::runone(long iter,long thin,std::ostream * fileout){
  	double x=rnd();
	if(x<P_SAMS) {
		numSAMS++;
		accSAMS+=(int)moveSAMS();
	}else if(x<P_MERGESPLIT) {
		numMergeSplit++;
		accMergeSplit+=(int)moveMergeAndSplit();
	}else if(x<P_IND) {
		numIndiv++;
		accIndiv+=moveIndiv();
	}else {
		numHyper++;
		accHyper+=(int)moveHyper();
	}
	if(iter % thin==0 && fileout!=NULL) exportXmlIter(fileout,iter);
}

void InfMCMC::metropolis(long prevints,long numints, long thin,std::ostream * fileout,long totallength)
{
   if(totallength<0)totallength=numints;
	for(long c1=prevints;c1<prevints+numints;c1++) {
		if (totallength>50 && (c1)%((totallength)/50)==0)
		{
			if (100l*c1/(totallength)<=1){
				cout<<"#  "<<100l*(double)c1/(numints)<<"%"<<flush;
			}else if (100l*c1/(totallength)<10)
				cout<<"\b\b\b\b#  "<<(int)(100l*(double)c1/(totallength))<<"%"<<flush;
			else
				cout<<"\b\b\b\b# "<<(int)(100l*(double)c1/(totallength))<<"%"<<flush;
		}
		if (c1+1==totallength)
			cout<<"\b\b\b\b# 100%"<<endl<<flush;
		runone(c1,thin,fileout);
	}
}


void InfMCMC::hillClimb(long prevints,long numints, long thin,std::ostream * fileout)
{
	for(long c1=0;c1<numints;c1++) {
		if (numints>50 && (c1)%((numints)/50)==0)
		{
			if (100l*c1/(numints)<=1){
				cout<<"#  "<<100l*(double)c1/(numints)<<"%"<<flush;
			}else if (100l*c1/(numints)<10)
				cout<<"\b\b\b\b#  "<<(int)(100l*(double)c1/(numints))<<"%"<<flush;
			else
				cout<<"\b\b\b\b# "<<(int)(100l*(double)c1/(numints))<<"%"<<flush;
		}
		if (c1+1==numints)
			cout<<"\b\b\b\b# 100%"<<endl<<flush;
		double x=rnd();
		if(x<P_SAMS) {
			numSAMS++;
			accSAMS+=(int)moveSAMS(true);
		}else if(x<P_MERGESPLIT) {
			numMergeSplit++;
			accMergeSplit+=(int)moveMergeAndSplit(true);
		}else if(x<P_IND) {
			numIndiv++;
			accIndiv+=moveIndiv(true);
		}else {
			numHyper++;
			accHyper+=(int)moveHyper(true);
		}
		if(c1 % thin==0 && fileout!=NULL) exportXmlIter(fileout,c1+prevints);
	}
}

void InfMCMC::exportXmlHead(std::ostream * fileout,double c,long x,long y,long z)
{
	*fileout<<"<?xml version = '1.0' encoding = 'UTF-8'?>"<<endl<<"<outputFile>"<<endl;
  *fileout<< "<header>"<<endl;
  *fileout<< "<runtype>MCMC</runtype>"<<endl;
  *fileout<<setprecision(12);
  if(c>=0) *fileout<< "<inflation>"<<c<< "</inflation>"<<endl;
  if(x>=0) *fileout<< "<burnin>"<<x<< "</burnin>"<<endl;
  if(y>=0) *fileout<< "<mcmclength>"<<y<< "</mcmclength>"<<endl;
  if(z>=0) *fileout<< "<skip>"<<z<< "</skip>"<<endl;
  *fileout<<"<datafilename>"<<data->getFileName()<<"</datafilename>"<<endl;
  *fileout<< "<copymodel>"<<state->getModelType()<< "</copymodel>"<<endl;    
  *fileout<<"</header>"<<endl;	
}

void InfMCMC::exportXmlTail(std::ostream * fileout)
{
	*fileout<<"</outputFile>"<<endl;
}

void InfMCMC::exportXmlIter(std::ostream * fileout,int iter)
{
//	state->printBeta(&cout);
	*fileout<<"<Iteration>"<<endl;
	state->setprint(fileout);
	*fileout<<setprecision (32)<<"<Posterior>"<<logstateprob<<"</Posterior>"<<setprecision (6)<<endl;
	*fileout<<"<K>"<<state->getP()<<"</K>"<<endl;
	*fileout<<"<alpha>"<<state->getAlpha()<<"</alpha>"<<endl;
	*fileout<<"<beta>"<<state->getSumBeta()<<"</beta>"<<endl;
	*fileout<<"<delta>";
	state->printDelta(fileout);
	*fileout<<"</delta>"<<endl;
	*fileout<<"<F>";
	state->printBetaF(fileout);
	*fileout<<"</F>"<<endl;
	if(uselengths){
	*fileout<<"<lengthalpha0>"<<state->getpriorLengthsParam(0)<<"</lengthalpha0>"<<endl;
	*fileout<<"<lengthbeta0>"<<state->getpriorLengthsParam(1)<<"</lengthbeta0>"<<endl;	
	*fileout<<"<lengthdeltaalpha>"<<state->getpriorLengthsParam(2)<<"</lengthdeltaalpha>"<<endl;	
	*fileout<<"<lengthdeltabeta>"<<state->getpriorLengthsParam(3)<<"</lengthdeltabeta>"<<endl;		
	}
	if(usesums){
	*fileout<<"<meanmualpha>"<<state->getpriorNumChunksParam(0)<<"</meanmualpha>"<<endl;
	*fileout<<"<meanmubeta>"<<state->getpriorNumChunksParam(1)<<"</meanmubeta>"<<endl;	
	*fileout<<"<meanmugamma>"<<state->getpriorNumChunksParam(2)<<"</meanmugamma>"<<endl;	
	}
	double ns=0,ni=0,nms=0,nh=0;
	if (numSAMS>0) ns=((double)accSAMS)/numSAMS;
	if (numIndiv>0) ni=((double)accIndiv)/numIndiv;
	if (numMergeSplit>0) nms=((double)accMergeSplit)/numMergeSplit;
	if (numHyper>0) nh=((double)accHyper)/numHyper;
	*fileout<<"<AccSAMS>"<<ns<<"</AccSAMS>"<<endl;
	*fileout<<"<AccMS>"<<nms<<"</AccMS>"<<endl;
	*fileout<<"<AccIndiv>"<<ni<<"</AccIndiv>"<<endl;
	*fileout<<"<AccHyper>"<<nh<<"</AccHyper>"<<endl;	*fileout<<"<Number>"<<iter<<"</Number>"<<endl;
	*fileout<<"</Iteration>"<<endl;
}

InfMCMC::~InfMCMC()
{
	if(state!=NULL) delete(state);
}


} // end namespace fines
