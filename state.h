#ifndef STATE_H
#define STATE_H
#include "data.h"
#include "fines.h"
#include "fsxml.h"
//#include "pcadata.h"
#include "safegetline.h"

#include <vector>
#include <string>
#include <iostream>
#include <math.h>
#include <float.h>

namespace fines
{

#define NUMHYPERPARAMLENGTH 4
#define NUMHYPERPARAMTOTLENGTH 3

/**
    @brief population structure class
*/
class State
{
public:

// ***** Creation and deletion of states

    State(Data *d,int npop,std::vector<double> bvec,double a,int betamodel,double corfactor,Data *d2=NULL,int datainference=INFDATA_COUNTS,int modeltype=MODELTYPE_FINESTRUCTURE);///<Constructs a population structure for nind individuals and npop populations.  Initialises beta with priorparam and alpha with a
    State(Data *d,FsXml *infile,std::vector<double> bvec,double a,int betamodel,double corfactor,bool readp=false,Data *d2=NULL,int datainference=INFDATA_COUNTS,int modeltype=MODELTYPE_FINESTRUCTURE);///<Gets the final iteration from a valid output file (using specified parameters)
//    State(Data *d,FsXml *infile,int betamodel,double corfactor);///<Gets from a valid output file and uses the parameters found
    State(Data *d,std::vector< std::vector<double> > *vecin,bool mergerule, double mergeval,std::vector<double> bvec,double a,int betamodel,double corfactor,Data *d2=NULL,int datainference=INFDATA_COUNTS,int modeltype=MODELTYPE_FINESTRUCTURE);///<Create a state from an average coincidence matrix
    State(Data *d,std::string stringin,std::vector<double> bvec,double a,int betamodel,bool isfile,double corfactor,Data *d2=NULL,int datainference=INFDATA_COUNTS,int modeltype=MODELTYPE_FINESTRUCTURE);///<Construct a population structure from our output format string
    State(State *in);///<Copies the state from another
    inline double getCorfactor(){
      return(corfactor);
    };///<returns the correlation factor
    void createParams(double a);///<Creates parameters and population copy matrix
    std::string getFromFile(FsXml *infile,bool readparams=false);///<Gets a state from an output file
    void readState(std::string statein);///< Reads the state from a string
    ~State();
    inline void setData(Data *d){
      data=d;
    }
    inline void setdataused(){
	usecounts=false; uselengths=false; usesums=false;
	if(datainference==INFDATA_COUNTS || datainference==INFDATA_ALLNOTLENGTHS || datainference==INFDATA_ALL) usecounts=true;
	if(datainference==INFDATA_LENGTHS || datainference==INFDATA_ALL) uselengths=true;
	if(datainference==INFDATA_TOTALLENGTHS || datainference==INFDATA_ALLNOTLENGTHS || datainference==INFDATA_ALL) usesums=true;
    }
// ***** Functions to do with merging and splitting
    void merge(int a,int b);// Merges a and b
    //void split(int a);// Splits a: with random allocation
    double splitSAMS(int i, int j,bool greedy, State * relstate=NULL);/// performs the SAMS split on the individuals i and j, which must be in the same population
    double probOfSplitSAMS(int a, int b, State * relstate=NULL);///< Calculates the probability of the named split having been performed by a merge.

    void splitSAMS(int a, State * relstate);/// splits a population using the SAMS algorithm
    void splitSAMSgreedy(int a,int cmax, State * relstate=NULL);/// splits a population using the SAMS algorithm, accepting only the highest possible split
    int addEmptyPop();///<Adds a new population, updating beta,popX and popsize, returning the index of the population
    //int movePcaAutoInd(int testi,int i, int j,PcaData *pca=NULL);/// Moves ind if it is "close enough" to i or j
//    void movePcaAuto(int opindex,int i, int j,double targd,PcaData *pca);
//    bool testPcaAuto(int a, int b,int i, int j,double targd,PcaData *pca);
    bool isObtuseAngle(double d1, double d2, double d3);///< For calculating whether a point is "between" two others
    
// ***** Functions to do with the posterior evaluation
    void setIntendedSplit(int i,int j);
    double posteriorProb(State * relstate=NULL);///< evaluates the posterior probability of a state
    double posteriorSetProb(int a);///evaluates the posterior of a single set
    double posteriorSetProb(int a,State * relstate);///evaluates the posterior of a single set
    double posteriorSetProbPartial(int a,State * relstate=NULL);///<evaluates a partial set probability, exclucing set rem

    double posteriorLengthSetProb(int a);///<evaluates the posterior probability of the length data for set a
    double posteriorLengthSetProb(int a, int b);///<evaluates the posterior probability of the length data for pop a copying from pop b

// ***** Alternative posteriors
    double posteriorIndSetProb(int a);///< evaluates the posterior of a single set, treating each ind as a pop but with a prior similar to above
// ***** Functions to do with the prior for chunk counts
    inline double sumPriorX(int i, int pop)
    {
	if((int)copyPrior[i].size()<data->getDim()) {std::cout<<"sumPriorX:"<<copyPrior[i].size()<<"<"<<data->getDim()<<"!"<<std::endl;throw(std::string("copyprior is of wrong size!  Is it initialized?"));}
	double s=0;
	std::vector<int> jlist=getIndInPop(pop);
	for(unsigned int j=0;j<jlist.size();j++) s+=copyPrior[i][jlist[j]];
	return(s);
    }
    inline double calcSumPriorXab(int a, int b)
    {
	if((int)copyPrior.size()<data->getDim()) {std::cout<<"sumPriorXab:"<<copyPrior.size()<<"<"<<data->getDim()<<"!"<<std::endl;throw(std::string("copyprior is of wrong size!  Is it initialized?"));}
	double s=0;
	std::vector<int> ilist=getIndInPop(a);
	for(unsigned int i=0;i<ilist.size();i++) {s+=sumPriorX(ilist[i],b);}
	return(s);
    }
    inline double getBetaV(int a)
    {
	return(sumXa(a)/(data->allsum()/corfactor+1.0));
    }
    inline double getBetaF(int a)
    {
	if(betamodel==BETAMOD_F2 || betamodel==BETAMOD_F2_COPYMAT) return(betaF[0]);
	else return(betaF[a]);
    }
    inline void setBetaF(int a,double newd)
    {
	if(betamodel==BETAMOD_F2|| betamodel==BETAMOD_F2_COPYMAT) betaF[0]=newd;
	else betaF[a]=newd;
    }
    inline int getHyperParamLength(){
	if(betamodel==BETAMOD_F2|| betamodel==BETAMOD_F2_COPYMAT) return(1);
	else if (betamodel==BETAMOD_F) return(delta.size());
	else return(0);
    }
    inline double getDelta(int a) {
	if(delta.size()==0) return(-1);// not a used parameter
	if(betamodel==BETAMOD_F2|| betamodel==BETAMOD_F2_COPYMAT) return(delta[0]);
	else return(delta[a]);
    }
    inline void printDelta(std::ostream * out) {
	if(delta.size()>0) {
		for(unsigned int i=0;i<delta.size()-1;i++) *out<<delta.at(i)<<",";
		*out<<delta.at(delta.size()-1);
	}
    }
    inline void printBetaF(std::ostream * out) {
	if(betaF.size()>0) {
	for(unsigned int i=0;i<betaF.size()-1;i++) *out<<betaF.at(i)<<",";
 	*out<<betaF.at(betaF.size()-1);
	}
    }
    inline void setDelta(int a,double newd) {
	if(betamodel==BETAMOD_F2|| betamodel==BETAMOD_F2_COPYMAT) delta[0]=newd;
	else delta[a]=newd;
    }
    inline double getBeta(int a,int b)
    {
	if(betamodel==BETAMOD_EQUI) {// equipartition: 1/K
	  if(beta.size()==0) throw(std::string("Incorrect beta size!"));
	  return(1.0/beta.size());
	}else if(betamodel==BETAMOD_CONST) {
	  if(b>=0) return(sumbeta/beta.size());
     	  else throw(std::string("Error in getBeta!"));
	}else if(betamodel==BETAMOD_F || betamodel==BETAMOD_F2 || betamodel==BETAMOD_F2_COPYMAT) {
	  double beta;
//	std::cout<<"a="<<a<<" b="<<b<<" beta="<<"(1.0-"<<getBetaF(b)<<")/"<<getBetaF(b)<<")*"<<getBetaV(b)<<std::endl;
//std::cout<<"geta="<<getBetaV(b)<<"*"<<(1.0-getBetaF(a))<<"/"<<getBetaF(a)<<std::endl;
	  if(a!=b) beta=getBetaV(b)*((1.0-getBetaF(a))/getBetaF(a)) * (getN()/(getN()-1.0));
 	  else beta=getBetaV(b)*(1.0+getDelta(a))*((1.0-getBetaF(a))/getBetaF(a))*((psize[a]-1.0)/psize[a]);
	  /* *** MAJOR VERSION CHANGE NOTE: the above line should read:
 	  else beta=getBetaV(b)*(1.0+getDelta(a))*((1.0-getBetaF(a))/getBetaF(a))*((psize[a]-1.0)/psize[a])* (getN()/(getN()-1.0));
	  // However, this gets absorbed into delta and doesn't importantly affect inference. 
	  I'm therefore saving this bug for a potential major version change, as it will affect backwards compatibility.*/
	  if(betamodel==BETAMOD_F2_COPYMAT) beta+=calcSumPriorXab(a,b);
//std::cout<<"beta["<<a<<","<<b<<"]="<<getBetaV(b)*(1.0+getDelta(a))*(1.0-getBetaF(a))/getBetaF(a)<<" + "<<calcSumPriorXab(a,b)<<std::endl;
	  return(beta);
	}else if(betamodel==BETAMOD_COPYMAT) {
	  return(calcSumPriorXab(a,b));
	}else throw(std::string("Error in getBeta: invalid beta model"));
    }
   inline double getSumBeta(int a=-1){
	if(betamodel==BETAMOD_EQUI) {// equipartition: 1/K
	  return(1.0);
	}else if(betamodel==BETAMOD_CONST) {
	  return(sumbeta);
	}else if(betamodel==BETAMOD_F || betamodel==BETAMOD_F2 || betamodel==BETAMOD_F2_COPYMAT) {
	  double ret=0.0;
	  for(unsigned int b=0;b<psize.size();b++) ret+=getBeta(a,b);
	  return(ret);
	}else if(betamodel==BETAMOD_COPYMAT) {
	  double ret=0.0;
	  for(unsigned int b=0;b<psize.size();b++) ret+=calcSumPriorXab(a,b);
	  return(ret);
	}else throw(std::string("Error in getBeta: invalid beta model"));
    }
    inline std::vector<double> getBetaVector(int a){
	std::vector<double> b;
	for(unsigned int c1=0;c1<beta.size();c1++) b.push_back(getBeta(a,c1));
    return(b);}
    inline double getAlpha(){return(alpha);}
    inline std::vector<double> getHyperPrior(){return(hyperprior);}
    inline int getBetaModel(){return betamodel;}
    inline void assigncopyPrior(std::vector<std::vector<double> > *vin){
	if((betamodel==BETAMOD_COPYMAT || betamodel==BETAMOD_F2_COPYMAT) && (int)vin->size()==getN()) {
	copyPrior=*vin;
	std::cout<<"copyprior.size="<<copyPrior.size()<<std::endl;
	}else throw(std::string("Invalid assignCopyPrior matrix!"));
    }///< Assigns the prior model with a specified vector

// ***** IO

    void setVerbose(bool nv){verbose=nv;};///<Sets verbosity
    void print(std::ostream * out,bool chat=true);
    void printBeta(std::ostream * out);
    void printX(std::ostream * out,bool perindiv);
    void setprint(std::ostream * out);///< Prints in a bracketed format
    inline void iterPrint(std::ostream * out){
	*out<<"<Iteration>"<<std::endl;
	setprint(out);
	*out<<"</Iteration>"<<std::endl;
    }///< Prints an iteration fronted set

// ***** IO from strings/outputfiles

    std::vector<double> readOutputVectorD(std::string str);///< Returns a double vector from an output string
    inline void setBetaFromString(std::string res){
	beta=readOutputVectorD(res);
    }
    inline void setAlphaFromString(std::string res){
	std::vector<double> tmp=readOutputVectorD(res);
	alpha=tmp[0];
	if(tmp.size()>1)throw(std::string("Error: too many alphas specified!"));
    }
    inline void setDeltaFromString(std::string res){
	delta=readOutputVectorD(res);
    }
    inline void setFFromString(std::string res){
	betaF=readOutputVectorD(res);
    }

// ***** Functions for accessing individuals and populations

    inline int getPop(int i)
    {
	if(i<0|| i>=(int)ind.size()) {std::cerr<<"Individual "<<i<<" not recognised in getPop!"<<std::endl;throw(std::string("Invalid individual requested in getPop"));}
        return ind[i];
    } ///<Get population of the individual specified
    inline void setInd(int i,unsigned int c)
    {
	if((int) c>=getP()) throw(std::string("Invalid population"));
        ind[i]=c;
    } ///<Set accessor to the individual population
    inline int getP()
    {
        return psize.size();
    }///<Returns the number of populations
    inline int getPsize(int i)
    {
        return psize[i];
    }///<Returns the number of individuals in a population
    inline int getPlength(int i)
    {
        return indinp[i].size();
    }///<Returns the number of superindividuals in a population
    inline int getN(){return(data->getN());}///< gets the number of individuals
    inline int getDim(){return(data->getDim());}///< gets the number of superindividuals

    void moveInd(int i,int popto);///<Moves individual i to population popto
    std::vector<int> getIndInPop(int i);///< Gets all the individuals in a population
    void setIndInPop(std::vector<int> v,int i);///< Assigns a vector of individuals to a population

// ***** Functions for accessing data in the count matrix

    inline std::vector<std::vector<double> >*  getPopCopyMatrix() {
 	return(&popX);
    }///<Returns a pointer to the population copying matrix
    inline double X(int i, int j,bool cf=true) {
	if(cf) return(data->get(i,j)/corfactor);
	else return(data->get(i,j));
    }/// Convenience accessor for the data
    double sumXab(int a, int b,bool cf=true) ;/// Sum of copy counts between two populations
    double calcSumXab(int a, int b,bool cf=true);/// CALCULATE Sum of copy counts between two populations
    double sumX(int i, int pop,bool cf=true);///<returns the sum of the number of copies to individual i from population pop
    double sumY(int pop, int i,bool cf=true);///<returns the sum of the number of copies from population pop to individual i 
    double sumLx(int i, int pop);///<returns the sum of the length of copies to individual i from population pop
    double sumLy(int pop, int i);///<returns the sum of the length of copies from population pop to individual i 
    double sumlogLpix(int i, int pop);///<returns the sum of the log pi = x*log(l*x)
    double sumlogLpiy(int pop, int i);///<returns the sum of the length of copies from population pop to individual i 
    double sumlogLgammax(int i, int pop);///<returns the sum of the log gammas
    double sumlogLgammay(int pop, int i);///<returns the sum of the log gammas
    inline double sumXa(int a) {
	std::vector<int> ilist=getIndInPop(a);
	double sXa=0.0;
	for(unsigned int i=0;i<ilist.size();i++) sXa+=sumXall(ilist[i]);
	return(sXa);
    }///<Sum of the copies to a population a
    inline double sumXall(int i) {
	return(data->rowsum(i)/corfactor);
    }/// returns the sum of the number of copies from to individual i from all pops
    inline void testsumXab(){
        bool test=true;
  for(unsigned int c1=0;c1<popX.size();c1++)for(unsigned int c2=0;c2<popX[c1].size();c2++) if(sumXab(c1,c2)!=calcSumXab(c1,c2)) test=false;
      if(!test) {
  for(unsigned int c1=0;c1<popX.size();c1++)for(unsigned int c2=0;c2<popX[c1].size();c2++) cout<<"["<<c1<<","<<c2<<"]="<<sumXab(c1,c2)<<" - "<<calcSumXab(c1,c2)<<endl;
	throw(string("sumXab invalid!"));
      }
    }
// ******** Sampling Distributions

    double sampleBetaP(int param, int modelpart);///< samples from the hyperprior Gamma distributions (0=delta, 1=F,2-5=L priors) from the 3 current models (0=chunks, 1=lengths, 2=mean sizes)
    double LGammaDistProb(double alpha,double beta, double y);///<LogGamma distribution probability
    double LDirichletProb(std::vector<double> prior,std::vector<double> post);///<Probability of the post vector when the parameters are "prior"

// ******** Utility functions

    void mergeOnCoincidence(std::vector< std::vector<double> > *vecin,bool mergerule, double mergeval);///< Creates a state consistent with the pairwise coincidence and the merge rule
    void addCoincidence(std::vector< std::vector<double> > *vecin);///< Adds this state to the coincidence matrix at *vecin
    double getDistanceSq(std::vector< std::vector<double> > *meanX);///< Calculates the distance from this state to meanX
    double popDist(std::vector< std::vector<double> > *meanX,int pop,int indiv,double penalty=1.0);///<The distance change when moving indiv out of pop
    double minDistanceMove(std::vector< std::vector<double> > *meanX,double penalty=1.0);///<Performs a single state update on the distance, moving to the one that minimises it
    inline std::vector<int> allIndInOrder(){
	std::vector<int> allindiv=getIndInPop(0);
	for(int c1=1;c1<getP();c1++) {
		std::vector<int> tmp=getIndInPop(c1);
		for(unsigned int c2=0;c2<tmp.size();c2++) {
			allindiv.push_back(tmp[c2]);
		}
	}
	return(allindiv);
    }///< Gets the list of individuals, in order of populations
    void reorderPop(int from, int to);///< reorder the population list
    void reorderIndiv(int pop,int from, int to);///< reorder the list of individuals within a population
 
// ******** Things to do with the likelihood extraction
    double getApproxLogLikelihood();///< Get the approximated log-likelihod using the MOM estimator of P_{ab}
    double getApproxLogLikelihood(int a);///< Get the approximated log-likelihod for a set
    
// ******** Things to do with the admixture model

    double admixtureSetLogLikelihood(int a,int b,std::vector<std::vector<double> > * Q,std::vector<std::vector<double> > * P, std::vector<double> *qColSums);///< Loglikelihood of a given set
    double admixtureLogLikelihood(std::vector<std::vector<double> > * Q,std::vector<std::vector<double> > * P, std::vector<double> *qColSums);///< Loglikelihood under the admixture model
    double admixtureLogLikelihoodIndiv(int a, std::vector<std::vector<double> > * Q,std::vector<std::vector<double> > * P, std::vector<double> *qColSums);///<Loglikelihood under the admixture model for a single individual
    double admixtureApproxSetLogLikelihood(int set,std::vector<std::vector<double> > * Q,std::vector<std::vector<double> > * P);///< Approximate ll for all indivs in a pop
    double admixtureApproxIndSetLogLikelihood(int ind,int pop,std::vector<std::vector<double> > * Q,std::vector<std::vector<double> > * P);///< Approximate ll for an indiv from a pop
    double admixtureApproxIndLogLikelihood(int ind,std::vector<std::vector<double> > * Q,std::vector<std::vector<double> > * P);///< Approximate ll for an indiv
    inline double calcColSum(std::vector<std::vector<double> > *M,int a){
	double r=0;
	for(unsigned int c1=0;c1<M->size();c1++) r+=M->at(c1)[a];
	return(r);
    }

// ***** Super individuals

    void setSuperIndivRule(bool super);///< If setSuperIndivRule =true, diag is removed from the likelihood via the superindividual "max" self copying rule
    
    bool checkPopForIgnoredSuper(int a);///< returns true if there is a super-individual who should be ignored in the population
// ***** FUNCTIONS TO DO WITH THE LENGTH MATRIX
    inline double L(int i, int j) {
	return(data->get(i,j) * dlength->get(i,j));
    }/// Convenience accessor for the length data
    inline double logLpi(int i, int j) {
	if(i==j || L(i,j)==0) return(0);
	return((data->get(i,j)-1.0) * log(L(i,j)));
    }/// Convenience accessor for the length data
    inline double sumLab(int a, int b){
	return(popL[a][b]);
    }/// Convenience accessor for the length data
    inline double sumlogLpi(int a, int b) {
	return(poplogLpi[a][b]);
    }/// Convenience accessor for the length data
    inline double sumlogLgamma(int a, int b) {
	return(poplogLgamma[a][b]);
    }/// Convenience accessor for the length data
    inline double calcSumLab(int a, int b){
	vector<int> pa=getIndInPop(a);
	vector<int> pb=getIndInPop(b);
	double ret=0.0;
	for(unsigned int c1=0;c1<pa.size();c1++)for(unsigned int c2=0;c2<pb.size();c2++)ret+=L(pa[c1],pb[c2]);
	return(ret);
    }/// Convenience accessor for the length data
    inline double calcSumlogLpi(int a, int b) {
	vector<int> pa=getIndInPop(a);
	vector<int> pb=getIndInPop(b);
	double ret=0.0;
	for(unsigned int c1=0;c1<pa.size();c1++)for(unsigned int c2=0;c2<pb.size();c2++)ret+=logLpi(pa[c1],pb[c2]);
	return(ret);
    }/// Convenience accessor for the length data
    inline double calcSumlogLgamma(int a, int b) {
	vector<int> pa=getIndInPop(a);
	vector<int> pb=getIndInPop(b);
	double ret=0.0;
	for(unsigned int c1=0;c1<pa.size();c1++)for(unsigned int c2=0;c2<pb.size();c2++){if(L(pa[c1],pb[c2])!=0){ cout<<"R"<<data->get(pa[c1],pb[c2])<<endl;  ret+=mylgamma(data->get(pa[c1],pb[c2])); }};
	return(ret);
    }/// Convenience accessor for the length data

// ***** FUNCTIONS TO DO WITH THE LENGTH MATRIX PRIOR
    inline double getLalpha(int a, int b) {
	if(a!=b) return(priorLengths[0]);
	return(priorLengths[0]*(1+priorLengths[2]));
    }
    inline double getLbeta(int a, int b) {
	if(a!=b) return(priorLengths[1]);
	return(priorLengths[1]*(1+priorLengths[3]));
    }
    inline std::vector<double> getLengthHyperPrior(){
	if(hyperprior.size()<NUMHYPERPARAMLENGTH*2) throw(string("Error in getLengthHyperPrior: too short!"));
	int endoffset=0;
	if(usesums) endoffset=NUMHYPERPARAMLENGTH;
	std::vector<double> ret(hyperprior.end()-NUMHYPERPARAMLENGTH*2-endoffset,hyperprior.end()-endoffset);
	return(ret);
    }
    inline double getpriorLengthsParam(int i){
	return(priorLengths[i]);
    }
    inline void setpriorLengthsParam(int i,double s){
	priorLengths[i]=s;
    }

// Functions to do with the average length vector
    
    inline double calcMu(int i) {
	//double mu=0.0;
//	for(int j=0;j<data->getN();j++) {mu+=data->get(i,j) * dlength->get(i,j);}
//	return(mu/data->rowsum(i));
//	for(int j=0;j<data->getN();j++) {mu+=data->get(i,j);}
	return(data->allsum()/data->rowsum(i));
    }
    inline double mu(int i) {
	if(i<0 || i>= (int)muvec.size()) throw(string("Incorrect mu entry requested!"));
	return(muvec[i]);
    }
/*    inline double mualpha(int i=0){
	if(i>0) throw(string("Request for Unimplemented mualpha>0"));
	return(0.000001);
    }
    inline double mubeta(int i=0){
	if(i>0) throw(string("Request for Unimplemented mualpha>0"));
	return(0.000001);
    }
    inline double effcountsprior(){
	return(effcounts);
    }*/

    double posteriorNumChunksSetProb(int a);

    inline std::vector<double> getTotalChunksHyperPrior(){
	if(hyperprior.size()<NUMHYPERPARAMTOTLENGTH) throw(string("Error in getTotalChunksHyperPrior: too short!"));
	int endoffset=0;
	std::vector<double> ret(hyperprior.end()-NUMHYPERPARAMTOTLENGTH*2-endoffset,hyperprior.end()-endoffset);
	return(ret);
    }
    inline double getpriorNumChunksParam(int i){
	if(i<0 || i>= (int)priorMeanSizes.size()) {cerr<<"i="<<i<<endl;throw(string("Invalid parameter requested!"));}
	return(priorMeanSizes[i]);
    }
    inline void setpriorNumChunksParam(int i,double s){
	if(i<0 || i>= (int)priorMeanSizes.size()) {cerr<<"i="<<i<<endl;throw(string("Invalid parameter requested!"));}
	priorMeanSizes[i]=s;
    }

    inline int getModelType(){
      return(modeltype);
    }
protected:
    Data *data, *dlength;///< the data. dlength is the chunk lengths data
    std::vector<int> ind;///<population of individual
    std::vector<int> psize;///<number of individuals in a population
    std::vector<double> beta;///<For betamodels 1 and 2 this is simply beta.  It is not used for model 3 or 4
    std::vector<double> betaF;///< Only for model 3, this is the F vector
    std::vector<double> delta;///< Only for model 3, this is the self copying vector
    std::vector<std::vector<double> > copyPrior;///< Only for model 5, this is the prior copying matrix
    std::vector<std::vector<double> > popX;///< data of copies for the populations
    std::vector<std::vector<double> > popL;///< data of copies for the populations
    std::vector<std::vector<double> > poplogLpi;///< data of copies for the populations
    std::vector<std::vector<double> > poplogLgamma;///< data of copies for the populations
    std::vector<std::vector<int> > indinp;///< list of individuals in a given population
    std::vector<double> hyperprior;///< Hyperprior parameters : either beta, or the vector of (k_F,theta_F,k_delta,theta_delta) for the F model.  Optionally, also (k_alpha0,k_beta0,k_alpha,k_beta,beta_alpha0,beta_beta0,beta_alpha,beta_beta) when lengths are added, and (k_alphamu, k_betamu, beta_alphamu,beta_betamu) for the average chunk lengths
    std::vector<double> priorLengths;///<priors for the chunk length model (alpha_0,beta_0,delta_alpha,delta_beta)
    std::vector<double> muvec;///< vector of average lengths 
    std::vector<double> priorMeanSizes;///<priors for the mean chunk length (alpha_mu,beta_mu,gamma_mu)

    double sumbeta;
    double alpha;///<The parameter for equality in set sizes
    int betamodel;///< The model for beta
    double corfactor;///< The denominator for the likelihood to correct for correlations (2.0 corrects for copying to implying copying from; appropriate for ordinary mcmc; 1.0 is appropriate for admixture since we have fixed populations in this case)

    std::vector<double> diagmod;///< Diagonal modifier; each population has reduced diagonal as self copying within accepted populations is forbidden
    std::vector<double> epopsize;///< Effective population size; each population has a reduced effective population size on the diagonal

    double mylgamma(double z);///<Needed for the prior
    void rpermute(std::vector<int> * in);///< Changes *in to a random permutation of *in
    void rpermute(std::vector<int>* in1,std::vector<int>* in2);///< Changes *in1 & *in2 (same length) to *the same* random permutation of themselves
    void copyState(State *in);///<Copies a state into this one
    void removePop(int lose,int keep=-1);///< Removes a population and places individuals in a another specified population
    bool verbose;
    inline std::string removeNumbers(std::string s) {
	size_t found=s.find_first_of("0123456789");
	while (found!=std::string::npos)
  	{
    	  s.erase(found,1);
	  found=s.find_first_of("0123456789",found);
	}
	return(s);
    }
    int datainference;
    bool usecounts,uselengths,usesums;
    int mergepreva,mergeprevb,mergenewab;
    int modeltype;
};

} // end namespace fines
#endif
