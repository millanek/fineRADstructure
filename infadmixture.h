#ifndef INFADMIXTURE_H
#define INFADMIXTURE_H
#include <vector>
#include <iomanip>
#include "state.h"
#include "fines.h"
#include "rng.h"
#include "math.h"

namespace fines
{

/**
    @brief Inference algorithm for the admixture "bolt on" model
*/
class InfAdmixture
{
public:
    InfAdmixture(Data *d,State *s,double alpha=1.0,bool v=false,double initialq=0.001,bool test=false);///<Constructor if everything already known
    ~InfAdmixture();

    inline void resetCounters(){
	numQs=0;numPs=0;
	accQs=0;accPs=0;
    }///<Resets the acceptance rate counters 
    void samplePs(bool mean=true);///< samples the Ps either with the mean or from the correct Dirichet
    void initialQs(double spread);///< creates the initial Q's from the current state, smearing an amount spread from the correct population over the rest
    void printPs(std::ostream * out);///< Prints the p;s
    void printQs(std::ostream * out);///< Prints the Q;s

    void exportXmlHead(std::ostream * fileout);
    void exportXmlTail(std::ostream * fileout);
    void exportXmlIter(std::ostream * fileout,int iter);
    inline void exportXmlComment(std::ostream * fileout,std::string comment){
	    *fileout<< "<comment>"<<comment<<"</comment>"<<std::endl;
    };///<Print a comment tag  
    void metropolis(int prevints,int numints, int thin=1,std::ostream * fileout=NULL);///< Runs the mcmc on the admixture model
    int moveProw(int a);///< Moves a single row of P by MH
    int moveQrow(int a);///< Moves a single row of Q by MH
    int movePs();
    int moveQs();
    inline std::vector<double> getAlphaVec(){
	return(alphavec);
    }
    inline double getAlpha(){return(alphavec[0]);}
    inline void setAlpha(double a){
	while(alphavec.size()<P.size())alphavec.push_back(a);
	for(unsigned int c1=0;c1<P.size();c1++) alphavec[c1]=a;
    }
    double priorForQ(int a);///< Gets the prior for Q row a
    double priorForP(int a);///< Gets the prior for P row a

    void merge(int a, int b);///< Merges two populations
    inline std::vector<std::vector<double> > * getP(){return(&P);};
    inline std::vector<std::vector<double> > * getQ(){return(&Q);};
    inline std::vector<double> *  getQsums(){return(&QcolSums);};
    void test();///<Runs tests for Q
protected:
    bool atest;
    Data *data;
    State *state;
    std::vector<std::vector<double> > P;
    std::vector<std::vector<double> > Q;
    std::vector<double> QcolSums;
    std::vector<double> alphavec;
    bool verbose;
    double logstateprob;
    long numQs,numPs;
    long accQs,accPs;
    double P_QS;
};

} // end namespace fines
#endif
