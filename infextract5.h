#ifndef INFEXTRACT5_H
#define INFEXTRACT5_H
#include <vector>
#include "state.h"
#include "data.h"
#include "rng.h"
#include "node.h"
#include "fsxml.h"
#include "infadmixture.h"

namespace fines
{

/**
    @brief Inference algorithm: extraction of likelihood distribution
*/
class InfExtract5
{
public:
    InfExtract5(Data *d,FsXml *fs,std::vector<double> bvec,double a,int betamodel=BETAMOD_CONST,double corfactor=1.0,long reps=1,bool v=false);///<Constructor 
    ~InfExtract5();
    inline std::vector<double> getLikelihoods(){return(likelihoods);};///< gets the likelihoods
    double calcLikelihood(std::vector<std::vector<double> > * P);///< Calculate the likelihood
protected:
    Data *data;
    State *state;
    bool verbose;
    std::vector<double> likelihoods;
    double corfactor;
    long counts;
};

} // end namespace fines
#endif
