#ifndef INFEXTRACT4_H
#define INFEXTRACT4_H
#include <vector>
#include "state.h"
#include "data.h"
#include "rng.h"
#include "node.h"
#include "fsxml.h"

namespace fines
{

/**
    @brief Inference algorithm: extraction of posterior distribution
*/
class InfExtract4
{
public:
    InfExtract4(Data *d,FsXml *fs,std::vector<double> bvec,double a,int betamodel=BETAMOD_CONST,double corfactor=1.0,bool v=false);///<Constructor 
    ~InfExtract4();
    inline std::vector<double> getPosteriors(){return(posteriors);};///< gets the posterior probabilities
protected:
    Data *data;
    State *state;
    bool verbose;
    std::vector<double> posteriors;
    long counts;
};

} // end namespace fines
#endif
