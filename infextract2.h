#ifndef INFEXTRACT2_H
#define INFEXTRACT2_H
#include <vector>
#include "state.h"
#include "data.h"
#include "rng.h"
#include "node.h"
#include "fsxml.h"

namespace fines
{

/**
    @brief Inference algorithm: extraction of highest posterior state
*/
class InfExtract2
{
public:
    InfExtract2(Data *d,FsXml *fs,std::vector<double> bvec,double a,int betamodel=BETAMOD_CONST,double corfactor=1.0,bool v=false,Data *dlength=NULL,int datainference=INFDATA_COUNTS,int modeltype=MODELTYPE_FINESTRUCTURE);///<Constructor 
    ~InfExtract2();
    inline State * getState(){return(state);}///< Returns the state

protected:
    Data *data;
    State *state;
    bool verbose;
};

} // end namespace fines
#endif
