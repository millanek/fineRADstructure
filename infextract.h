#ifndef INFEXTRACT_H
#define INFEXTRACT_H
#include <vector>
#include "state.h"
#include "data.h"
#include "rng.h"
#include "node.h"
#include "fsxml.h"

namespace fines
{

/**
    @brief Inference algorithm: extraction of mean state
*/
class InfExtract
{
public:
    InfExtract(Data *d,FsXml *fs,bool v=false);///<Constructor 
    void printMeanX(std::ostream * out);///
    inline std::vector<std::vector<double> > *  getMeanX(){return(&meanX);}
    ~InfExtract();
    inline State * getState(){return(state);}///< Returns the state
    void makeMinSquaresState(double pen=1.0, State * startstate=NULL);///< Creates a minimum distance state from the mean matrix
    void reorder(std::vector<int> allvec);
protected:
    Data *data;
    State *state;
    bool verbose;
    std::vector<std::vector<double> > meanX;
    std::vector<int> order;
};

} // end namespace fines
#endif
