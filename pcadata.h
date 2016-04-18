#ifndef PCADATA_H
#define PCADATA_H
#include <vector>

using namespace std;
#undef min
#undef max
#include <limits>
#include <Eigen/QR>
#include <Eigen/Core>
#include <Eigen/Geometry>

using namespace Eigen;

namespace fines
{

/**
    @brief Eigen vector/value decomposition using Singular Value Decomposition
*/
class PcaData
{
public:
  PcaData(vector<vector<double> > *datin,bool calcsvd=true,bool calcdistances=false);
  ~PcaData(){};
  void svdOfData();
  void calcDistances();
  
  inline vector<vector<double> > getMraw(){
      return(svdM);
  };
  inline vector<vector<double> > *getM(){
      return(&svdM);
  }
  inline double getM(long i, long j){
    if(i<0 || i>= (long)svdM.size()) throw(string("Error: attempting to get illegal i element of SVD matrix"));
    if(j<0 || j>= (long)svdM[i].size()) throw(string("Error: attempting to get illegal j element of SVD matrix"));
    return(svdM[i][j]);
  }
  inline vector<double> *getE(){
    return(&eval);
  }
  inline vector<double> getEraw(){
    return(eval);
  }
  inline double getE(long i){
    if(i<0 || i>= (long)eval.size()) throw(string("Error: attempting to get illegal i element of eigenvalue vector"));
    return(eval[i]);
  }
  inline long getSize(){
    return(eval.size());
  }
  inline double getDist(long i,long j){
    if(i<0 || i>= (long)distances.size()) throw(string("Error: attempting to get illegal i element of distances matrix"));
    if(j<0 || j>= (long)distances[i].size()) throw(string("Error: attempting to get illegal j element of distances matrix"));
    return(distances[i][j]);
  }
protected:
  vector<vector<double> > *rawdata;
  vector<vector<double> > svdM;
  vector<double> eval;
// things to do with distances
  vector<vector<double> > distances;

};

} // end namespace fines


#endif
