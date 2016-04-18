#include "pcadata.h"

using namespace std;
namespace fines
{

PcaData::PcaData(vector<vector<double> > *rawdata,bool calcsvd,bool calcdistances){
  this->rawdata=rawdata;
  if(calcsvd){
    svdOfData();
    if(calcdistances)calcDistances();
  }   
}

void PcaData::svdOfData()
{
  int n=rawdata->size();

  MatrixXf A = MatrixXf(n,n);
  for(int c1=0;c1<n;c1++){
    if((long)rawdata->at(c1).size()!=n) throw(string("Error creating SVD: Data is not square"));
    for(int c2=0;c2<n;c2++) A(c1,c2) = rawdata->at(c1)[c2];
  }
  
// no need to deal with continents , they are normalised already
// Move all the zeroes to the last column
  for(int c1=0;c1<n;c1++){
    double tmean=0;
    for(int c2=0;c2<n;c2++){
      tmean+=A(c1,c2);
    }
    tmean/=(n-1);
    //double tmp=A(c1,c1); // should be 0 if not a continent
    A(c1,c1)=tmean;
    //A(c1,c1)=A(c1,n-1);
    //A(c1,n-1)=tmp;
  }
  /*for(int c1=0;c1<n;c1++){for(int c2=0;c2<n;c2++){
    if(c2!=c1)A(c1,c1)+=A(c1,c2)/(double)(n-1);
  }}*/
// subtract row means
  vector<double> rowmeans(n,0);
  for(int c1=0;c1<n;c1++){for(int c2=0;c2<n;c2++){
    rowmeans[c1]+=A(c1,c2); // *** Warning: should this be c1 or c2?
  }}
/*  cout<<"RM:";
  for(int c1=0;c1<n;c1++)cout<<rowmeans[c1]<<",";
  cout<<endl;*/
  for(int c1=0;c1<n;c1++){
    //if(A(c1,n-1)>0)  rowmeans[c1]= rowmeans[c1]/(double)n;
    if(A(c1,c1)>0)  rowmeans[c1]= rowmeans[c1]/(double)n;
    else rowmeans[c1]= rowmeans[c1]/(double)(n-1);
  }
  for(int c1=0;c1<n;c1++){for(int c2=0;c2<n;c2++){
    if(A(c1,c2)>0) A(c1,c2)-=rowmeans[c1]; // *** Warning: should this be c1 again? It should match above!
  }}
/*  cout<<"Name="<<std::string(inputdata->getName(0).mb_str())<<endl;
  cout<<"F:";
  for(int c1=0;c1<n;c1++){for(int c2=0;c2<n;c2++){
    if(c1==0)cout<< A(c1,c2)<<",";
  }if(c1==0)cout<< endl;}*/
  
  MatrixXf tA=A.transpose();
  MatrixXf AtA= A * tA; // Simon : best to look for eigenvectors in this matrix (not tA * A)
// Calculate SVD  
  Eigen::EigenSolver<MatrixXf>  ES(AtA);
  MatrixXf AV = ES.pseudoEigenvectors().transpose();
  MatrixXf Aeval = ES.pseudoEigenvalueMatrix ();
  svdM.clear();
  for(int c1=0;c1<n;c1++){svdM.push_back(vector<double>(n,0));}
  for(int c1=0;c1<n;c1++){
    for(int c2=0;c2<n;c2++){
    svdM[c1][c2]=AV(c1,c2);
  }}
  eval=vector<double>(n,0);
// sort
  vector<int> neworder;
  for(int c1=0;c1<n;c1++) neworder.push_back(c1);

  for(int c1=1;c1<n;c1++) {
    double val=Aeval(neworder[c1],neworder[c1]);
    int c2=c1-1;
    while(c2>=0){
      if(Aeval(neworder[c2],neworder[c2])>val){
	neworder[c2+1]=neworder[c2];
	c2--;
      }else break;
    }
    neworder[c2+1]=c1;
  }
// assign  
  for(int c1=0;c1<n;c1++) {
    eval[c1]=Aeval(neworder[n-c1-1],neworder[n-c1-1]);
    for(int c2=0;c2<n;c2++){
      svdM[neworder[n-c1-1]][neworder[n-c2-1]]=AV(neworder[n-c1-1],neworder[n-c2-1]);
    }
  }

}

void PcaData::calcDistances(){
  if(getSize()==0) return;
  distances = vector<vector<double> >(svdM.size(),vector<double>(svdM.size(),0));
  for(unsigned int c1=0;c1<distances.size();c1++){
//    distances[c1]= vector<double>(svdM.size(),0);
    for(unsigned int c2=0;c2<distances[c1].size();c2++){
      for(unsigned int c3=0;c3<svdM[c1].size();c3++){
      //distances[c1][c2] += eval[c3] * eval[c3] * (svdM[c3][c1] -svdM[c3][c2])*(svdM[c3][c1] -svdM[c3][c2]);
      if(eval[c3]>1) distances[c1][c2] += eval[c3] * (svdM[c3][c1] -svdM[c3][c2])*(svdM[c3][c1] -svdM[c3][c2]);
      }
      distances[c1][c2] = sqrt(distances[c1][c2]);      
    }
  }
}

}