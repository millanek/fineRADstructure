#ifndef __RNG_H__
#define __RNG_H__

#include <gsl/gsl_rng.h>
#include <ctime>
#include<iostream>
#include <cstring>
#include <fstream>
#include <vector>

#define ABS(x) ((x) >= 0 ? (x) : -(x))
#define myMin(a,b) ((a) <= (b) ? (a) : (b))
#define myMax(a,b) ((a) >= (b) ? (a) : (b))
#define OVERFLO 1e100
#define UNDERFLO 1e-100


namespace fines
{

extern gsl_rng * rng;

unsigned long makerng(bool fast=false);
unsigned long seedrng(unsigned long seed=0,bool verbose=true);
int saverng(std::string fname);
int loadrng(std::string fname);
void freerng();

double RandomReal(double low, double high);
int RandomInteger(int low, int high);
double rnd();
double RGamma(double n,double lambda);
void RDirichlet(const std::vector<double> * a, std::vector<double>  * b);
long RPoisson(double mu);
double RExpon(double av);
double RNormal(double mu,double sd) ;
double fsign( double num, double sign );
double sexpo(void);
 double snorm();
 double genexp(double av);   
 long ignpoi(double mean);  
 long ignuin(int low, int high);   
 double genunf(double low, double high);   
 long   Binomial(int n, double p);
 long   Binomial1(int n, double p);
 double BinoProb(int n, double p,int i);
 void LogRDirichlet (const double *a, const int k, double *b,double *c);


} // end namespace fines
#endif // __RNG_H__
