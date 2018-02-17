#include "rng.h"
#include <ctime>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

using namespace std;
namespace fines {

gsl_rng * rng;

unsigned long makerng(bool fast)
{
    const gsl_rng_type *rng_type;
    gsl_rng_env_setup();
    rng_type = gsl_rng_default;
    rng = gsl_rng_alloc (rng_type);
    unsigned long seed=0;//Just to create a valid rng
    if(fast) {gsl_rng_set(rng,seed);
	return((unsigned long)seed);
    }else return(seedrng(-1));
    
}

unsigned long seedrng(unsigned long seed,bool verbose)
{
    unsigned long tseed;
    FILE *devrandom;
    if(seed==0) {
    if ((devrandom = fopen("/dev/random","r")) == NULL)
    {
        tseed = (unsigned long) time(NULL);
        if(verbose) printf("Got seed %lu from time()\n",tseed);
    }
    else
    {
        size_t fread_res=fread(&tseed,sizeof(tseed),1,devrandom);
// Note: on e.g. cygwin, sizeof(tseed) != fread_res, but the RNG is initiated correctly.
// This is possibly an unsigned long conversion issue? 
        if (fread_res != sizeof(tseed) && verbose) {cerr<<"Warning: RNG read of "<<tseed<<" is not of expected size. This is usually not a problem and occurs in some environments including cygwin."<<endl;}
	if(verbose) printf("Got seed %lu from /dev/random\n",tseed);
        fclose(devrandom);
    }
    }else {tseed=seed;
	if(verbose) printf("Using specified seed %lu\n",tseed);
    }
    gsl_rng_set(rng,tseed);
    return(tseed);
}

void freerng()
{
  gsl_rng_free(rng);
}

int saverng(std::string fname)
{
    FILE *stream;
    if( (stream = fopen(fname.c_str(),"w")) == NULL)
    {
	return(1);//error;
    }
    int res=gsl_rng_fwrite (stream, rng);
    fclose(stream);
    return(res);
}

int loadrng(std::string fname)
{
    const gsl_rng_type *rng_type;
    gsl_rng_env_setup();
    rng_type = gsl_rng_default;
    rng = gsl_rng_alloc (rng_type);
    FILE *stream;

    if( (stream = fopen(fname.c_str(),"r")) == NULL)
    {
	makerng();
	return(1);//error;
    }
    int res=gsl_rng_fread (stream, rng);
    if(res!=0){//no success
	makerng();
    }
    fclose(stream);
    return(res);
}



double RandomReal(double low, double high)
/* Get a random number between low and high */
{
  return(gsl_rng_uniform(rng)*(high-low)+low);
}
/*-------------------------------------*/
int RandomInteger(int low, int high)
/* Get a random integer between low and high INCLUSIVE*/
{
  return (low + gsl_rng_uniform_int(rng,high+1-low));
}

/*=======================================================*/
/*  Uniform(0,1) random number generation*/

double rnd()
{
  double value;

  do
    value = RandomReal(0.0,1.0);
  while ((value==0.0)||(value==1.0));

  return value;
}
/*-----------Gamma and dirichlet from Matt.----------*/
  /* gamma random generator from Ripley, 1987, P230 */


double RGamma(double n,double lambda)
{
  double aa;
  double w;
  //  int i;

	double x=0.0;
	if(n<1)
	{
		const double E=2.71828182;
		const double b=(n+E)/E;
		double p=0.0;
		one: 
		p=b*rnd();
		if(p>1) goto two;
		x=exp(log(p)/n);
		if(x>-log(rnd())) goto one;
		goto three;
		two: 
		x=-log((b-p)/n);
		if (((n-1)*log(x))<log(rnd())) goto one;
		three:;	
	}
	else if(n==1.0)

	  /* exponential random variable, from Ripley, 1987, P230  */	
	{
		double a=0.0;
		double u,u0,ustar;
	ten:
		u=rnd();
		u0=u;
	twenty:
		ustar=rnd();
		if(u<ustar) goto thirty;
		u=rnd();
		if(u<ustar) goto twenty;
		a++;
		goto ten;
	thirty:
		return (a+u0)/lambda;
	}
	else
	{
		double static nprev=0.0;
		double static c1=0.0;
		double static c2=0.0;
		double static c3=0.0;
		double static c4=0.0;
		double static c5=0.0;
		double u1;
		double u2;
		if(n!=nprev)
		{
			c1=n-1.0;
			aa=1.0/c1;
			c2=aa*(n-1/(6*n));
			c3=2*aa;
			c4=c3+2;
			if(n>2.5) c5=1/sqrt(n);
		}
		four:
		u1=rnd();
		u2=rnd();
		if(n<=2.5) goto five;
		u1=u2+c5*(1-1.86*u1);
		if ((u1<=0) || (u1>=1)) goto four;
		five:
		w=c2*u2/u1;
		if(c3*u1+w+1.0/w < c4) goto six;
		if(c3*log(u1)-log(w)+w >=1) goto four;
		six:
		x=c1*w;		
		nprev=n;
	}	

	return x/lambda;
}


/*
double
LogRGamma (double n, double lambda)
{
  //double aa;
  //  double w;
  //  int i;
  double logx;
  //  return log(RGamma(n, lambda));
  if (n < 1)
  //this is the case we need to worry about underflow
  //copied code from down below but work with logx
  //instead of x
    {
      const double E = 2.71828182;
      const double b = (n + E) / E;
      double p = 0.0;
    one:
      p = b * rnd ();
      if (p > 1)
        goto two;
      logx =  log (p) / n;
      if (logx > log(-log (rnd ())))
        goto one;
      goto three;
    two:
      logx = log(-log (b - p)) -log(n);

      if (((n - 1) * logx) < log (rnd ()))
        goto one;
    three:
return logx-log(lambda);
}
else
//otherwise log the standard version 
return log(RGamma(n,lambda));
}*/




//Melissa's version, adapted from an algorithm on wikipedia.  January 08
double LogRGamma(double n, double lambda) {
  double v0, v[3], E=2.71828182, em, logem, lognm;
  int i;
  if (lambda!=1.0) {printf("lambda=%e!\n", lambda); exit(-1);}
  if (n >= 1.0) return log(RGamma(n, lambda));
  v0 = E/(E+n);
  while (1) {
    for (i=0; i<3; i++) v[i] = rnd();
    if (v[0] <= v0) {
      logem = 1.0/n*log(v[1]);
      em = exp(logem);
      lognm = log(v[2])+(n-1)*logem;
    }
    else {
      em = 1.0-log(v[1]);
      logem = log(em);
      lognm = log(v[2]) - em;
    }
    if (lognm <= (n-1)*logem - em)
      return logem - log(lambda);
  }
  }



/*--------------------------------------*/

/* Dirichlet random generator
   a and b are arrays, containing doubles.
   a is the array of parameters
   b is the output array, where b ~ Dirichlet(a)  
  modified for C++ by Dan Lawson
   */

void RDirichlet(const std::vector<double> * a, std::vector<double>  * b)
{
unsigned int i,k=a->size();
	if(b->size()!=k) throw(string("Invalid a and b sizes in RDirichlet!"));
	double sum=0.0;
	for(i=0;i<k;i++)
	{
		b->at(i)=RGamma(a->at(i),1);
		sum += b->at(i);
	}
	for(i=0;i<k;i++)
	{
		b->at(i) /= sum;
	}
}


/*This function returns both a logged and unlogged version
of the dirichlet function. Designed to cope with
underflows in the RGamma function.
made by Daniel
b is the output array and c is a logged version of b*/

void
LogRDirichlet (const double *a, const int k, double *b,double *c)
{
  int i;
  double sum = 0.0;
  double sum2;
  for (i = 0; i < k; i++)
    {
      c[i] = LogRGamma (a[i], 1);
      b[i]=exp(c[i]);
      sum += b[i];
    }
  
    /* patch added May 2007 to set gene frequencies equal if all draws from the Gamma distribution are very low. Ensures that P and logP remain defined in this rare event */
  if(sum<UNDERFLO){
    for(i=0;i<k;i++){
      b[i] = 1.0/(double)(k);
      c[i] = log(b[i]);
    }
    }else{
    sum2=log(sum);
    for (i = 0; i < k; i++)
      {
	c[i]-=sum2;
	b[i]/=sum;
      }
    }
  
}


/*---------------------------------------*/

long RPoisson(double mu)
/*
**********************************************************************
     long RPoissondouble mu)
                    GENerate POIsson random deviate
                              Function
     Generates a single random deviate from a Poisson
     distribution with mean AV.
                              Arguments
     av --> The mean of the Poisson distribution from which
            a random deviate is to be generated.
     RExpon <-- The random deviate.
                              Method
     Renames KPOIS from TOMS as slightly modified by BWB to use RANF
     instead of SUNIF.

     ----substituted rnd for ranf--JKP, 11/98------

     For details see:
               Ahrens, J.H. and Dieter, U.
               Computer Generation of Poisson Deviates
               From Modified Normal Distributions.
               ACM Trans. Math. Software, 8, 2
               (June 1982),163-179
**********************************************************************
**********************************************************************
                                                                      
                                                                      
     P O I S S O N  DISTRIBUTION                                      
                                                                      
                                                                      
**********************************************************************
**********************************************************************
                                                                      
     FOR DETAILS SEE:                                                 
                                                                      
               AHRENS, J.H. AND DIETER, U.                            
               COMPUTER GENERATION OF POISSON DEVIATES                
               FROM MODIFIED NORMAL DISTRIBUTIONS.                    
               ACM TRANS. MATH. SOFTWARE, 8,2 (JUNE 1982), 163 - 179. 
                                                                      
     (SLIGHTLY MODIFIED VERSION OF THE PROGRAM IN THE ABOVE ARTICLE)  
                                                                      
**********************************************************************
      INTEGER FUNCTION RPOISSONIR,MU)
     INPUT:  IR=CURRENT STATE OF BASIC RANDOM NUMBER GENERATOR
             MU=MEAN MU OF THE POISSON DISTRIBUTION
     OUTPUT: IGNPOI=SAMPLE FROM THE POISSON-(MU)-DISTRIBUTION
     MUPREV=PREVIOUS MU, MUOLD=MU AT LAST EXECUTION OF STEP P OR B.
     TABLES: COEFFICIENTS A0-A7 FOR STEP F. FACTORIALS FACT
     COEFFICIENTS A(K) - FOR PX = FK*V*V*SUM(A(K)*V**K)-DEL
     SEPARATION OF CASES A AND B
*/
{
extern double fsign( double num, double sign );
static double a0 = -0.5;
static double a1 = 0.3333333;
static double a2 = -0.2500068;
static double a3 = 0.2000118;
static double a4 = -0.1661269;
static double a5 = 0.1421878;
static double a6 = -0.1384794;
static double a7 = 0.125006;
static double muold = 0.0;
static double muprev = 0.0;
static double fact[10] = {
    1.0,1.0,2.0,6.0,24.0,120.0,720.0,5040.0,40320.0,362880.0
};
static long ignpoi,j,k,kflag,l,m;
static double b1,b2,c,c0,c1,c2,c3,d,del,difmuk,e,fk,fx,fy,g,omega,p,p0,px,py,q,s,
    t,u,v,x,xx,pp[35];

    if(mu == muprev) goto S10;
    if(mu < 10.0) goto S120;
/*
     C A S E  A. (RECALCULATION OF S,D,L IF MU HAS CHANGED)
*/
    muprev = mu;
    s = sqrt(mu);
    d = 6.0*mu*mu;
/*
             THE POISSON PROBABILITIES PK EXCEED THE DISCRETE NORMAL
             PROBABILITIES FK WHENEVER K >= M(MU). L=IFIX(MU-1.1484)
             IS AN UPPER BOUND TO M(MU) FOR ALL MU >= 10 .
*/
    l = (long) (mu-1.1484);
S10:
/*
     STEP N. NORMAL SAMPLE - SNORM(IR) FOR STANDARD NORMAL DEVIATE
*/
    g = mu+s*snorm();
    if(g < 0.0) goto S20;
    ignpoi = (long) (g);
/*
     STEP I. IMMEDIATE ACCEPTANCE IF IGNPOI IS LARGE ENOUGH
*/
    if(ignpoi >= l) return ignpoi;
/*
     STEP S. SQUEEZE ACCEPTANCE - Srnd(IR) FOR (0,1)-SAMPLE U
*/
    fk = (double)ignpoi;
    difmuk = mu-fk;
    u = rnd();  /*was ranf -- JKP*/
    if(d*u >= difmuk*difmuk*difmuk) return ignpoi;
S20:
/*
     STEP P. PREPARATIONS FOR STEPS Q AND H.
             (RECALCULATIONS OF PARAMETERS IF NECESSARY)
             .3989423=(2*PI)**(-.5)  .416667E-1=1./24.  .1428571=1./7.
             THE QUANTITIES B1, B2, C3, C2, C1, C0 ARE FOR THE HERMITE
             APPROXIMATIONS TO THE DISCRETE NORMAL PROBABILITIES FK.
             C=.1069/MU GUARANTEES MAJORIZATION BY THE 'HAT'-FUNCTION.
*/
    if(mu == muold) goto S30;
    muold = mu;
    omega = 0.3989423/s;
    b1 = 4.166667E-2/mu;
    b2 = 0.3*b1*b1;
    c3 = 0.1428571*b1*b2;
    c2 = b2-15.0*c3;
    c1 = b1-6.0*b2+45.0*c3;
    c0 = 1.0-b1+3.0*b2-15.0*c3;
    c = 0.1069/mu;
S30:
    if(g < 0.0) goto S50;
/*
             'SUBROUTINE' F IS CALLED (KFLAG=0 FOR CORRECT RETURN)
*/
    kflag = 0;
    goto S70;
S40:
/*
     STEP Q. QUOTIENT ACCEPTANCE (RARE CASE)
*/
    if(fy-u*fy <= py*exp(px-fx)) return ignpoi;
S50:
/*
     STEP E. EXPONENTIAL SAMPLE - SEXPO(IR) FOR STANDARD EXPONENTIAL
             DEVIATE E AND SAMPLE T FROM THE LAPLACE 'HAT'
             (IF T <= -.6744 THEN PK < FK FOR ALL MU >= 10.)
*/
    e = sexpo();
    u = rnd();  /*was ranf--JKP*/
    u += (u-1.0);
    t = 1.8+fsign(e,u);
    if(t <= -0.6744) goto S50;
    ignpoi = (long) (mu+s*t);
    fk = (double)ignpoi;
    difmuk = mu-fk;
/*
             'SUBROUTINE' F IS CALLED (KFLAG=1 FOR CORRECT RETURN)
*/
    kflag = 1;
    goto S70;
S60:
/*
     STEP H. HAT ACCEPTANCE (E IS REPEATED ON REJECTION)
*/
    if(c*fabs(u) > py*exp(px+e)-fy*exp(fx+e)) goto S50;
    return ignpoi;
S70:
/*
     STEP F. 'SUBROUTINE' F. CALCULATION OF PX,PY,FX,FY.
             CASE IGNPOI .LT. 10 USES FACTORIALS FROM TABLE FACT
*/
    if(ignpoi >= 10) goto S80;
    px = -mu;
    py = pow(mu,(double)ignpoi)/ *(fact+ignpoi);
    goto S110;
S80:
/*
             CASE IGNPOI .GE. 10 USES POLYNOMIAL APPROXIMATION
             A0-A7 FOR ACCURACY WHEN ADVISABLE
             .8333333E-1=1./12.  .3989423=(2*PI)**(-.5)
*/
    del = 8.333333E-2/fk;
    del -= (4.8*del*del*del);
    v = difmuk/fk;
    if(fabs(v) <= 0.25) goto S90;
    px = fk*log(1.0+v)-difmuk-del;
    goto S100;
S90:
    px = fk*v*v*(((((((a7*v+a6)*v+a5)*v+a4)*v+a3)*v+a2)*v+a1)*v+a0)-del;
S100:
    py = 0.3989423/sqrt(fk);
S110:
    x = (0.5-difmuk)/s;
    xx = x*x;
    fx = -0.5*xx;
    fy = omega*(((c3*xx+c2)*xx+c1)*xx+c0);
    if(kflag <= 0) goto S40;
    goto S60;
S120:
/*
     C A S E  B. (START NEW TABLE AND CALCULATE P0 IF NECESSARY)
*/
    muprev = 0.0;
    if(mu == muold) goto S130;
    muold = mu;
    m = myMax(1L,(long) (mu));
    l = 0;
    p = exp(-mu);
    q = p0 = p;
S130:
/*
     STEP U. UNIFORM SAMPLE FOR INVERSION METHOD
*/
    u = rnd();  /*was ranf here-- JKP*/
    ignpoi = 0;
    if(u <= p0) return ignpoi;
/*
     STEP T. TABLE COMPARISON UNTIL THE END PP(L) OF THE
             PP-TABLE OF CUMULATIVE POISSON PROBABILITIES
             (0.458=PP(9) FOR MU=10)
*/
    if(l == 0) goto S150;
    j = 1;
    if(u > 0.458) j = myMin(l,m);
    for(k=j; k<=l; k++) {
        if(u <= *(pp+k-1)) goto S180;
    }
    if(l == 35) goto S130;
S150:
/*
     STEP C. CREATION OF NEW POISSON PROBABILITIES P
             AND THEIR CUMULATIVES Q=PP(K)
*/
    l += 1;
    for(k=l; k<=35; k++) {
        p = p*mu/(double)k;
        q += p;
        *(pp+k-1) = q;
        if(u <= q) goto S170;
    }
    l = 35;
    goto S130;
S170:
    l = k;
S180:
    ignpoi = k;
    return ignpoi;
}

/*-----------------------------------*/
double RNormal(double mu,double sd) 
     /* Returns Normal rv with mean mu, variance sigsq.   
        Uses snorm function of Brown and Lovato.  By JKP*/
{

  return (mu + sd*snorm());

}
/*
**********************************************************************
                                                                      
                                                                      
     (STANDARD-)  N O R M A L  DISTRIBUTION                           
                                                                      
                                                                      
**********************************************************************
**********************************************************************
                                                                      
     FOR DETAILS SEE:                                                 
                                                                      
               AHRENS, J.H. AND DIETER, U.                            
               EXTENSIONS OF FORSYTHE'S METHOD FOR RANDOM             
               SAMPLING FROM THE NORMAL DISTRIBUTION.                 
               MATH. COMPUT., 27,124 (OCT. 1973), 927 - 937.          
                                                                      
     ALL STATEMENT NUMBERS CORRESPOND TO THE STEPS OF ALGORITHM 'FL'  
     (M=5) IN THE ABOVE PAPER     (SLIGHTLY MODIFIED IMPLEMENTATION)  
                                                                      
     Modified by Barry W. Brown, Feb 3, 1988 to use RANF instead of   
     SUNIF.  The argument IR thus goes away.                          
                                                                      
**********************************************************************
     THE DEFINITIONS OF THE CONSTANTS A(K), D(K), T(K) AND
     H(K) ARE ACCORDING TO THE ABOVEMENTIONED ARTICLE
*/
double snorm()    /*was snorm(void) -- JKP*/
{
static double a[32] = {
    0.0,3.917609E-2,7.841241E-2,0.11777,0.1573107,0.1970991,0.2372021,0.2776904,
    0.3186394,0.36013,0.4022501,0.4450965,0.4887764,0.5334097,0.5791322,
    0.626099,0.6744898,0.7245144,0.7764218,0.8305109,0.8871466,0.9467818,
    1.00999,1.077516,1.150349,1.229859,1.318011,1.417797,1.534121,1.67594,
    1.862732,2.153875
};
static double d[31] = {
    0.0,0.0,0.0,0.0,0.0,0.2636843,0.2425085,0.2255674,0.2116342,0.1999243,
    0.1899108,0.1812252,0.1736014,0.1668419,0.1607967,0.1553497,0.1504094,
    0.1459026,0.14177,0.1379632,0.1344418,0.1311722,0.128126,0.1252791,
    0.1226109,0.1201036,0.1177417,0.1155119,0.1134023,0.1114027,0.1095039
};
static double t[31] = {
    7.673828E-4,2.30687E-3,3.860618E-3,5.438454E-3,7.0507E-3,8.708396E-3,
    1.042357E-2,1.220953E-2,1.408125E-2,1.605579E-2,1.81529E-2,2.039573E-2,
    2.281177E-2,2.543407E-2,2.830296E-2,3.146822E-2,3.499233E-2,3.895483E-2,
    4.345878E-2,4.864035E-2,5.468334E-2,6.184222E-2,7.047983E-2,8.113195E-2,
    9.462444E-2,0.1123001,0.136498,0.1716886,0.2276241,0.330498,0.5847031
};
static double h[31] = {
    3.920617E-2,3.932705E-2,3.951E-2,3.975703E-2,4.007093E-2,4.045533E-2,
    4.091481E-2,4.145507E-2,4.208311E-2,4.280748E-2,4.363863E-2,4.458932E-2,
    4.567523E-2,4.691571E-2,4.833487E-2,4.996298E-2,5.183859E-2,5.401138E-2,
    5.654656E-2,5.95313E-2,6.308489E-2,6.737503E-2,7.264544E-2,7.926471E-2,
    8.781922E-2,9.930398E-2,0.11556,0.1404344,0.1836142,0.2790016,0.7010474
};
static long i;
static double snorm,u,s,ustar,aa,w,y,tt;
    u = rnd();   /* was ranf--JKP*/
    s = 0.0;
    if(u > 0.5) s = 1.0;
    u += (u-s);
    u = 32.0*u;
    i = (long) (u);
    if(i == 32) i = 31;
    if(i == 0) goto S100;
/*
                                START CENTER
*/
    ustar = u-(double)i;
    aa = *(a+i-1);
S40:
    if(ustar <= *(t+i-1)) goto S60;
    w = (ustar-*(t+i-1))**(h+i-1);
S50:
/*
                                EXIT   (BOTH CASES)
*/
    y = aa+w;
    snorm = y;
    if(s == 1.0) snorm = -y;
    return snorm;
S60:
/*
                                CENTER CONTINUED
*/
    u = rnd();                /*was ranf--JKP*/
    w = u*(*(a+i)-aa);
    tt = (0.5*w+aa)*w;
    goto S80;
S70:
    tt = u;
    ustar = rnd();                /*was ranf--JKP*/
S80:
    if(ustar > tt) goto S50;
    u = rnd();               /*was ranf--JKP*/
    if(ustar >= u) goto S70;
    ustar = rnd();               /*was ranf--JKP*/
    goto S40;
S100:
/*
                                START TAIL
*/
    i = 6;
    aa = *(a+31);
    goto S120;
S110:
    aa += *(d+i-1);
    i += 1;
S120:
    u += u;
    if(u < 1.0) goto S110;
    u -= 1.0;
S140:
    w = u**(d+i-1);
    tt = (0.5*w+aa)*w;
    goto S160;
S150:
    tt = u;
S160:
    ustar = rnd();               /*was ranf--JKP*/
    if(ustar > tt) goto S50;
    u = rnd();               /*was ranf--JKP*/
    if(ustar >= u) goto S150;
    u = rnd();               /*was ranf--JKP*/
    goto S140;
}

/*
**********************************************************************
     double RExpon(double av)
                    GENerate EXPonential random deviate
                              Function
     Generates a single random deviate from an exponential
     distribution with mean AV.
                              Arguments
     av --> The mean of the exponential distribution from which
            a random deviate is to be generated.
                              Method
     Renames SEXPO from TOMS as slightly modified by BWB to use RANF
     instead of SUNIF.
     For details see:
               Ahrens, J.H. and Dieter, U.
               Computer Methods for Sampling From the
               Exponential and Normal Distributions.
               Comm. ACM, 15,10 (Oct. 1972), 873 - 882.
**********************************************************************
*/
double RExpon(double av)
{
static double RExpon;

    RExpon = sexpo()*av;
    return RExpon;
}

/*
**********************************************************************
                                                                      
                                                                      
     (STANDARD-)  E X P O N E N T I A L   DISTRIBUTION                
                                                                      
                                                                      
**********************************************************************
**********************************************************************
                                                                      
     FOR DETAILS SEE:                                                 
                                                                      
               AHRENS, J.H. AND DIETER, U.                            
               COMPUTER METHODS FOR SAMPLING FROM THE                 
               EXPONENTIAL AND NORMAL DISTRIBUTIONS.                  
               COMM. ACM, 15,10 (OCT. 1972), 873 - 882.               
                                                                      
     ALL STATEMENT NUMBERS CORRESPOND TO THE STEPS OF ALGORITHM       
     'SA' IN THE ABOVE PAPER (SLIGHTLY MODIFIED IMPLEMENTATION)       
                                                                      
     Modified by Barry W. Brown, Feb 3, 1988 to use RANF instead of   
     SUNIF.  The argument IR thus goes away.                          
                                                                      
**********************************************************************
     Q(N) = SUM(ALOG(2.0)**K/K!)    K=1,..,N ,      THE HIGHEST N
     (HERE 8) IS DETERMINED BY Q(N)=1.0 WITHIN STANDARD PRECISION
*/
double sexpo(void)
{
static double q[8] = {
    0.6931472,0.9333737,0.9888778,0.9984959,0.9998293,0.9999833,0.9999986,1.0
};
static long i;
static double sexpo,a,u,ustar,umin;
static double *q1 = q;
    a = 0.0;
    u = rnd();  /* was ranf--JKP */
    goto S30;
S20:
    a += *q1;
S30:
    u += u;
    if(u <= 1.0) goto S20;
    u -= 1.0;
    if(u > *q1) goto S60;
    sexpo = a+u;
    return sexpo;
S60:
    i = 1;
    ustar = rnd();
    umin = ustar;
S70:
    ustar = rnd();  /* was ranf--JKP */
    if(ustar < umin) umin = ustar;
    i += 1;
    if(u > *(q+i-1)) goto S70;
    sexpo = a+umin**q1;
    return sexpo;
}

/*------------------------------------*/
double fsign( double num, double sign )
/* Transfers sign of argument sign to argument num */
{
if ( ( sign>0.0f && num<0.0f ) || ( sign<0.0f && num>0.0f ) )
    return -num;
else return num;
}

/*------------------------------------*/
double genexp(double av)
{
  return RExpon(av);
}
/*------------------------------------*/
long ignpoi(double mean)
{
  return RPoisson(mean);
}
/*------------------------------------*/
long ignuin(int low, int high)
{
  return RandomInteger(low,high);
}
/*-------------------------------------*/
double genunf(double low, double high)
{
  return RandomReal(low,high);
}
/*-------------------------------------*/
long Binomial(int n, double p)
/*returns a binomial random number, for the number of successes in n trials
  with prob of sucess p.  There's probably a qicker algorithm than this, but I
  can't see how to write the cumulative prob in a simple form*/
{
  int i,sofar;

  sofar = 0;
  for (i=0; i<n; i++)
    if (rnd() < p) sofar++;
  return sofar;
  
}
/*-------------------------------------*/
long Binomial1(int n, double p)
/*returns a binomial random number, for the number of successes in n
trials with prob of sucess p.  There's probably a qicker algorithm
than this, but I can't see how to write the cumulative prob in a
simple form.  This more complicated algorithm, which involves summing
the probabilities appears to be substantially slower than the
simple-minded approach, above.*/
{
  double cum = 0.0;
  int up,down; 
  //  double upvalue,downvalue;
  double rv;
  //  double q = 1 - p;

  if (p<=0.0) return 0;  /*trivial cases*/
  if (p>=1.0) return 0;
  if (n<1) return 0;
  
  rv = rnd();            /*random number in (0,1)*/
  up = n*p;              /*start at mean and work out, adding probs to the total (cum)*/
  down = up;
  
  do
    {
      if (up <= n)
	{
	  cum += BinoProb(n,p,up);
	  if (rv <= cum) return up;
	  up++;
	}
      down--;
      if (down >= 0)
	{	  
	  cum += BinoProb(n,p,down);
	  if (rv <= cum) return down;
	}
    }
  while ((up <=n ) || (down >= 1));

  return Binomial(n,p);  /*in case of reaching no result...possibly due to underflow(?)*/
}
/*-------------------------------------*/
double BinoProb(int n, double p,int i)
/*returns the prob of i successes in n trials with prob of sucess p.*/
{

  double logsum = 0.0;
  double runningtotal = 1.0;
  int j;

  if (i>(n-i))  /*figure out the n-choose-i part*/
    {
      for (j=2; j <= (n-i); j++)
	{
	  runningtotal /= j;
	  if (runningtotal<UNDERFLO)
	    {
	      logsum += log(runningtotal);
	      runningtotal = 1.0;
	    }
	}
      for (j=i+1; j <= n; j++)
	{
	  runningtotal *= j;
	  if (runningtotal>OVERFLO)
	    {
	      logsum += log(runningtotal);
	      runningtotal = 1.0;
	    }
	}
    }
  else
    {
      for (j=2; j <= i; j++)
	{
	  runningtotal /= j;
	  if (runningtotal<UNDERFLO)
	    {
	      logsum += log(runningtotal);
	      runningtotal = 1.0;
	    }
	}
      for (j=n-i+1; j <= n; j++)
	{
	  runningtotal *= j;
	  if (runningtotal>OVERFLO)
	    {
	      logsum += log(runningtotal);
	      runningtotal = 1.0;
	    }
	}
    }
  logsum += log(runningtotal);
  logsum += i*log(p);
  logsum += (n-i)*log(1-p);
  
  return exp(logsum);
}

} // End namespace fines


