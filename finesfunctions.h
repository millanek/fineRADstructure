#include <cstring>
#include <fstream>
#include "fines.h"
#include "data.h"
#include "inf1.h"
#include "infextract.h"
#include "infextract2.h"
#include "infextract3.h"
#include "infextract4.h"
#include "infextract5.h"
#include "infmcmc.h"
#include "infadmixture.h"
#include "rng.h"
#include "fsxml.h"
#include "fines.h"

#define TREETYPE_USEOBSSTATE 1
#define TREETYPE_USEMERGESTATE 2
#define TREETYPE_USEHILLCLIMBSTATE 3

using namespace std;
namespace fines
{

std::vector<double> getBvec(int betamod,int datainference,double corfactor,vector<double> bvec=vector<double>(0),
			      string betapriorstring=string(""),long ignorelines=0,bool xhead=true,bool yhead=true);
Inf1 mergeTree(int treetype, Data *d, string fs,long testmax,long hcsteps, double corfactor,double betamod,vector<double> bvec,
	       int datainference=INFDATA_COUNTS, int modeltype=MODELTYPE_FINESTRUCTURE,State *startstate=NULL, Data *dlength=NULL, bool havefullxmlinput=true,bool fixK=false,int  treescale=0,bool verbose=false);
Inf1 GlobalReadTree(string filename,Data *d,double alpha,double corfactor,double betamod,vector<double> bvec,
	      int datainference=INFDATA_COUNTS,int modeltype=MODELTYPE_FINESTRUCTURE,bool verbose=false);
InfMCMC GlobalRunMCMC(Data *d, State *initstate,ostream *os,long burnin,long additional,long thinin,string comment,
		      int datainference=INFDATA_COUNTS,Data *dlength=NULL,double pcaprob=0,bool fixK=false,bool verbose=false);
bool getHeader(string filename, double &cval,long &burnin, long &mcmclength,long &cvalstr, string &datafilestr);
int compareDataFiles(string f1, string f2);
} // end namespace fines
