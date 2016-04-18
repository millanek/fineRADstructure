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
#include "finesfunctions.h"

using namespace std;
namespace fines
{

ProgramOptions& opt() {
	static ProgramOptions po;	// define a single instance of ProgramOptions per process.
	return po;
}


} // end namespace fines

using namespace fines;
// main function goes outside the weakarg namespace so that
// the linker can find it

// The f model proper has been removed:
//			3/f/F:  F model of Falush et al 2003 (RJ MOVE NOT IMPLEMENTED).
//			(default: -b -0.001,-0.001,-1,-1 for model 3,.

static const char * help=
    "\
    Usage: finestructure [OPTIONS] datafile <initialpopfile> outputfile\n\
    	Datafile is a matrix of copy counts.\n\
	initialpopfile (optional) is a population state e.g. an outputfile.\n\
	outputfile is the destination.\n\
    -m <method>		Method to use.  Default: oMCMC.\n\
			<method> is either MCMCwithTree, oMCMC (MCMC without tree), \n\
			SplitMerge, Tree, or admixture, or a contraction of any.\n\
			Note that admixture does not infer the population\n\
 			so should be provided with a good one in \"initialpopfile\".\n\
    -I <x>          	Initial number of populations.  <x> is either a number\n\
			or \"n\" for the number of individuals, or \"l\" for label detected \n\
			populations.  Default is 1.\n\
    -s <s>		Sets the RNG seed to s (>0)\n\
    -i <i>		Ignores the first i lines of the input file\n\
    -x <num>		Number of burn in iterations for MCMC method.\n\
    -y <num>		Number of sample iterations for MCMC method.\n\
    -z <num>		Thin interval in the output file, for MCMC method.\n\
    -t <num>		Maximum number of tree comparisons for splitting/merging.\n\
    -K                  Fix the number of populations to whatever you started with.\n\
                        This would be set by '-I' or by an initial state file.\n\
    -l <filename>	Specify the average copy length datafile.  -i,-X,-Y options\n\
			*preciding* this file will affect this read; you can set different\n\
			options for the copy rate datafile by specifying these -i,-X,-Y again\n\
			after the -l option.\n\
    -u <datatype>	Use a data inference method; one of :\n\
			counts: use only the copy counts data. (default if -l not specified)\n\
			lengths: use only the copy length data (still needs valid counts data!)\n\
			totallengths: use the mean length of chunk sizes \n\
			all: use all data (careful: this may not be statistically valid).\n\
			default: use counts and totallengths (default with -l specified).\n\
    -a <num>		Set alpha, the prior of the number of parameters\n\
			(default: 1.0).\n\
    -c <num>		Set the likelihood correction factor: L_{used}=L^{1/<corfactor>}.\n\
			(default: 1.0)\n\
    -B <model>		Choose a model for beta:\n\
			1/e/E:	Equipartition model of Pella and Masuda.\n\
			2/c/C:	Constant model.\n\
			4/o/O:  F model of Falush et al 2003 with a single parameter\n\
				for all populations (default).\n\
    -b <num>(,<num>,..)	Hyperparameters for ALL models, in the order COUNTS,LENGTHS,MEANS.  \n\
			COUNTS: *must* be included, even if count matrix not used!\n\
			For model 1, there are no parameters.\n\
			For model 2, set the prior of the distribution of\n\
			population sizes (each population has beta_i=<num>).\n\
			(default: 1.0).\n\
			For model 4, set the hyperprior of the distribution of\n\
			delta and F. Parameters are \n\
			(k_f,k_delta,theta_f,theta_delta) for the parameters of the\n\
			gamma distribution F~Gamma(k_f,theta_f), \n\
			and delta~Gamma(k_delta,theta_delta)\n\
			(default: -b 2,2,0.01,0.01).\n\
			LENGTHS: 8 parameters:\n\
			(k_alpha0,k_beta0,k_alpha,k_beta,beta_alpha0,beta_beta0,beta_alpha,beta_beta)\n\
			MEANS: 6 parameters:\n\
			(k_betamu, k_alphamu, k_kappa, beta_alphamu,beta_betamu,beta_kappa)\n\
			Set K parameters negative for fixed =|k|\n\
			e.g. when finding a tree given the mean parameters.\n\
    -M <modeltype>	Specify the type of inference model for chunk counts.  \n\
			<modeltype> accept contractions and lower case, and can be:\n\
			  1 or Finestructure: standard finestructure model (default).\n\
			  2 or Normalised: Normalise data row and columns within a population.\n\
			  3 or MergeOnly: As 2, but only compare populations being merged or split.\n\
			  4 or Individual: Prior is placed on individual rows instead of \n\
					  population rows. (slowest model).\n\
    -e <name>		Extract details from a state; can be (a unique contraction of):\n\
			beta: the parameter matrix\n\
			X: the copying data matrix for populations\n\
			X2: the normalised copying matrix\n\
			maxstate: maximum observed posterior probability state\n\
			meancoincidence: the mean coincidence matrix\n\
			merge<:value><:split>: create a merge(or split)\n\
			  population from the mean coincidence.\n\
			admixture: gets the population as an admixture matrix.\n\
			Pmatrix: gets the P matrix for the admixture.\n\
  			range:<from>:<to> gets the iterations in the specified range.\n\
			thin:<step>: thins the output by step.\n\
			probability: get the posterior probability of the data\n\
			given the conditions of the outputfile.\n\
			likelihood: samples the likelihood of the data given the conditions\n\
			in the outputfile.\n\
			tree: extract the tree in newick format and print it to a FOURTH file\n\
    -F <filename>	Fix the populations specified in the file.  They should be specified as\n\
			population format, i.e. PopA(ind1,ind2) would fix the data rows ind1 and ind2\n\
			to always be in the same population (they form a 'super individual')\n\
			called PopA. Continents are specified with a * before the name, and are treated\n\
			specially in the tree building phase,  i.e. *ContA(ind1,ind2).  Continents\n\
			are not merged with the rest of the tree.\n\
    -T <type>		When using a merge tree, initialisation can be set to the following:\n\
			1:	Use the initial state \"as is\".\n\
			2:	Perform merging to get to best posterior state.\n\
			3:	Perform full range of moves to to get to best posterior state.\n\
				This is the default.  Set number of attempts with -x <num>.\n\
			4:	As 1, but don't flatten maximum copy rates for the main tree.\n\
			5:	As 2, but don't flatten maximum copy rates for the main tree.\n\
			6:	As 3, but don't flatten maximum copy rates for the main tree.\n\
			7:	As 1, but maximise hyperparameters between merges.\n\
			8:	As 2, but maximise hyperparameters between merges.\n\
			9:	As 3, but maximise hyperparameters between merges.\n\
    -o <name>		File containing a state to use for ordering, if not the main file.\n\
    -k <num>		Change the tree building algorithm.\n\
			0:	Discard all ordering and likelihood information (default).\n\
			1:	Maintain ordering.\n\
			2:	Maintain ordering and likelihood.\n\
    -X			Specifies that there are row names in the data (not necessary for \n\
			ChromoPainter or ChromoCombine style files.)\n\
    -Y			Specifies that there are column names in the data file (as -X, not necessary.)\n\
    -v          	Verbose mode\n\
    -V          	Print Version info\n\
    -h          	This help message\n\
	\n\
    Examples:\n\
    finestructure -X -Y -m omcmc -i 2 -B 4 -b 2,2,0.01,0.01 -s 1 -x 100000 -y 100000 \n\
	-z 1000 datafile.csv resfile.xml\n\
			Infers population structure (-m omcmc) from datafile.csv which\n\
			contains 2 irrelevent lines (-i 2) with row (-X) and column (-Y)\n\
			names, using the F model with a global F and Delta (-B 4) using\n\
			Gamma(2,0.01) distributions. 100000 burn in steps are used (-x)\n\
			and 100000 further iterations are sampled (-y) keeping every\n\
			1000th sample (-z).\n\
    finestructure -X -Y -i 2 -e min datafile.csv resfile.xml resmsfile.xml\n\
			Create a minimum state file from the MCMC output.\n\
    finestructure -X -Y -i 2 -m T -t 20000 -B 4 -b -0.003,-0.94,-1,-1 \n\
	datafile.csv resmsfile.xml restree.xml\n\
			Create a tree (-m T) from the minimum state using the inferred\n\
			values of F (0.003) and Delta (0.94), allowing 20000 (-t 20000).\n\
			different trees to be examined per merge attempt (slow!).\n\
    finestructure -X -Y -i 2 -m admixture -B 4 -b -0.003,-0.94,-1,-1 -x 100000 -y 100000 \n\
	-z 1000 datafile.csv resmsfile.xml admixfile.xml\n\
			Perform admixture (-m admixiture) MCMC using the minimum state \n\
			and parameters found as above.\n\
    finestructure -X -Y -i 2 -e admixture datafile.csv resmsfile.xml admixstate.csv\n\
			Extract the admixture matrix Q for a state in csv format.  This\n\
			is useful for making comparisons to the observed admixture matrix.\n\
    ";

/*    -p <num>		Use the PCA enhanced merge/split proposals.  <num> is between zero and one\n\
			and controls the probability with which PCA proposals are used.\n\
			PCA proposals are much faster but may increase autocorrelation time.\n\
			In most cases, the increased speed is worth it, but it may cause problems\n\
			locating the posterior mode if there are a lot of individuals.  Default: 0.\n\*/
// ********************** PCA MOVE CURRENTLY DOES NOTHING
    
    
string getVersion(){
	string ret;
	//#ifdef PACKAGE_STRING // Omitted due to version matching issues
	//ret.append(PACKAGE_STRING);
	//#else
	ret.append("finestructure");
	//#endif
	ret.append(" build date "); 
	ret.append(__DATE__);
	ret.append(" at ");
	ret.append(__TIME__);
	return(ret);
}

void printVersion(){
	cout<<getVersion()<<endl;
}

std::vector<double> &splitasdouble(const std::string &s, char delim, std::vector<double> &elems) {
    std::stringstream ss(s);
    std::string item;
    while(std::getline(ss, item, delim)) {
    	elems.push_back(atof(item.c_str()));
    }
    return elems;
}

State setState(string fs2,Data * d,std::vector<double> bvec,int betamod, double corfactor,Data *dlength=NULL,int datainference=INFDATA_COUNTS,int modeltype=MODELTYPE_FINESTRUCTURE)
{
	if(opt().verbose) cout<<"Reordering on "<<fs2<<"."<<endl;
	FsXml infile(fs2);
	InfExtract2 iext2r(d,&infile,bvec,1.0,betamod,corfactor,false,dlength,datainference,modeltype);
	return(iext2r.getState());
}

int getModelType(string modelarg){
  std::transform(modelarg.begin(), modelarg.end(),modelarg.begin(), ::toupper);
  if(modelarg.substr(0,1).compare("F")==0 || modelarg.substr(0,1).compare("1")==0) return(MODELTYPE_FINESTRUCTURE);
  if(modelarg.substr(0,1).compare("N")==0 || modelarg.substr(0,1).compare("2")==0) return(MODELTYPE_NORMALISED);
  if(modelarg.substr(0,1).compare("M")==0 || modelarg.substr(0,1).compare("3")==0) return(MODELTYPE_NORMALISEDMERGEONLY);
  if(modelarg.substr(0,1).compare("I")==0 || modelarg.substr(0,1).compare("4")==0) return(MODELTYPE_INDIVIDUAL);
  return(1);
}


int main(int argc, char *argv[])
{
    string comment="Command line: ";
    for(int c1=0;c1< argc;c1++) {comment.append(argv[c1]);comment.append(" ");}
    comment.append("\nVersion: ");
    comment.append(getVersion());
    std::stringstream ss;
    std::string tmp;
    string fs;// The input file
    string fs2;// The input file for ordering
    string fstree;// The *input* file for the tree (only for -e TREE)
    string betapriorstring;// the input file for the beta prior, if betamodel=5
    makerng(true);
    int c;
    char *tmpchar;
    Data *d, *dlength=NULL;
    double alpha=1.0;
    int betamod=BETAMOD_F2;
    double corfactor=-1.0;
    int initpop=-1;
    int ignorelines=0;
    vector<double> bvec; // vector for initialisation of beta
    bool xhead=false,yhead=false;
    unsigned long seed=0;
    bool extract=false;string ext;
    bool havefullxmlinput=false;
    int datainference=INFDATA_ALLNOTLENGTHS;
    int treetype=TREETYPE_USEHILLCLIMBSTATE;
    int modeltype=MODELTYPE_FINESTRUCTURE;
    int treemodification=TREEMOD_FLATTEN;
    bool initpopset=false;
    string fixfile;
    while ((c = getopt (argc, argv, "x:y:z:s:l:u:I:i:m:p:T:XYa:b:M:F:B:e:t:o:c:vKk:hV")) != -1)
        switch (c)
        {
        case('x'):opt().burnin=atoi(optarg);break;
        case('y'):opt().additional=atoi(optarg);break;
        case('z'):opt().thinin=atoi(optarg);break;
        case('s'):seed=strtoul(optarg,NULL,10);break;
	case('l'):dlength=new Data(string(optarg),ignorelines,xhead,yhead);break;
	case('u'):opt().usedata=string(optarg);;break;
        case('I'):if(*optarg!='n' && *optarg!='N') {
			if(*optarg=='l' || *optarg=='L') {initpop=-2;break;}
			initpop=atoi(optarg);
		 }else initpop=-1;
		 initpopset=true;
		 break;
	case('i'):ignorelines=atoi(optarg); break;
	case('e'):extract=true; if(!initpopset) initpop=-1;ext=string(optarg);opt().method.assign("");break;
	case('o'):fs2=string(optarg);break;
	case('p'):opt().pcaprob=atof(optarg);break;
	case('m'):opt().method.assign(optarg);	std::transform(opt().method.begin(), opt().method.end(),opt().method.begin(), ::toupper);break;
	case('M'):modeltype=getModelType(string(optarg));break;
	case('X'):xhead=true;break;
	case('Y'):yhead=true;break;
	case('a'):alpha=atof(optarg);break;
	case('b'):tmp=string(optarg); splitasdouble(tmp,',',bvec);break;
	case('K'):opt().fixK=true;break;
	case('k'):opt().treescale=atoi(optarg);break;
	case('B'):tmp=string(optarg);
		if(tmp.compare("1")==0 || tmp.substr(0,1).compare("e")==0 || tmp.substr(0,1).compare("E")==0) {betamod=BETAMOD_EQUI;
		}else if(tmp.compare("2")==0 || tmp.substr(0,1).compare("c")==0 || tmp.substr(0,1).compare("C")==0) {betamod=BETAMOD_CONST;
		}else if(tmp.compare("3")==0 || tmp.substr(0,1).compare("f")==0 || tmp.substr(0,1).compare("F")==0) {betamod=BETAMOD_F;
		}else if(tmp.compare("4")==0 || tmp.substr(0,1).compare("o")==0 || tmp.substr(0,1).compare("O")==0) {betamod=BETAMOD_F2;
		}else if(tmp.substr(0,1).compare("5")==0 || tmp.substr(0,1).compare("x")==0 || tmp.substr(0,1).compare("X")==0) {betamod=BETAMOD_COPYMAT;tmpchar=strtok(optarg,":"); 
		tmpchar=strtok(NULL,":");if(tmpchar==NULL){ cerr<<"Must specify a datafile for beta matrix with this prior!"<<endl;exit(0);}; betapriorstring=string(tmpchar);
		}else if(tmp.substr(0,1).compare("6")==0 || tmp.substr(0,1).compare("p")==0 || tmp.substr(0,1).compare("P")==0) {betamod=BETAMOD_F2_COPYMAT;tmpchar=strtok(optarg,":"); 
		tmpchar=strtok(NULL,":");if(tmpchar==NULL){ cerr<<"Must specify a datafile for beta matrix with this prior!"<<endl;exit(0);}; betapriorstring=string(tmpchar);}
		break;
	case('T'):treetype=atoi(optarg);
		if(treetype>3 && treetype<7) {treemodification=TREEMOD_NOFLATTEN;treetype-=3;}
		if(treetype>=7) {treemodification=TREEMOD_FLATTENHILLCLIMB;treetype-=6;}
		break;
	case('c'):corfactor=atof(optarg);break;
	case('F'):fixfile=optarg;break;
 	case('t'):opt().test_max=atoi(optarg);break;
        case('v'):opt().verbose=true;break;
	case('h'):cout<<help<<endl;return 0;
        case('V'):printVersion();  return 0;
        case '?':cout<<"Wrong arguments: did not recognise "<<c<<" "<<optarg<<endl<<help<<endl;return 1;
        default:
	    cout<<help<<endl;
            abort ();
        }
    if (argc-optind<1) {cout<<help<<endl;return 0;}
    seed=seedrng(seed,opt().verbose);// 0 means use /dev/random or clock.
    comment.append("\nSeed: ");
    ss<<seed;
    comment.append(ss.str());
// process the main arguments
// read data
    if (argc-optind>=1) {
	if(opt().verbose) cout<<"Opening data file: "<<argv[optind]<<endl;
	try{d=new Data(argv[optind++],ignorelines,xhead,yhead);
	}catch(std::string x){cout<<x<<endl;exit(0);}
    }else {cout<<"Need data first!"<<endl<<help<<endl; return 0;}
    if (argc-optind<1) { 
	cout<<"Need an output file!"<<endl<<help<<endl;return 0;
    }
    if(fixfile.size()>0) {
      if(opt().verbose) cout<<"Using fixed population file "<<fixfile<<endl;
	try{d->makeSuperFromFile(fixfile);
	}catch(std::string x){cout<<x<<endl;exit(0);}
      
    }
    
// assign the corfactor
    if(corfactor>0 && d->getCfactor()<0) {cerr<<"WARNING: You have specified 'C' and provided a datafile that does not contain it."<<endl<<"If you know what you are doing you can ignore this warning."<<endl<<"Otherwise you are advised to use the 'chromocombine' tool to estimate it."<<endl<<"See www.paintmychromosomes.com under 'ChromoCombine' for details."<<endl;
	}else if(d->getCfactor()<0) {cerr<<"WARNING: 'C' NOT READ FROM DATA INPUT FILE."<<endl<<"This means that fineSTRUCTURE may expect to see the wrong variance"<<endl<<"in the data and you will probably experience poor clustering."<<endl<<"You are advised to use the 'chromocombine' tool to estimate 'C'."<<endl<<"See www.paintmychromosomes.com under 'ChromoCombine' for details."<<endl;
	}
    if(corfactor<0 && d->getCfactor()>0) corfactor=d->getCfactor(); // now also 1.0 for everything else
    if(corfactor<0) corfactor=1.0;// backup in case it is not validly set anywhere
      
// decide which data to use
/*    if(dlength==NULL) datainference=INFDATA_COUNTS;// no data for the lengths
    else {*/
	string sdat=opt().usedata;
	std::transform(sdat.begin(), sdat.end(),sdat.begin(), ::toupper);
	if(sdat.substr(0,1).compare("L")==0) datainference=INFDATA_LENGTHS;
	else if(sdat.substr(0,1).compare("C")==0) datainference=INFDATA_COUNTS;
	else if(sdat.substr(0,1).compare("T")==0) datainference=INFDATA_TOTALLENGTHS;
	else if(sdat.substr(0,1).compare("A")==0) datainference=INFDATA_ALL;
	else if(sdat.substr(0,1).compare("D")==0) datainference=INFDATA_ALLNOTLENGTHS;
	else datainference=INFDATA_ALLNOTLENGTHS;
	// otherwise we use all INFDATA_ALLNOTLENGTHS, the default
//    }
// assign the bvec hyperprior for betamod (lengths are added later)
    bvec=getBvec(betamod,datainference,corfactor,bvec,betapriorstring,ignorelines,xhead,yhead);
    
    State * state;
// read state if appropriate, or create a new one
    if (argc-optind>=2) {
	fs=string(argv[optind++]);
	if(opt().verbose) cout<<"Reading state file: "<<fs.c_str()<<endl;
	FsXml infile(fs);
	try{
		if(infile.gotoLineContaining("<Iteration>",true)<0) {// is a population file
			state=new State(d,fs,bvec,alpha,betamod,true,corfactor,dlength,datainference,modeltype);
		}else{// is an xml output file
//cout<<"datasize="<<d->getDim()<<" fn="<<fs<<" ";
//for(unsigned int c1=0;c1<bvec.size();c1++) cout<<bvec[c1]<<",";
//cout<<endl<<" alpha="<<alpha<<" betamod="<<betamod<<" cf="<<corfactor<<" di="<<datainference<<" mt="<<modeltype<<endl;		  
			state=new State(d,&infile,bvec,alpha,betamod,corfactor,false,dlength,datainference,modeltype);
			havefullxmlinput=true;
		}
	;}catch(std::string x){cout<<help<<endl<<"ERROR: "<<x<<endl;exit(0);}

    }else {
	if(opt().verbose) cout<<"Creating state from data."<<endl;
	try{state=new State(d,initpop,bvec,alpha,betamod,corfactor,dlength,datainference,modeltype);}catch(std::string x){cout<<help<<endl<<"ERROR: "<<x<<endl;exit(0);}
    }
    
//
    if(extract && (ext.substr(0,2).compare("TR")==0 )){ // we must have an extra argument, the third line being the tree we READ IN, the 4th being the output file
      if (argc-optind>=2)  {
	fstree=argv[optind++];
	cout<<"Reading tree file "<<fstree<<endl;
      }else{
	cerr<<"For extract -e TREE, need 4 files in the order data mcmc tree output"<<endl;
	exit(0);
      }
    }
// create output file
    filebuf fb;
    try{
    fb.open (argv[optind++],ios::out);
    }catch(std::string x){
	cerr<<"Error opening file!"<<endl<<x<<endl<<help<<endl; return 0;}
    ostream os (&fb);
// Extract only
    if(extract){
	std::transform(ext.begin(), ext.end(),ext.begin(), ::toupper);
	if(ext.substr(0,2).compare("X2")==0) {// normalised copying matrix
		if(fs2.size()>0) setState(fs2,d,bvec,betamod,corfactor,dlength,datainference,modeltype).printX(&os,true);
		else state->printX(&os,true);
	}else if(ext.substr(0,1).compare("X")==0) {// unnormalised copying matrix
		if(fs2.size()>0) setState(fs2,d,bvec,betamod,corfactor,dlength,datainference,modeltype).printX(&os,false);
		else state->printX(&os,false);
	}else if(ext.substr(0,1).compare("B")==0) {// beta matrix
		if(fs2.size()>0) setState(fs2,d,bvec,betamod,corfactor,dlength,datainference,modeltype).printBeta(&os);
		else state->printBeta(&os);
	}else if(ext.substr(0,3).compare("MEA")==0) {// meancoincidence (pairwise) matrix
	FsXml *infile=new FsXml(fs);
	InfExtract iext(d,infile,opt().verbose);
	delete(infile);
	infile=new FsXml(fs);
	InfExtract2 iext2(d,infile,bvec,1.0,betamod,corfactor,false,dlength,datainference,modeltype);
	if(fs2.size()>0){
		State ords=setState(fs2,d,bvec,betamod,corfactor,dlength,datainference,modeltype);
		try{iext.reorder(ords.allIndInOrder());
		}catch(string x){cout<<"ERROR:reorder: "<<x<<endl;exit(0);}
	}
//	iext.getMeanX();
	iext.printMeanX(&os);
	}else if(ext.substr(0,2).compare("MI")==0){// minimum distance state
		FsXml infile(fs);
		size_t firstcolon=ext.find_first_of(":");
		double pen=1.0;
		if(firstcolon!=string::npos) { pen=atof(ext.substr(firstcolon+1).c_str());}
		InfExtract iext(d,&infile,opt().verbose);
		try{state=new State(d,iext.getMeanX(),true, 0.0,bvec,alpha,betamod,corfactor,dlength,datainference,modeltype);}catch(std::string x){cout<<help<<endl<<"ERROR: "<<x<<endl;exit(0);}
		iext.makeMinSquaresState(pen,state);//iext.getState()
		iext.getState()->setprint(&os);	
	}else if(ext.substr(0,3).compare("MER")==0) {
		double mergeval=0.95;
		bool mergerule=true;
		size_t firstcolon=ext.find_first_of(":");
		string opts;
		if(firstcolon!=string::npos) {
			opts=ext.substr(firstcolon+1);
			size_t secondcolon=opts.find_first_of(":");
			mergeval=atof(opts.substr(0,secondcolon).c_str());
			if(secondcolon!=string::npos){
			opts=opts.substr(secondcolon+1);
			if(opts.find_first_of("sS0t")!=string::npos) mergerule=false;
			}
		}
		FsXml infile(fs);
		InfExtract iext(d,&infile,opt().verbose);
		try{state=new State(d,iext.getMeanX(),mergerule, mergeval,bvec,alpha,betamod,corfactor,dlength,datainference,modeltype);}catch(std::string x){cout<<help<<endl<<"ERROR: "<<x<<endl;exit(0);}
		state->setprint(&os);
	}else if(ext.substr(0,3).compare("MAX")==0) {//maximum posterior state
		FsXml infile(fs);
		InfExtract2 iext(d,&infile,bvec,1.0,betamod,corfactor,opt().verbose,dlength,datainference,modeltype);
		state=iext.getState();
		state->setprint(&os);
	}else if(ext.substr(0,1).compare("A")==0 ){// extract admixture matrix
		State * adstate;
		State tmpstate=setState(fs2,d,bvec,betamod,corfactor,dlength,datainference,modeltype);
		if(fs2.size()>0) adstate= new State(&tmpstate);
		else adstate=new State(state);
		InfAdmixture infad(d,adstate,1.0,opt().verbose);
		infad.printQs(&os);
	}else if(ext.substr(0,2).compare("PM")==0 ){// extract P matrix
		State * adstate;
		State tmpstate=setState(fs2,d,bvec,betamod,corfactor,dlength,datainference,modeltype);
		if(fs2.size()>0) adstate= new State(&tmpstate);
		else adstate=new State(state);
		InfAdmixture infad(d,adstate,1.0,opt().verbose);
		infad.printPs(&os);
	}else if(ext.substr(0,2).compare("TH")==0 ){//thin
		size_t firstcolon=ext.find_first_of(":");
		string opts;
		long thin=1;
		if(firstcolon!=string::npos) {
			thin=atoi(ext.substr(firstcolon+1).c_str());
		}
		if(opt().verbose)cout<<"Thinning of "<<fs<<" by "<<thin<<endl;
		FsXml infile(fs);
		infile.thin(thin,&os);
	}else if(ext.substr(0,1).compare("R")==0 ){//extract range
		size_t firstcolon=ext.find_first_of(":");
		string opts;
		long from=0,to=0;
		if(firstcolon!=string::npos) {
			opts=ext.substr(firstcolon+1);
			size_t secondcolon=opts.find_first_of(":");
			from=atof(opts.substr(0,secondcolon).c_str());
			if(secondcolon!=string::npos){
			to=atoi(opts.substr(secondcolon+1,string::npos).c_str());
			}
		}
		if(opt().verbose)cout<<"Range extraction of "<<fs<<" from "<<from<<" to "<<to<<endl;
		FsXml infile(fs);
		infile.rangeextract(from,to,&os);
	}else if(ext.substr(0,1).compare("L")==0 ) {// Likelihood sample
		FsXml infile(fs);
		try{InfExtract5 iext5(d,&infile,bvec,1.0,betamod,corfactor,1,false);
		vector<double> liks=iext5.getLikelihoods();
		for(unsigned long c1=0;c1<liks.size()-1;c1++) os<<setprecision(12)<<liks[c1]<<", ";
 		if(liks.size()>0)os<<setprecision(12)<<liks[liks.size()-1]<<endl;	
		}catch(string x){cerr<<"Error in extraction of likelihoods:"<<x<<endl;exit(0);}
	}else if(ext.substr(0,2).compare("PR")==0 ){// Posterior Probability
		FsXml infile(fs);
		try{InfExtract4 iext4(d,&infile,bvec,1.0,betamod,corfactor,false);
		vector<double> posteriors=iext4.getPosteriors();
		for(unsigned int c1=0;c1<posteriors.size()-1;c1++) os<<setprecision(12)<<posteriors[c1]<<", ";
 		if(posteriors.size()>0)os<<setprecision(12)<<posteriors[posteriors.size()-1]<<endl;	
		}catch(string x){cerr<<"Error in extraction of probabilities:"<<x<<endl;exit(0);}
	}else if(ext.substr(0,2).compare("TR")==0 ){// tree in newick format
	      Inf1 *tmptree = new Inf1(GlobalReadTree(fstree,d,1.0,corfactor,betamod,bvec,
	      datainference,modeltype,opt().verbose));
	      tmptree->printTree(&os,false);
	}else {cerr<<"Error: invalid extraction."<<endl<<help<<endl;}
	// end of extract functions
/// MCMC MODEL
    }else if(opt().method.compare(0,1,"M")==0|| opt().method.compare(0,1,"O")==0) {
      InfMCMC infMCMC=GlobalRunMCMC(d,state,&os,opt().burnin,opt().additional,opt().thinin,comment,datainference,dlength,opt().pcaprob,opt().fixK,opt().verbose);
/*
	try{InfMCMC infMCMC(d,state,dlength,datainference,opt().verbose);
	if(opt().verbose) cout<<"BURN IN PHASE"<<endl;
	infMCMC.metropolis(0,opt().burnin);
	if(opt().verbose) cout<<"MCMC PHASE"<<endl;
	infMCMC.resetCounters();
	infMCMC.metropolis(opt().burnin,opt().additional,opt().thinin,&os);
*/
	try{
	State * state2=new State(infMCMC.getState());
	if(opt().method.compare(0,1,"M")==0 ) {
	  Inf1 inf2(d,state2,dlength,datainference,opt().verbose,opt().test_max,opt().treescale);
	  if(opt().verbose) cout<<"TREE CREATION PHASE"<<endl;
	  try{inf2.mergeHillClimb(&os,false,treemodification);}catch(std::string x){cout<<x<<endl;exit(0);}
	}
	infMCMC.exportXmlTail(&os);
	}catch(std::string x){cout<<"Tree creation error: "<<x<<endl;}
/// (SEMI) DETERMINISTIC MODELS
    }else if (opt().method.compare(0,1,"S")==0) {// split tree
	Inf1 inf1(d,initpop,dlength,datainference,modeltype,opt().verbose,opt().test_max);
	inf1.exportXmlHead(&os,fs,string("SplitTree"),corfactor,opt().burnin);
	inf1.exportXmlComment(&os,comment);
	
	if(initpop>0) {
	  if(opt().verbose) cout<<"SPLIT PHASE"<<endl;
	  try{inf1.splitHillClimb(true);}catch(std::string x){cout<<x<<endl;exit(0);}
	}

	State * state2=new State(inf1.getState());

	Inf1 inf2(d,state2,dlength,datainference,opt().verbose,opt().test_max);
	if(opt().verbose) cout<<"MERGE PHASE"<<endl;
	try{inf2.mergeHillClimb(&os,false,treemodification);}catch(std::string x){cout<<x<<endl;exit(0);}
	inf1.exportXmlTail(&os);
/// MAIN TREE METHOD
    }else if (opt().method.compare(0,1,"T")==0) {// merge tree
      /*
	if(opt().verbose) cout<<"MERGE PHASE"<<endl;
	FsXml *infile=new FsXml(fs);
	InfExtract iext(d,infile,opt().verbose);
	Inf1 * inf1_i=NULL;
	delete(infile);
	infile=new FsXml(fs);
	InfExtract2 * iext2=NULL;
	InfMCMC * infHillClimb=NULL;
	State * state2;
	if(havefullxmlinput){
		iext2=new InfExtract2(d,infile,bvec,1.0,betamod,corfactor,opt().verbose,dlength,datainference);
		if(treetype==1){
			state2=new State(iext2->getState());
			inf1_i=new Inf1(d,iext2->getState(),dlength,datainference,opt().verbose,opt().test_max);
		}else if(treetype==2) {
			inf1_i=new Inf1(d,iext2->getState(),dlength,datainference,opt().verbose,opt().test_max);
			try{inf1_i->mergeHillClimb(NULL,true,false);}catch(std::string x){cout<<x<<endl;exit(0);}
			state2=new State(inf1_i->getState());

		}else if(treetype==3) {
			try{infHillClimb=new InfMCMC(d,iext2->getState(),dlength,datainference,opt().verbose);
			//infHillClimb->hillClimb(0,opt().burnin);
			infHillClimb->metropolis(0,opt().burnin);
			state2=new State(infHillClimb->getState());
			}catch(std::string x){cout<<x<<endl;exit(0);}
		}else {cerr<<"Invalid tree type (-T option)."<<endl;exit(0);}
		delete(infile);
	}else {
		if(treetype==1){
			inf1_i=new Inf1(d,state,dlength,datainference,opt().verbose,opt().test_max);
			state2=new State(state);
		}else if(treetype==2) {
			inf1_i=new Inf1(d,state,dlength,datainference,opt().verbose,opt().test_max);
			try{inf1_i->mergeHillClimb(NULL,true,false);}catch(std::string x){cout<<x<<endl;exit(0);}
			state2=new State(inf1_i->getState());

		}else if(treetype==3){
			try{infHillClimb=new InfMCMC(d,state,dlength,datainference,opt().verbose);
			infHillClimb->hillClimb(0,opt().burnin);
			state2=new State(infHillClimb->getState());
			}catch(std::string x){cout<<x<<endl;exit(0);}
		}else {cerr<<"Invalid tree type (-T option)."<<endl;exit(0);}
	}
	*/
	long newx=-1,newy=-1,newz=-1;
	string olddatafile;
	try{
	double newc=-1;
	if(getHeader(fs,newc,newx,newy,newz,olddatafile)) corfactor=newc; // do it this way to make sure we only overwrite corfactor if it doesn't exist in the xml file
	}catch(string x){
	  cerr<<x<<endl;
	}
	int cpval=0;
	if((cpval=compareDataFiles(olddatafile, d->getFileName()))>0){
	   cerr<<"WARNING! Cannot confirm data file is the same as the MCMC was run on!"<<endl;
	}else if(cpval<0){
	  cerr<<"WARNING!  You are trying to build a tree from a differently named datafile than the one used for the MCMC!  This might be due to running it on a different system or might imply that the file is incorrect. This may result in strange behaviour!"<<endl;
	}
//cout<<"treetype="<<treetype<<" (data) (fs) treetestmax="<<opt().test_max<<" hcs="<<opt().burnin<<" cf="<<corfactor<<" betamod="<<betamod<<" (bvec) ";
//cout<<" di="<<datainference<<" mt="<<modeltype<<" (state) (dlength) hfxml="<<havefullxmlinput<<" fixK="<<opt().fixK<<" ts="<<opt().treescale<<endl;
      Inf1 inf1=mergeTree(treetype,d, fs,opt().test_max,opt().burnin,corfactor,betamod,bvec,
			  datainference,modeltype, state, dlength, havefullxmlinput,opt().fixK,opt().treescale,opt().verbose);
//	Inf1 inf1(state2,dlength,datainference,opt().verbose,opt().test_max);
	inf1.exportXmlHead(&os,fs,string("MergeTree"),corfactor,opt().burnin);
	inf1.exportXmlComment(&os,comment);
    InfMCMC* tmcmc=new InfMCMC(d,inf1.getState(),NULL,INFDATA_COUNTS,0,false);//Used to print the state
	tmcmc->exportXmlIter(&os,0); // export the iteration
//	inf1.getState()->iterPrint(&os);
	try{inf1.mergeHillClimb(NULL,false,treemodification);}catch(std::string x){cout<<x<<endl;exit(0);}
	if(opt().verbose) cout<<"Assigning certainty"<<endl;
	FsXml *infile=new FsXml(fs);
	InfExtract3 iext3(d,infile,inf1.getNodes(),opt().verbose);
	delete(infile);
	/*
	// NOTE: Diagonalise disabled due to errors with force files.
	if(opt().verbose) cout<<"Diagonalise tree"<<endl;
	inf1.diagonaliseOnM(d->getMatrix(),false);
*/
	if(opt().verbose) cout<<"Finish up"<<endl;

/*	if(inf1_i!=NULL) {
		inf1.reorderState(inf1_i->getState());
		inf1_i->getState()->iterPrint(&os);
	}else if(infHillClimb!=NULL){
		inf1.reorderState(infHillClimb->getState());
		infHillClimb->getState()->iterPrint(&os);
	}
*/
	inf1.printTree(&os);
	inf1.exportXmlTail(&os);	
    }else if(opt().method.compare(0,1,"A")==0) {// admixture model
	bool atest=false;
	if(opt().method.compare(0,10,"ADMIXTURET")==0) atest=true;
	if(opt().verbose &&!atest) cout<<"Admixture model..."<<endl;
	else if(opt().verbose &&atest) cout<<"Admixture model test only..."<<endl;
	//state->setprint(&cout);
	InfAdmixture infad(d,state,2,opt().verbose,atest);
//	infad.printPs(&cout);
	try{
	if(opt().verbose) cout<<"BURN IN PHASE"<<endl;
	infad.metropolis(0,opt().burnin);
	infad.exportXmlHead(&os);
	infad.exportXmlComment(&os,comment);
	if(opt().verbose) cout<<"MCMC PHASE"<<endl;
	infad.resetCounters();
	infad.metropolis(opt().burnin,opt().additional,opt().thinin,&os);
	infad.exportXmlTail(&os);
	}catch(std::string x){cerr<<"Error in admixture:"<<endl<<x<<endl;}
    }else {
	cerr<<"Invalid method."<<endl<<help<<endl;
    }

    freerng();
    if(d!=NULL) delete(d);
    fb.close();
    return 0;
}


