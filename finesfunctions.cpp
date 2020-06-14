#include "finesfunctions.h"

namespace fines
{

  std::vector<double> getBvec(int betamod,int datainference,double corfactor,vector<double> bvec,
			      string betapriorstring,long ignorelines,bool xhead,bool yhead) {
    int bvecmin=1;
    if(bvec.size()==0 && betamod==BETAMOD_CONST) bvec.push_back(1.0);
    if(betamod==BETAMOD_F) {
      while(bvec.size()<4) bvec.push_back(-0.001); 
      bvecmin=4;
   }
    else if(betamod==BETAMOD_F2 || betamod==BETAMOD_F2_COPYMAT) {
	while(bvec.size()<2) bvec.push_back(2);
	while(bvec.size()<4) bvec.push_back(0.01);
	bvecmin=4;
    }
    if(betamod==BETAMOD_COPYMAT || betamod==BETAMOD_F2_COPYMAT) {
	Data *d2;
        try{d2=new Data(betapriorstring,ignorelines,xhead,yhead);
	//bvec.clear();//=vector<double>(d2->getN()*d2->getN(),0);
	for(int c1=0;c1<d2->getDim();c1++) {
		    for(int c2=0;c2<d2->getDim();c2++) bvec.push_back(d2->get(c1,c2)/corfactor);
	}
	}catch(string x){cerr<<"Error in betamodel creation:"<<x<<endl;exit(0);}
	bvecmin+=d2->getDim()*d2->getDim();
	delete(d2);	
    }
    int bvecmin2=0;
    if(datainference==INFDATA_LENGTHS || datainference==INFDATA_ALL){
	bvecmin2=NUMHYPERPARAMLENGTH*2;
	while((long)bvec.size()<bvecmin+NUMHYPERPARAMLENGTH) bvec.push_back(2);
	while((long)bvec.size()<bvecmin+NUMHYPERPARAMLENGTH) bvec.push_back(0.01);
    }
    if(datainference==INFDATA_TOTALLENGTHS || datainference==INFDATA_ALLNOTLENGTHS || datainference==INFDATA_ALL){
	while((long)bvec.size()<bvecmin+bvecmin2+NUMHYPERPARAMTOTLENGTH) bvec.push_back(2);
	while((long)bvec.size()<bvecmin+bvecmin2+2*NUMHYPERPARAMTOTLENGTH) bvec.push_back(0.01);
    }
    return(bvec);
}

Inf1 mergeTree(int treetype, Data *d, string fs,long testmax,long hcsteps, double corfactor,double betamod,vector<double> bvec,
	       int datainference,int modeltype, State *startstate, Data *dlength, bool havefullxmlinput,bool fixK,int treescale,bool verbose) {
	if(verbose) cout<<"MERGE PHASE"<<endl;

	FsXml *infile=new FsXml(fs);
//	InfExtract iext(d,infile,verbose);
//	delete(infile);
//	infile=new FsXml(fs);
	Inf1 * inf1_i=NULL;
	InfExtract2 * iext2=NULL;
	InfMCMC * infHillClimb=NULL;
	State * state2;
	if(havefullxmlinput){
		iext2=new InfExtract2(d,infile,bvec,1.0,betamod,corfactor,verbose,dlength,datainference,modeltype);
		if(treetype==TREETYPE_USEOBSSTATE){
			state2=new State(iext2->getState());
			inf1_i=new Inf1(d,iext2->getState(),dlength,datainference,verbose,testmax,treescale);
		}else if(treetype==TREETYPE_USEMERGESTATE) {
			inf1_i=new Inf1(d,iext2->getState(),dlength,datainference,verbose,testmax,treescale);
			try{inf1_i->mergeHillClimb(NULL,true,false);}catch(std::string x){cout<<x<<endl;exit(0);}
			state2=new State(inf1_i->getState());

		}else if(treetype==TREETYPE_USEHILLCLIMBSTATE) {
			try{infHillClimb=new InfMCMC(d,iext2->getState(),dlength,datainference,0,verbose);
			//infHillClimb->hillClimb(0,opt().burnin);
			if(fixK)  infHillClimb->fixK();
			infHillClimb->metropolis(0,hcsteps);
			state2=new State(infHillClimb->getState());
			}catch(std::string x){cout<<x<<endl;exit(0);}
		}else {cerr<<"Invalid tree type (-T option): "<<treetype<<endl;exit(0);}
		delete(infile);
	}else {
		if(treetype==TREETYPE_USEOBSSTATE){
			inf1_i=new Inf1(d,startstate,dlength,datainference,verbose,testmax,treescale);
			state2=new State(startstate);
		}else if(treetype==TREETYPE_USEMERGESTATE) {
			inf1_i=new Inf1(d,startstate,dlength,datainference,verbose,testmax,treescale);
			try{inf1_i->mergeHillClimb(NULL,true,false);}catch(std::string x){cout<<x<<endl;exit(0);}
			state2=new State(inf1_i->getState());

		}else if(treetype==TREETYPE_USEHILLCLIMBSTATE){
			try{infHillClimb=new InfMCMC(d,startstate,dlength,datainference,0,verbose);
			infHillClimb->hillClimb(0,hcsteps);
			state2=new State(infHillClimb->getState());
			}catch(std::string x){cout<<x<<endl;exit(0);}
		}else {cerr<<"Invalid tree type (-t option)."<<endl;exit(0);}
	};
	//	state2->iterPrint(&cout);
	Inf1* inf1=new Inf1(d,state2,dlength,datainference,verbose,testmax,treescale);
	//	inf1.getState()->iterPrint(&cout);
	return(*inf1);
/*	inf1.exportXmlHead(&os);
	inf1.exportXmlComment(&os,comment);
	try{inf1.mergeHillClimb(NULL,false,treesuper);}catch(std::string x){cout<<x<<endl;exit(0);}
	if(opt().verbose) cout<<"Assigning certainty"<<endl;
	infile=new FsXml(fs);
	InfExtract3 iext3(d,infile,inf1.getNodes(),opt().verbose);
	delete(infile);
	if(opt().verbose) cout<<"Diagonalise tree"<<endl;
	inf1.diagonaliseOnM(d->getMatrix(),false);
*/
}

Inf1 GlobalReadTree(string filename,Data *d, double alpha,double corfactor,double betamod,vector<double> bvec,
	       int datainference, int modeltype,bool verbose)
{
  FsXml infile(filename);
  
  Data *d2=NULL;
  string pop=infile.getParam("Pop");
  string newick=infile.getParam("Tree");
//  cout<<"READING TREE:"<<newick<<endl;
  State *state = new State(d,pop,bvec,alpha,betamod,false,corfactor,d2,datainference,modeltype);
  Inf1* ret=new Inf1(d,state,newick,verbose);
  return(*ret);
}

InfMCMC GlobalRunMCMC(Data *d, State *initstate,ostream *os,long burnin,long additional,long thinin,string comment,
		      int datainference,Data *dlength,double pcaprob,bool fixK,bool verbose)
{

  InfMCMC* infMCMC=new InfMCMC(d,initstate,dlength,datainference,pcaprob,verbose);

    if(fixK) infMCMC->fixK();
  try{
    infMCMC->exportXmlHead(os,initstate->getCorfactor(),burnin,additional,thinin);
	infMCMC->exportXmlComment(os,comment);
	if(verbose==1) cout<<"BURN IN PHASE"<<endl;
	infMCMC->metropolis(0,burnin);
	if(verbose==1) cout<<"MCMC PHASE"<<endl;

	infMCMC->resetCounters();	

	for(long c1=0;c1<additional;c1++){infMCMC->metropolis(c1,1,thinin,os,additional);}
	if(additional % thinin==0) infMCMC->exportXmlIter(os,additional);

  }catch(string x){
    cerr<<"Error in GlobalRunMCMC:"<<x<<endl;
    throw(x);
  }	
	//infMCMC.metropolis(0,additional,thinin,os);

  return(*infMCMC);

}

bool getHeader(string filename, double &cval,long &burnin, long &mcmclength,long &cvalstr, string &datafilestr){
    FsXml *infile=new FsXml(filename);
    infile->gotoLineContaining("<header>");
    try{
      cval=atof(infile->getParam("inflation").c_str());
      burnin=atoi(infile->getParam("burnin").c_str());
      mcmclength= atoi(infile->getParam("mcmclength").c_str());
      cvalstr= atoi(infile->getParam("skip").c_str());
      try{
	datafilestr= infile->getParam("datafilename");
      }catch(string x){cerr<<"WARNING: No data file name stored in the output file."<<endl;}
    }catch(string x){
      cerr<<"Headers missing from output file."<<endl;
      delete(infile);
      return(false);
    }
    delete(infile);
    return(true);
}

int compareDataFiles(string f1, string f2) {
  	if(f1.size()==0 || f2.size()==0){
//	   cerr<<"WARNING! Cannot confirm data file is the same as the MCMC was run on!"<<endl;
	  return(1);
	}else if(false){
// see if they have the same file name but different directories	  
	  return(2);
	}else if(f1!=f2){
//	  cerr<<"WARNING!  You are trying to build a tree from a differently named datafile than the one used for the MCMC!  This might be due to running it on a different system or might imply that the file is incorrect. This may result in strange behaviour!"<<endl;
	  return(-1);
	}
    return(0);
}

}
