#include "infadmixture.h"

using namespace std;
namespace fines
{

InfAdmixture::InfAdmixture(Data *d,State *s,double alpha,bool v,double initialq,bool test)
{
	atest=test;
	data=d;
	state=s;
	verbose=v;
	resetCounters();
	double numvars=100;
	P_QS=100.0/numvars;
	//P_PS =1-sum(these)
	for(int c1=0;c1<data->getDim();c1++) {
		Q.push_back(vector<double>(state->getP(),0.0));
	}
	QcolSums=vector<double>(state->getP(),0.0);
	for(int c1=0;c1<state->getP();c1++) {
		P.push_back(vector<double>(state->getP(),0.0));
	}
	samplePs(true);
	initialQs(initialq);
	setAlpha(alpha);
	logstateprob=state->admixtureLogLikelihood(&Q,&P,&QcolSums);
}

void InfAdmixture::initialQs(double spread)
{
	for(unsigned int c1=0;c1<Q.size();c1++) {
		Q[c1][state->getPop(c1)]=1.0-spread;
		QcolSums[state->getPop(c1)]+=Q[c1][state->getPop(c1)];
		for(unsigned int c2=0;c2<Q[c1].size();c2++) if((int)c2!=state->getPop(c1)){
			Q[c1][c2]=spread/(Q[c1].size()-1.0);
			QcolSums[c2]+=Q[c1][c2];
		}
	}
}

void InfAdmixture::samplePs(bool mean)
{
	if(mean){
	  for(int c1=0;c1<state->getP();c1++) {
	    for(int c2=0;c2<state->getP();c2++) {
		P[c1][c2]=(state->sumXab(c1,c2)) + state->getBeta(c1,c2);
	    }
	}}else{
	  std::vector <double> beta(state->getP(),1.0);
	
	  for(int c1=0;c1<state->getP();c1++) {
		for(int c2=0;c2<state->getP();c2++) beta[c2]= state->getBeta(c1,c2);
		RDirichlet(&beta,&(P[c1]));
	  }
	}
	for(int c1=0;c1<state->getP();c1++) {
		double rs=0;
		for(int c2=0;c2<state->getP();c2++) rs+=P[c1][c2];
		for(int c2=0;c2<state->getP();c2++) P[c1][c2]/=rs;
	}
}

void InfAdmixture::printPs(std::ostream * out)
{
	*out<<"P,";
	for(unsigned int i=0;i<P[0].size()-1;i++) *out<<"Pop"<<i<<", ";
	*out<<"Pop"<<P[0].size()-1<<endl;
	for(unsigned int i=0;i<P.size();i++) {
	   *out<<"Pop"<<i<<", ";
	  for(unsigned int j=0;j<P[i].size()-1;j++) *out<<P[i][j]<<",";
	  *out<<P[i][P[i].size()-1]<<endl;
	}
}

void InfAdmixture::printQs(std::ostream * out)
{
	*out<<"Q,";
	for(unsigned int i=0;i<Q[0].size()-1;i++) *out<<"Pop"<<i<<", ";
	*out<<"Pop"<<Q[0].size()-1<<endl;
	for(unsigned int i=0;i<Q.size();i++) {
	  *out<<data->getnames(i)<<", ";
	  for(unsigned int j=0;j<Q[i].size()-1;j++) *out<<Q[i][j]<<",";
	  *out<<Q[i][Q[i].size()-1]<<endl;
	}
}

int InfAdmixture::moveProw(int a)
{
	int r1=RandomInteger(0,P[a].size()-1),r2=r1;
	while(r2==r1) r2=RandomInteger(0,P[a].size()-1);
	double x=rnd();// the proportion of r1 to move to r2
	double oldloglik=logstateprob;
	double oldpriorP=priorForP(a);
	double op1=P[a][r1],op2=P[a][r2];
	P[a][r1]=(op1+op2) * x;
	P[a][r2]=(op1+op2) * (1.0-x);
	logstateprob=state->admixtureLogLikelihood(&Q,&P,&QcolSums);
	if(log(rnd())>logstateprob + priorForP(a) - oldloglik - oldpriorP) {
		//rejected
		P[a][r1]=op1;	
		P[a][r2]=op2;
		logstateprob=oldloglik;
		return(0);
	}
	// else accept
	return(1);
}

/*int InfAdmixture::moveQrow(int a)
{
	int r1=RandomInteger(0,Q[a].size()-1),r2=r1;
	while(r2==r1) r2=RandomInteger(0,Q[a].size()-1);
	double x=rnd();// the proportion of r1+r2 to move to r2
	double oldloglik=state->admixtureLogLikelihoodIndiv(a,&Q,&P,&QcolSums);
//	double oldloglik=state->admixtureLogLikelihood(&Q,&P,&QcolSums);
	double oldpriorQ=priorForQ(a);
	double oq1=Q[a][r1],oq2=Q[a][r2];

	double sr1=0;
	for(int c1=0;c1<Q.size();c1++) sr1+=Q[c1][r1];
	if(fabs(sr1-QcolSums[r1])>0.0005) throw(string("Colsums error!"));
	
	Q[a][r1]=(oq1+oq2) * x;
	Q[a][r2]=(oq1+oq2) * (1.0-x);
	QcolSums[r1]+= Q[a][r1] - oq1;
	QcolSums[r2]+= Q[a][r2] - oq2;

	//double newloglik=state->admixtureLogLikelihood(&Q,&P,&QcolSums);
	double newloglik=state->admixtureLogLikelihoodIndiv(a,&Q,&P,&QcolSums);
//	cout<<"p="<<logstateprob<<" + "<<priorForQ(a)<<" - "<<oldloglik<<" - "<<oldpriorQ<<endl;
	if(log(rnd())>newloglik + priorForQ(a) - oldloglik - oldpriorQ) {
		//rejected
//		cout<<"rejected"<<endl;
		QcolSums[r1]-= Q[a][r1] - oq1;
		QcolSums[r2]-= Q[a][r2] - oq2;
		Q[a][r1]=oq1;
		Q[a][r2]=oq2;
		return(0);
	}
	// else accept
//	cout<<"accepted!"<<endl;
	logstateprob += newloglik-oldloglik;// update our running total of the likelihood
	return(1);
}*/

int InfAdmixture::moveQrow(int a)
{
	int r1=RandomInteger(0,Q[a].size()-1),r2=r1;
	while(r2==r1) r2=RandomInteger(0,Q[a].size()-1);
	double x=rnd();// the proportion of r1+r2 to move to r2
//	double oldloglik=state->admixtureApproxSetLogLikelihood(r1,&Q,&P)+state->admixtureApproxSetLogLikelihood(r2,&Q,&P);
	double oldloglik=state->admixtureApproxIndLogLikelihood(a,&Q,&P);
//	double oldloglik=state->admixtureLogLikelihood(&Q,&P,&QcolSums);
	double oldpriorQ=priorForQ(a);
	double oq1=Q[a][r1],oq2=Q[a][r2];

	if(atest){
	Q[a][r1]=oq2;
	Q[a][r2]=oq1;
	}else{
	Q[a][r1]=(oq1+oq2) * x;
	Q[a][r2]=(oq1+oq2) * (1.0-x);
	}
	//double newloglik=state->admixtureLogLikelihood(&Q,&P,&QcolSums);
	//double newloglik=state->admixtureApproxIndLogLikelihood(a,&Q,&P);
	double newloglik=state->admixtureApproxIndLogLikelihood(a,&Q,&P);
//	double newloglik=state->admixtureApproxSetLogLikelihood(r1,&Q,&P)+state->admixtureApproxSetLogLikelihood(r2,&Q,&P);
//	cout<<"p="<<logstateprob<<" + "<<priorForQ(a)<<" - "<<oldloglik<<" - "<<oldpriorQ<<endl;
	if(log(rnd())>newloglik + priorForQ(a) - oldloglik - oldpriorQ) {
		//rejected
//		cout<<"rejected"<<endl;
		Q[a][r1]=oq1;
		Q[a][r2]=oq2;
		return(0);
	}
	// else accept
//	cout<<"accepted!"<<endl;
	logstateprob += newloglik-oldloglik;// update our running total of the likelihood
	return(1);
}

/*int InfAdmixture::moveQrow(int a)
{
	if(verbose)cout<<"Proposing to move Qs for individual "<<a<<":"<<flush;
	int opop=state->getPop(a);
	int r1=opop;
	while(r1==opop) r1=RandomInteger(0,Q[a].size()-1);	
	double oq1=Q[a][opop],oq2=Q[a][r1];
//	double oldPost=state->posteriorSetProb(opop)+state->posteriorSetProb(r1);
	double oldloglik=state->admixtureApproxIndLogLikelihood(a,&Q,&P);
//	double oldloglik=state->admixtureApproxSetLogLikelihood(opop,&Q,&P)+state->admixtureApproxSetLogLikelihood(r1,&Q,&P);
//	double oldloglik=state->admixtureLogLikelihoodIndiv(a,&Q,&P,&QcolSums);
//	double oldloglik=state->admixtureLogLikelihood(&Q,&P,&QcolSums);
	double oldPost=oldloglik+priorForQ(a);
	QcolSums[opop]+= oq2 - oq1;
	QcolSums[r1]+= oq1 - oq2;
	Q[a][opop]=oq2;
	Q[a][r1]=oq1;

//	state->moveInd(a,r1);
//	samplePs(true);

//	double newPost=state->posteriorSetProb(opop)+state->posteriorSetProb(r1);
	double newloglik=state->admixtureApproxIndLogLikelihood(a,&Q,&P);
//	double newloglik=state->admixtureApproxSetLogLikelihood(opop,&Q,&P)+state->admixtureApproxSetLogLikelihood(r1,&Q,&P);
//	double newloglik=state->admixtureLogLikelihoodIndiv(a,&Q,&P,&QcolSums);
//	double newloglik=state->admixtureLogLikelihood(&Q,&P,&QcolSums);
//cout<<"oll="<<oldloglik<<" nll="<<newloglik<<endl;
	double newPost=newloglik+priorForQ(a);
//	cout<<"postdiff="<<(newPost-oldPost)<<endl;
	if(log(rnd())>newPost-oldPost) {
		//rejected
		if(verbose)cout<<"rejected"<<endl;
		QcolSums[opop]-= oq2 - oq1;
		QcolSums[r1]-= oq1 - oq2;
		Q[a][opop]=oq1;
		Q[a][r1]=oq2;
//		state->moveInd(a,opop);
//		samplePs(true);
		return(0);
	}
	// else accept
	if(verbose)cout<<"accepted!"<<endl;
	logstateprob += newloglik-oldloglik;// update our running total of the likelihood
	return(1);
}*/

int InfAdmixture::moveQs()
{
	int a=0;
	for(unsigned int c1=0;c1<Q.size();c1++) a+=moveQrow(c1);
	return(a);
}

int InfAdmixture::movePs()
{
	int a=0;
	for(unsigned int c1=0;c1<P.size();c1++) a+=moveProw(c1);
	return(a);
}

void InfAdmixture::merge(int a, int b){
// merge population b into population a
	int keep=min(a,b);
	int lose=max(a,b);
	// update Q
	for(unsigned int c1=0;c1<Q.size();c1++){
		Q[c1][keep]+=Q[c1][lose];
		Q[c1].erase(Q[c1].begin()+lose);
	}
	// and its column sums
	QcolSums[keep]+=QcolSums[lose];
	QcolSums.erase(QcolSums.begin()+lose);
	// update the prior vector alpha
	alphavec.erase(alphavec.begin()+lose);
	// update P
	for(unsigned int c1=0;c1<P.size();c1++) {
		if((int)c1!=lose){
		P[c1][keep]+=P[c1][lose];
		P[keep][c1]+=P[lose][c1];
		}else P[keep][keep]+=P[lose][lose];
	}
	for(unsigned int c1=0;c1<P.size();c1++) {
		P[c1].erase(P[c1].begin()+lose);
	}
	P.erase(P.begin()+lose);
	// update the state, too
	state->merge(keep,lose);
	test();
}



double InfAdmixture::priorForQ(int a)
{
	return(state->LDirichletProb(getAlphaVec(),Q[a]));
}

double InfAdmixture::priorForP(int a)
{
	return(state->LDirichletProb(state->getBetaVector(a),P[a]));
}

void InfAdmixture::metropolis(int prevints,int numints, int thin,std::ostream * fileout)
{
	for(long c1=0;c1<numints;c1++) {
		if (numints>50 && (c1)%((numints)/50)==0)
		{
			if (100l*c1/(numints)<=1){
				cout<<"#  "<<100l*(double)c1/(numints)<<"%"<<flush;
			}else if (100l*c1/(numints)<10)
				cout<<"\b\b\b\b#  "<<(int)(100l*(double)c1/(numints))<<"%"<<flush;
			else
				cout<<"\b\b\b\b# "<<(int)(100l*(double)c1/(numints))<<"%"<<flush;
		}
		if (c1+1==numints)
			cout<<"\b\b\b\b# 100%"<<endl<<flush;
		double x=rnd();
		if(x<P_QS) {
			numQs+=Q.size();
			accQs+=moveQs();
		}else {
			numPs+=P.size();
			accPs+=movePs();
		}
		if(c1 % thin==0 && fileout!=NULL) {exportXmlIter(fileout,c1+prevints);}
	}
}

void InfAdmixture::test()
{
	for(int i=0;i<(int) Q[0].size();i++) {
		double qcolsum=0;
		for(int j=0;j<(int) Q.size();j++) {
			qcolsum+=Q[j][i];
		}
		if(fabs(qcolsum-QcolSums[i])>0.0005) {
			cout<<"colsum Q["<<i<<"]="<<qcolsum<<" but we have "<<QcolSums[i]<<" stored!"<<endl;
			throw(string("Colsums for Q dont match"));
		}
	}
}

void InfAdmixture::exportXmlHead(std::ostream * fileout)
{
	*fileout<<"<?xml version = '1.0' encoding = 'UTF-8'?>"<<endl<<"<outputFile type='admixture'>"<<endl;
}

void InfAdmixture::exportXmlTail(std::ostream * fileout)
{
	*fileout<<"</outputFile>"<<endl;
}

void InfAdmixture::exportXmlIter(std::ostream * fileout,int iter)
{
//	state->printBeta(&cout);
	*fileout<<"<Iteration>"<<endl;
	state->setprint(fileout);
	*fileout<<"<Likelihood>"<<setprecision (10) <<logstateprob<<"</Likelihood>"<<endl;
	double lq=0;
	for(int c1=0;c1<(int)Q.size();c1++) lq+=priorForQ(c1);
	*fileout<<"<PriorForQ>"<<lq<<"</PriorForQ>"<<endl;
	*fileout<<"<Q>"<<endl;
	printQs(fileout);
	*fileout<<"</Q>"<<endl;
	if(P_QS<1.0) {
		*fileout<<"<P>"<<endl;
		printPs(fileout);
		*fileout<<"</P>"<<endl;
	}
	double np=0,nq=0;
	if (numPs>0) np=((double)accPs)/numPs;
	if (numQs>0) nq=((double)accQs)/numQs;
	*fileout<<"<AccPs>"<<np<<"</AccPs>"<<endl;
	*fileout<<"<AccQs>"<<nq<<"</AccQs>"<<endl;
	*fileout<<"<Number>"<<iter<<"</Number>"<<endl;
	*fileout<<"</Iteration>"<<endl;
}


InfAdmixture::~InfAdmixture()
{
}

} // end namespace fines
