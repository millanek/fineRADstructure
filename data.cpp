#include "data.h"

namespace fines
{

Data::Data(string filename,int ignore,bool yhead,bool xhead)
{
    datacfactor=-1;
    this->filename=filename;
    numignore=0;
    int lineon=0;
    int indrowon=-1;
    string line;
    ifstream file;
    size_t found;
    file.open(filename.data());//Open file
    if(!file.good()) throw(std::string("Invalid file!"));
    while (1)
    {
	safeGetline(file,line);//Read next line from file
	if(line.substr(0,8).compare(string("#Cfactor"))==0) {
	  addDataCfactor(line);
	  lineon++;
	  continue;
	}else if (line.substr(0,1).compare(string("#"))==0) {
	  lineon++;
	  continue;
	}
        if(lineon<ignore) {lineon++;continue;//ignore the first ignore lines
	}else{
	  string teststr=line.substr(0,9); // this is what we need to know to figure out if we have a header line
	  std::transform(teststr.begin(), teststr.end(),teststr.begin(), ::toupper);
	  // detect headers and footers
	  bool isheader=false,hassname=false;
	  stringstream sstest;
	  int inttest;
	  sstest<<line.substr(0,1);
	  if(!(sstest>>inttest)) hassname=true;
	  if((indrowon==-1 &&yhead) || teststr.compare(string("RECIPIENT"))==0)isheader=true;
	  else if((xhead)) hassname=true;
	  if(isheader) {// read the names
	      int colon=0;
	      indrowon=0;
		  while (1)
		  {
		    found=line.find_first_of(",\t ");
		    string tdval=line.substr(0,found);
		    if(colon>0 || ((!xhead) && (teststr.compare(string("RECIPIENT"))!=0))) {
		      string dval=line.substr(0,found);
		      dval=removequotes(dval);
		      names.push_back(dval.c_str());
		    }
		    if(found==string::npos) line.clear();
		    else line=line.substr(found+1,line.length());
		    if(line.length()==0) break;
		    colon++;
		  }
	    lineon++;continue;
	  };// End reading the names
	  if (file.eof() && line.size()==0)
	      break;//Stop if end of file
	  if (line.size()==0 || line[0]=='#' )
	      continue;//Ignore empty lines or comments lines
	  else
	  {//read the row of data
	      int colon=0;
	      data.push_back(vector<double>());
		  while (1)
		  {//read each element of the matrix
		    found=line.find_first_of(",\t ");
		    string dval=line.substr(0,found);//this contains the element
		    dval=removequotes(dval);// get rid of any quotes around names
		    if(found==string::npos) line.clear();
		    else line=line.substr(found+1,line.length());
		    if(colon==0 && hassname) {
		      if(indrowon>=(int)names.size() || indrowon<0) throw(string("Invalid data. Possibly there are names in the data file without RECIPIENT in the top left... you should modify the data file to ChromoPainter format, but you might be able to get around this with the -X and -Y flags."));
 		      if(names[indrowon].compare(dval.c_str())!=0){ 
			cerr<<"ERROR: row and column name "<<lineon-ignore<<" are incompatible!"<<endl;
			cerr<<"row name="<<dval.c_str()<<", col name=";
			if((int)names.size()>lineon-ignore-(int)yhead) cerr<<names[lineon-ignore-(int)yhead]<<endl;
			else cerr<<"Invalid name (number "<<lineon-ignore-(int)yhead<<") requested!"<<endl;
			std::string tfilename=filename;
			tfilename.append(".testdata");
			cout<<"Writing data to "<<tfilename<<" check it manually!"<<endl;
			cout<<"This is the symptom if you have the wrong line endings."<<endl<<"Try converting them to that of your platform (or UNIX are usually handled correctly everywhere)"<<endl;
			cout<<"Alternatively, you may not be using a square matrix."<<endl;
			output(tfilename);			
			throw(std::string("Data error"));
		      }
		      colon++;
		    }else {// read the data value
		      stringstream sstest;
			  sstest.precision(12);
		      int inttest;
		      sstest<<dval;
		      if(!(sstest>>inttest)) throw(string("Invalid data.  Expected a number but received ").append(dval).append("... do you need to use -X and -Y?)"));		      
		      double dbltest=atof(dval.c_str()); // NOTE: a compiler error on cygwin rules out the "correct" sstest>>dbltest; solution
		      data.back().push_back(dbltest);
		      if(line.length()==0) break;
		      colon++;
		    }
		  }
	      lineon++;
	      indrowon++;
		  if (file.eof()) break;
	      continue;
	  }
	}
    }
    file.close();//Close file
    n=data.size();
    if(!testData()){
      throw(std::string("Data error"));
    }
    
    datarowsums.resize(data.size(),0);
    neff.resize(data.size(),1);// number of individuals is 1 per row at this stage
    ignoresuper=vector<bool>(data.size(),false);
    sumXall=0.0;
    for(unsigned int i=0;i<data.size();i++) {
//	if(names.size()>i) cout<<i<<"\t:"<<names[i]<<endl;
//	if((int)names.size()<i) names.push_back(itoa(i));
	for(unsigned int j=0;j<data[i].size();j++) {
		datarowsums[i]+=data[i][j];
	}
	sumXall+=datarowsums[i];
    }
}

Data::Data(std::vector < std::vector<double> > * inmat,std::vector <std::string> namesin)
{
    numignore=0;
    for(unsigned int c1=0;c1< inmat->size();c1++) {
	data.push_back(vector<double>());
	for(unsigned int c2=0;c2< inmat->at(c1).size();c2++) {
	data.back().push_back(inmat->at(c1)[c2]);
        }
    }
    n=data.size();
    testData();
    for(unsigned int c1=0;c1<namesin.size();c1++) {
	names.push_back(namesin[c1]);
    }// we don't need names for all individuals so its ok if this is a different size to above

    datarowsums.resize(data.size(),0);
    sumXall=0.0;
    for(unsigned int i=0;i<data.size();i++) {
//	if(names.size()>i) cout<<i<<"\t:"<<names[i]<<endl;
//	if((int)names.size()<i) names.push_back(itoa(i));
	for(unsigned int j=0;j<data[i].size();j++) {
		datarowsums[i]+=data[i][j];
	}
	sumXall+=datarowsums[i];
    }
}

Data::Data(Data * datin)
{
    n=datin->n;
    data=datin->data;
    names=datin->names;
    datarowsums=datin->datarowsums;
    neff=datin->neff;
    sumXall=datin->sumXall;
    ignoresuper=datin->ignoresuper;
    numignore=datin->numignore;
}

void Data::print(ostream * out)
{
    (*out)<<"NAME, ";
    for (unsigned int i=0;i<data.size()-1;i++)(*out)<<getnames(i)<<", ";
    (*out)<<getnames(data.size()-1)<<endl;
    for (unsigned int i=0;i<data.size();i++)
    {
		(*out)<<getnames(i)<<", "<<flush;
        for (int j=0;j<(((int)data[i].size())-1);j++)
        {
            (*out)<<data[i][j]<<", "<<flush;
        }
        if(data[i].size()!=00) (*out)<<data[i][data[i].size()-1];
		(*out)<<endl;
    }	
}

void Data::output(std::string filename)
{
  filebuf fb;
  fb.open (filename.c_str(),ios::out);
  ostream os(&fb);
  print(&os);
  fb.close();
/* ofstream file;
  file.open (filename.c_str());
    print(file);
  file.close();*/
}

void Data::applySuper(std::string name,std::vector<int> ind){
// sort into descending order so that removing is safe
  bool ignoresuperb=false;
  if(name[0]=='*') {
    ignoresuperb=true;
    name=name.substr(1);
    numignore++;
  }
  sort (ind.begin(), ind.end());
  if(names.size()<data.size())names.resize(data.size(),string());
  for(unsigned int i=1;i<ind.size();i++){
    for(unsigned int c1=0;c1<data.size();c1++){
      data[ind[0]][c1]+=data[ind[i]][c1];
      data[c1][ind[0]]+=data[c1][ind[i]];
    }
    neff[ind[0]]+=neff[ind[i]];
    datarowsums[ind[0]]+=datarowsums[ind[i]];
    ignoresuper[ind[0]]=ignoresuperb;
  }
//  for(int c1=0;c1<ind.size();c1++) cout<<ind[c1]<<",";
//  cout<<endl;
  for(unsigned int i=ind.size()-1;i>0;i--){
   for(unsigned int c1=0;c1<data.size();c1++) data[c1].erase(data[c1].begin() + ind[i]);
    data.erase(data.begin() + ind[i]);
    neff.erase(neff.begin() + ind[i]);
    names.erase(names.begin() + ind[i]);
    datarowsums.erase(datarowsums.begin() + ind[i]);
    ignoresuper.erase(ignoresuper.begin() + ind[i]);
  }
  names[ind[0]]=name;
}

void Data::makeSuper(std::string superstring){
  forceNames();
  superstring.erase (std::remove (superstring.begin(), superstring.end(), ' '), superstring.end());
  superstring.erase (std::remove (superstring.begin(), superstring.end(), '\t'), superstring.end());
  superstring.erase (std::remove (superstring.begin(), superstring.end(), '"'), superstring.end());
  
  const char* test = ",()";
  int gnum=0;
  vector<int> ind;
  string dname;
  size_t found=superstring.find_first_of(test), pos=0;
	while(found!=string::npos){
	  if(superstring.at(pos)==')'){// should be a name, or the end
	   applySuper(dname,ind);
	   dname=superstring.substr(pos+1,found-pos-1);
	  }else if(pos==0){// should be a name
	    dname=superstring.substr(pos,found-pos);
	  }
	  if(superstring.at(pos)=='(') {// new population
	    if(dname.size()==0){ // there is no name provided for this population
	      stringstream ss;
	      ss<<string("Group")<<gnum; 
	      dname=ss.str();
	    }
	    ind.clear();
	    if(found==0) {
	      cout<<"Invoking found=0"<<endl;
	      found=superstring.find_first_of(test,pos+1);
	    }
	    if (strchr(test, superstring.at(pos+1)) != NULL) {
	      cerr<<"Invalid fixed file structure: two control characters appear together (population: \""<<dname<<"\")."<<endl;
	      throw(std::string("State fixed file error: name not found!"));
	    }
	  }
	  if(superstring.at(pos)=='(' || superstring.at(pos)==',') { // expect to see the next in a list of individuals
	    if (strchr(test, superstring.at(pos+1)) != NULL) {
	      cerr<<"Invalid fixed file structure: two control characters appear together (population: \""<<dname<<"\")."<<endl;
	      throw(std::string("State fixed file error: name not found!"));
	    }
	    int indiv=getIndex(superstring.substr(pos+1,found-pos-1));
	    if(indiv<0) {
	      cerr<<"Individual \""<<superstring.substr(pos+1,found-pos-1).c_str()<<"\" not recognised (in population \""<<dname<<"\")"<<endl;
	      cerr<<"Options are: ";
	      for(unsigned int c1=0;c1<data.size();c1++)cerr<<getnames(c1)<<" ";
	      cerr<<endl;
	      throw(std::string("State fixed file error: name not found!"));
	    }
	    ind.push_back(indiv);
	  }
	  pos=found;
	  found=superstring.find_first_of(test,pos+1);
	}
	if(ind.size()>0) applySuper(dname,ind);
}

void Data::makeSuperFromFile(std::string filename){
      string line, superstring;
      ifstream file;
      file.open(filename.data());//Open file
      if(!file.good()) throw(std::string("Invalid super individual file!"));
      while(1){
	safeGetline(file,line);
	if(line.size()>8) if(line.substr(0,8).compare(string("#Cfactor"))==0) addDataCfactor(line);
	if(line.size()>0) if(line[0]!='#')superstring.append(line); // Ignore comment lines
	if (file.eof())
	      break;//Stop if end of file
      }
      file.close();
      makeSuper(superstring);
}
void Data::addDataCfactor(string line){
  size_t found=line.find_first_of(",\t ");
  string cstring=line.substr(found);
  datacfactor=atof(cstring.c_str());
}

Data::~Data()
{}

} // end namespace weakarg
