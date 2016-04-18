#include "fsxml.h"
#include <fstream>
#include <sstream>
#include <cstdlib>

using namespace std;
namespace fines
{

FsXml::FsXml(std::string f)
{
	fname=f;
	iterfile.open(fname.data());
        if (!iterfile) {
            cerr << "Can't open input file " << f << endl;
            throw(std::string("Unable to open file"));
        }
}

FsXml::~FsXml()
{
	if(iterfile.is_open()) iterfile.close();
}

streampos FsXml::gotoLineContaining(std::string str,bool getlast){
	streampos pos =iterfile.tellg (),it=-1;
	string res;
	size_t found;
	int len=str.length();
	while(1){
		pos =iterfile.tellg ();
		res=getLine();
		found=res.find_first_of('<');
		if(found!=string::npos) if(res.length()>=len+found) if(res.substr(found,len).compare(str)==0) {it=pos;}
		if (iterfile.eof()||!iterfile.good()||(!getlast && it>=0)) break;
	}
	iterfile.clear();
	iterfile.seekg(it);
	return(it);
}

std::streampos FsXml::getLineContaining(std::string str,bool getlast){
	std::streampos p0=iterfile.tellg();
	std::streampos r=gotoLineContaining(str,getlast);
	iterfile.seekg(p0);
	return(r);
}

std::string FsXml::getParam(string tag, streampos itstart)
{
	string res;
	string stag=tag,etag=tag;
	stag.insert(0,"<");
	etag.insert(0,"</");
	stag.append(">");
	etag.append(">");
	if(itstart>=0) iterfile.seekg(itstart);
	else itstart=iterfile.tellg ();
	streampos p=gotoLineContaining(stag.c_str(),false);
	if(p<0){
		cout<<"Tag "<<stag<<" not found in "<<fname<<endl;
		iterfile.clear();
		iterfile.seekg(itstart);		
		throw(string("Tag not found in file!"));
	}
    res=getLine();
	size_t found=res.find(stag.c_str());
	res.erase(found,stag.length());
	size_t found2=res.find(etag.c_str());
	res=res.substr(found,found2);
	iterfile.seekg(itstart);
	return(res);
}

std::string FsXml::getTree(streampos itstart)
{
	string res;
	size_t found;
	if(itstart>=0) iterfile.seekg(itstart);
	else itstart=iterfile.tellg ();
	getLineContaining("<Tree>",true);
	while((found=res.find_first_of('('))==string::npos) {
        res=getLine();
		if(iterfile.eof())throw(string("Error in readtree: <Tree> does not appear to contain a newick tree!"));
	}
	size_t found2=res.find("</Tree>");
	if(found2!=string::npos) found2-=6;
	res=res.substr(found,found2);
	iterfile.seekg(itstart);
	return(res);
}

void FsXml::rangeextract(long from, long to,std::ostream * out)
{
	clear();
	string res;
	size_t found;
	long iteron=0;
	if(from<0) from=0;
	if(to<0) to=0;
	bool putout=true;
	while(1){
		res=getLine();
		if((found=res.find("<Iteration>"))!=string::npos && (iteron<from || iteron>to))putout=false;
		if((found=res.find("<Iteration>"))!=string::npos && (iteron>=from && iteron<=to))putout=true;
		found=res.find("</Iteration>");
		size_t found2=res.find("</outputFile>");
		if(found2!=string::npos) putout=true;
		if(putout) *out<<res<<endl;
		if(found!=string::npos) iteron++;
		if (iterfile.eof()||!iterfile.good()) break;
	}
	iterfile.clear();
}

void FsXml::thin(long by,std::ostream * out)
{
	clear();
	string res;
	long iteron=0;
	if(by<=0) by=1;
	while(1){
		res=getLine();
		size_t found=res.find("</Iteration>");
		size_t found2=res.find("</outputFile>");
		if((iteron%by==0)||found2!=string::npos) *out<<res<<endl;
		if(found!=string::npos) iteron++;
		if (iterfile.eof()||!iterfile.good()) break;
	}
	iterfile.clear();

}

std::string FsXml::getLine(){ 
    std::string res;
    safeGetline(iterfile, res);
    return(res);
}

  /*
std::string FsXml::getLine(){
streampos tpos =iterfile.tellg ();
    std::string res;
    char c;
	while(1){
	    iterfile.get ( c );
        if (iterfile.eof()||!iterfile.good()) break;
	    if((c=='\n'||c=='\r') &&res.size()==0) continue;
        if(c=='\n'||c=='\r') break;
        res.push_back(c);
	}
	iterfile.seekg(tpos+(streampos)(res.size()+1));
	return(res);
}
*/

}//end namespace

