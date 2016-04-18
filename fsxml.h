#ifndef FSXML_H
#define FSXML_H
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <ios>

#include "safegetline.h" // for cross-platform safely getting lines

namespace fines
{

/**
    @brief Reads in the output file
*/
class FsXml
{
public:
	FsXml(std::string f);
	~FsXml();
	std::streampos getLineContaining(std::string str,bool getlast=false);
	std::streampos gotoLineContaining(std::string str,bool getlast=false);
	inline std::streampos getLastIteration(){
		return(getLineContaining("<Iteration>",true));
	};///<Gets the last iteration
	std::string getLine();
	inline std::streampos gotoNextLineContaining(std::string str,bool getlast=false){
		getLine();
		return(gotoLineContaining(str,getlast));
	}
	inline void seekg(std::streampos is){iterfile.seekg(is);};
	inline std::streampos tellg(){return(iterfile.tellg());};
	std::string getTree(std::streampos itstart=-1);
	std::string getParam(std::string tag, std::streampos itstart=-1);
	inline bool eof(){return(iterfile.eof());};
	inline void clear(){iterfile.clear();};
	void rewind(){
	  clear();
	  std::streampos pos=-1;
	  iterfile.seekg(pos);
	}
	inline std::string getLine(bool *ok){
		std::string res=getLine();
		if ((res.size() > 0) && (res[res.size() - 1] == '\r')) res.resize(res.size() - 1);
		if(ok!=NULL)if(iterfile.eofbit) *ok=false;
		return(res);
	}
	void rangeextract(long from, long to,std::ostream * out);
	void thin(long by,std::ostream * out);
	inline std::string getName(){return(fname);}
private:
	std::string fname;
	std::ifstream iterfile;
};

}//end namespace
#endif
