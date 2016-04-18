#include "safegetline.h"

std::istream& safeGetline(std::istream& is, std::string& t)
{
  if ( getline( is, t ) ) {
    if ( t.size() && t[t.size()-1] == '\r' ) {
      t = t.substr( 0, t.size() - 1 );
    }
    //t.erase (std::remove (t.begin(), t.end(), '\r'), t.end());
  }
    /*    if(!eofok) if(is.eof()) {
	throw(string("Error: End of file encountered unexpectedly."));
	}*/
    return is;
}
