#include <string>


namespace myFuncs
{
  //--------------------------------------------------------------------------------------------
  //getRunNum
  //********************************************************************************************
  //Get a filename (may include directories) and return the run number string.
	//This is without the partial run number.
  //--------------------------------------------------------------------------------------------
  inline std::string getRunNum(const std::string& filename)
	{
		std::string trimLeft = filename.substr(filename.find("00000", 0) + 5);
		return trimLeft.substr(0,trimLeft.find("_000"));	
	}
  
}