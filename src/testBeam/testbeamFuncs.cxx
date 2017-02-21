#include "testbeam/testbeamFuncs.h"

//Mine
#include "fileFuncs.h"
#include "testbeam/RunDB.h"

namespace myFuncs {
namespace testbeam{

//-----------------------------------------------------------
//getNominalMomentumString
//-----------------------------------------------------------
std::string getNominalMomentumString(const int runNum) {
	int momentum = static_cast<int>(RunDB::instance()[runNum].getNominalBeamMomentum()); //Convert to int to get rid of decimal place
	return std::to_string(momentum) + " MeV/c";
}

//-----------------------------------------------------------
//getCrystalString
//-----------------------------------------------------------
std::string getCrystalString(const int runNum) {
	
	Crystal crystal = RunDB::instance()[runNum].getCrystal();
	if(crystal == Crystal::CsI_Tl_Belle) return "CsI(Tl) Belle";
	else if(crystal == Crystal::CsI_Tl_Babar) return "CsI(Tl) Babar";
	else if(crystal == Crystal::CsI_Ukrainian) return "CsI AMCRYS";
	else if(crystal == Crystal::CsI_Chinese) return "CsI SICCAS";
	else return "Unknown crystal";
}

//-----------------------------------------------------------
//getSourceDistanceString
//-----------------------------------------------------------
std::string getSourceDistanceString(const int runNum) {	
	int distance = static_cast<int>(RunDB::instance()[runNum].getSourceDistance()); //Convert to int to get rid of decimal place
	if(distance < 0 ) return "no ^{60}Co";
	return "^{60}Co " + std::to_string(distance) + " cm away";
}

std::vector<std::string> getFilesRelatedToRun(const std::string pathToFiles, const int runNum, const std::string extension)
{
	const std::vector<std::string> allFilenames = myFuncs::getFilesEndingWith( pathToFiles, extension );
	std::vector<std::string> filenames;
	//Find relevant files
	const std::string runString = std::to_string( runNum );
	for( const auto filename : allFilenames) if (filename.find(runString) != std::string::npos) filenames.push_back(filename);
	
	return filenames;
}

}//testbeam namespace
}//myFuncs namespace
