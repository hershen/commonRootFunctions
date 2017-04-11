#pragma once

//STL
#include <memory> //for unique_ptr

//Boost
#include "boost/format.hpp"

class TGraphErrors;

namespace myFuncs {
	class PaveText;
namespace testbeam{

constexpr std::array<int, 20> representitiveRuns =  {591, 
																								 597, 
																								 599, 
																								 600, 
																								 589,
																								 
																								 627, //Chinese
																								 611,
																								 625,
																								 653,
																								 654,

																								 
																								 662, //Babar
																								 664,
																								 666,
																								 676,
																								 687,
																								 688,
																								 
// 																								 713, //Ukranian
																								 715,
// 																								 716,
																								 717,
																								 730,
																								 731};

//-----------------------------------------------------------
//Return Nominal beam momentum string with units
//-----------------------------------------------------------
std::string getNominalMomentumString(const int runNum);

//-----------------------------------------------------------
//Return crystal string
//-----------------------------------------------------------
std::string getCrystalString(const int runNum);

//-----------------------------------------------------------
//Return source distance string with units
//-----------------------------------------------------------
std::string getSourceDistanceString(const int runNum);

inline std::string getRunParamsString(const int runNum) {	
	return boost::str(boost::format("Run %1%, %2%, %3%, %4%")% runNum % getCrystalString(runNum) % getNominalMomentumString(runNum) % getSourceDistanceString(runNum)); 
}

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

//-----------------------------------------------------------
//Get all root files in pathToFiles with runNum in their name and ending in extension.
//-----------------------------------------------------------
std::vector<std::string> getFilesRelatedToRun(const std::string pathToFiles, const int runNum, const std::string = ".root");

inline bool isCrystalChannel(const int channel) {
	return (channel == 1 or channel == 15);
}

//Return TGraphErrors with voltage as function of time
std::unique_ptr<TGraphErrors> getWaveformGraph(const std::vector<double>& voltage, const std::vector<double>& errors = std::vector<double>()); 


//Get PaveText of run parameters, channel and event num
std::unique_ptr<myFuncs::PaveText> getRunChannelEventPaveText(const int runNum, const int channelNum, const int eventNum);

void drawWaveform(const std::vector<double>& voltages, const int runNum, const int channelNum, const size_t eventNum, const bool waitPrimitive = true);

}//testbeam namespace
}//myFuncs namespace
