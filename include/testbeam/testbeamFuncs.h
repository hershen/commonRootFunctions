#pragma once

#ifndef TESTBEAMFUNCS_H
#define TESTBEAMFUNCS_H

//Boost
#include "boost/format.hpp"


namespace myFuncs {
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
//Get all root files in pathToFiles with runNum in their name and ending in '.root'.
//-----------------------------------------------------------
std::vector<std::string> getFilesRelatedToRun(const std::string pathToFiles, const int runNum);


}//testbeam namespace
}//myFuncs namespace

#endif