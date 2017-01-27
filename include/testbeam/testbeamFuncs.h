#pragma once

//Boost
#include "boost/format.hpp"

//Mine
#include "testbeam/RunDB.h"



namespace myFuncs {
namespace testbeam{
	
	
constexpr double electronMass = 0.5109989461; //MeV/c^2
constexpr double muonMass = 105.6583745; //MeV/c^2
constexpr double chargedPionMass = 139.57018; //MeV/c^2

constexpr double lightSpeed = 0.299792458;  //meters / ns

constexpr double upstreamTOFtoS0Distance = 0.437; //in meters #(average of 43.1 cm, 43.43 cm, 44.56 cm)
constexpr double downstreamTOFtoUpstreamTOF = 3.092; //in meters
constexpr double DownstreamtoS0TOF = downstreamTOFtoUpstreamTOF - upstreamTOFtoS0Distance;

	
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
std::string getNominalMomentumString(const int runNum) {
	int momentum = static_cast<int>(RunDB::instance()[runNum].getNominalBeamMomentum()); //Convert to int to get rid of decimal place
	return std::to_string(momentum) + " MeV/c";
}

//-----------------------------------------------------------
//Return crystal string
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
//Return source distance string with units
//-----------------------------------------------------------
std::string getSourceDistanceString(const int runNum) {	
	int distance = static_cast<int>(RunDB::instance()[runNum].getSourceDistance()); //Convert to int to get rid of decimal place
	if(distance < 0 ) return "no ^{60}Co";
	return "^{60}Co " + std::to_string(distance) + " cm away";
}

std::string getRunParamsString(const int runNum) {	
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



}//testbeam namespace
}//myFuncs namespace