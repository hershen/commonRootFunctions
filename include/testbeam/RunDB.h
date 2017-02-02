#pragma once

#include <vector>
#include <unordered_map>

//Mine
#include "testbeam/constants.h"

namespace myFuncs {
namespace testbeam{
	
enum class Crystal {CsI_Tl_Belle, CsI_Tl_Babar, CsI_Ukrainian, CsI_Chinese};


//-----------------------------------------------------------
//Class representing one TB run and all associated parameters.
//-----------------------------------------------------------
class RunParams
{
public:
	constexpr RunParams(const int runNum, const Crystal crystal, const double sourceDistance, const double HV, const double nominalBeamMomentum, const double crystalFrontFaceToIncubatorSideWallDistance):
	m_runNum(runNum),
	m_crystal(crystal),
	m_sourceDistance(sourceDistance), 
	m_HV(HV),
	m_nominalBeamMomentum(nominalBeamMomentum),
	m_crystalFrontFaceToIncubatorSideWallDistance(crystalFrontFaceToIncubatorSideWallDistance)
	{}
	
	int getRunNum() const {return m_runNum;}
	
	Crystal getCrystal() const {return m_crystal;}
	
	// Negative means no source
	double getSourceDistance() const {return m_sourceDistance;}
	
	double getHV() const {return m_HV;}
	
	// Returns the nominal beam momentum - i.e. the momentum we thought the runs were at.
	//[MeV/c]. Negative means negative charged particles in beam
	double getNominalBeamMomentum() const {return m_nominalBeamMomentum;}
	
	//Returns the downstream TOF center to crystal center distance in meters.
	/*constexpr */double getDownstream2crystalCenterDistance() const {
		return c_downstreamCenter2incubatorWall + 
					 c_incubatorWallSideWidth + 
					 m_crystalFrontFaceToIncubatorSideWallDistance + //different per run
					 0.5 * c_crystalLength;
	}
	
private:
	int m_runNum;
	Crystal m_crystal;
	double m_sourceDistance; // Negative means no source
	double m_HV; // [V]
	double m_nominalBeamMomentum; // [MeV/c]. Negative means negative charged particles in beam
	double m_crystalFrontFaceToIncubatorSideWallDistance; // meters
};


//-----------------------------------------------------------
//Singleton class to access DB
//-----------------------------------------------------------
class RunDB
{
public:
	
	//-----------------------------------------------------------
	//Get instance to singleton
	//-----------------------------------------------------------
	static RunDB& instance()
	{
			static RunDB instance; // Guaranteed to be destroyed.
														   // Instantiated on first use.
			return instance;
	}
private:
	RunDB();

public:
	RunDB(RunDB const&)         = delete;
	void operator=(RunDB const&)  = delete;

	//access RunParams with run number runNum
	const RunParams& operator[](const int runNum) const;
	
private:
//This holds the DB	
const std::unordered_map<int, const RunParams> m_DB;


	
};

}//testbeam namespace
}//myFuncs namespace