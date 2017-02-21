#pragma once

#include <vector>
#include <memory>

//Root
#include "Rtypes.h"

//Rootana

//Prevent unused variable warning from this file
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter" 
#include "TRootanaEventLoop.hxx"
#pragma GCC diagnostic pop


class TChain;
class TDataContainer;

namespace myFuncs {
namespace testbeam{
	

//-----------------------------------------------------------
//
//Must setup the files using setupFiles() member function before using the class.
//
//-----------------------------------------------------------
class EventLoopBase : public TRootanaEventLoop {
public:
	
	EventLoopBase(const std::string midasFilesPath, const std::string TOFfilesPath, const int runNum);
	
	void setRunNum(const int runNum) {m_runNum = runNum;}
	inline int getRunNum() const {return m_runNum;}
	
	void setSkipMissingTimeEvents(const bool skipMissingTimeEvents) {m_skipMissingTimeEvents = skipMissingTimeEvents;}
	bool skipMissingTimeEvents() const {return m_skipMissingTimeEvents;}
			
	//Proccess only certain events
	//Actions:
	//Return false if GetEventId() != 1 - process only events with data.
	//Advance timing chain until TOFeventNumber >= midasEventNum.
	//Returns true (proccess event) if midas event serial number == timing chain serial number or if skipMissingTimeEvents() == false
	//Returns else otherwise
	//Dangerous - no bounds checking that m_timingChain does not overflow!!!
	bool PreFilter(TDataContainer& dataContainer) override final;

	
	//Once before running. 
	//Reads first entry in TOF chain.
	void Initialize(void) override final;
	
	void run();
	
private:
	int m_runNum;
	
	//Current entry in m_timingChain chain.
	Long64_t m_timingEntry;
	
	//Max entries in the m_timingChain chain.
	Long64_t m_maxEntries;
	
	//Midas files related to m_runNum
	std::vector<std::string> m_midasFilenames;
	
	//Timing related members
	std::shared_ptr<TChain> m_timingChain;
	
	std::shared_ptr<Long64_t> m_TOFeventNumber;
	std::shared_ptr<double> m_beta;
	std::shared_ptr<double> m_betaError;
	std::shared_ptr<double> m_timeAtCrystal_ns;
	std::shared_ptr<double> m_timeAtCrystalError_ns;
	
	//Use or ignore crystal events that don't have a corresponding TOF time.
	bool m_skipMissingTimeEvents; 
	

	
private:
	
	//Sets ups the TOF chain.
	//Looks for all .root files in TOFfilesPath containg getRunNum() in their filename.
	//There should be only one file
	void setupTOFchain(const std::string TOFfilesPath);
};


}//testbeam namespace
}//myFuncs namespace