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
class MidasLoop : public TRootanaEventLoop {
public:

	MidasLoop(const std::string midasFilesPath, const int runNum);

	void setRunNum(const int runNum) {m_runNum = runNum;}
	inline int getRunNum() const {return m_runNum;}

	inline void setSkipNoiseEvents(const bool skipNoiseEvents) {m_skipNoiseEvents = skipNoiseEvents;}
	inline bool isSkipNoiseEvents() const {return m_skipNoiseEvents;}

	//Used in case only noise events are necessary
	inline void setSkipSignalEvents(const bool skipSignalEvents) {m_skipSignalEvents = skipSignalEvents;}
	inline bool isSkipSignalEvents() const {return m_skipSignalEvents;}

	//Proccess only certain events
	//Actions:
	//Return false if GetEventId() != 1 - process only events with data.
	//Advance timing chain until TOFeventNumber >= midasEventNum.
	//Returns true (proccess event) if midas event serial number == timing chain serial number or if skipNoiseEvents() == false
	//Returns else otherwise
	//Dangerous - no bounds checking that m_timingChain does not overflow!!!
	bool PreFilter(TDataContainer& dataContainer) override final;

	void run(const std::string options = "");

	//Set maximum files to process.
	inline void setMaxFiles(const unsigned int maxFiles) {m_maxNumberOfFiles = maxFiles; }

	//Get Midas event number
	inline uint32_t getEventNum(const TDataContainer& dataContainer) const {return dataContainer.GetMidasEvent().GetSerialNumber();}

	//Files to skip before processing
	inline int getNumFilesToSkip() const {return m_filesToSkip;}
	inline void setNumFilesToSkip(const int filesToSkip) {m_filesToSkip = filesToSkip;}

private:
	//Midas files related to m_runNum
	std::vector<std::string> m_midasFilenames;

	//Maximum amplitude used in prefilter
	double m_maxAmplitude;

//Max entries in the m_timingChain chain.
	Long64_t m_maxEntries;

	int m_runNum;

	//Number of files to skip when running
	int m_filesToSkip;

	unsigned int m_maxNumberOfFiles;



	//Use or ignore crystal events that DON'T have a corresponding TOF time.
	bool m_skipSignalEvents;
	//Use or ignore crystal events that HAVE a corresponding TOF time. Used in case empty waveforms are necessary.
	bool m_skipNoiseEvents;
	bool m_timingValid;


private:
};


}//testbeam namespace
}//myFuncs namespace
