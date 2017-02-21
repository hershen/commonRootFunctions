#include "testbeam/EventLoopBase.h"

//std
#include <cassert>

//Root
#include "TChain.h"

//Mine
#include "testbeam/testbeamFuncs.h"
#include "fileFuncs.h"

using namespace myFuncs::testbeam;

EventLoopBase::EventLoopBase(const std::string midasFilesPath, const std::string TOFfilesPath, const int runNum):
TRootanaEventLoop(),
m_runNum(0),
m_timingEntry(0),
m_maxEntries(0),
m_TOFeventNumber(std::make_shared<Long64_t>(-1)),
m_beta(std::make_shared<double>(0.0)),
m_betaError(std::make_shared<double>(0.0)),
m_timeAtCrystal_ns(std::make_shared<double>(0.0)),
m_timeAtCrystalError_ns(std::make_shared<double>(0.0)),
m_skipMissingTimeEvents(true)
{
	setRunNum(runNum);
	
	//Setup the TOF chain
	setupTOFchain(TOFfilesPath);
	
	//Get midas files
	m_midasFilenames = getFilesRelatedToRun(midasFilesPath, getRunNum(), ".mid.gz");

	//Sort by filename so that event numbers are ordered
	std::sort( m_midasFilenames.begin(), m_midasFilenames.end() );
	
	std::cout << "Info::: EventLoopBase::EventLoopBase: Found " << m_midasFilenames.size() << " .mid.gz files in " << midasFilesPath << std::endl;
}

void EventLoopBase::setupTOFchain(const std::string TOFfilesPath) {
	
	//Get files with timing information
	const auto timingFilenames = getFilesRelatedToRun(TOFfilesPath, getRunNum(), ".root");
	
	//Print warning if there is more than one file
	if( timingFilenames.size() > 1 )
		std::cout << "Warning::: EventLoopBase::EventLoopBase: Found " << timingFilenames.size() << " files for run " << getRunNum() << " in path " << TOFfilesPath << ". Expecting 1!!" << std::endl;

	const std::vector<std::string> branchNames = {"eventNumber", "beta", "betaError", "timeAtCrystal_ns", "timeAtCrystalError_ns"};
	
	//Define vector of pointer addresses
	const std::vector<void*> pointers = {m_TOFeventNumber.get(), m_beta.get(), m_betaError.get(), m_timeAtCrystal_ns.get(), m_timeAtCrystalError_ns.get()};
	
	//Create chain of files
	m_timingChain = std::shared_ptr<TChain>(myFuncs::openChain_setBranch(timingFilenames , "tree", branchNames, pointers));
	
	m_timingEntry = 0;
	m_maxEntries = m_timingChain->GetEntries();
	
}

void EventLoopBase::Initialize(void) {
	if(m_maxEntries > 0) m_timingChain->GetEntry(0);
}

bool EventLoopBase::PreFilter(TDataContainer& dataContainer) {
	
	//Only process eventId = 1 - others are 15,100 - don't know what they are.
	if(dataContainer.GetMidasEvent().GetEventId() != 1) return false;
	
	//Advance timing chain until timing event number >= midas event number
	//Dangerous - no bounds checking!!! chain might overflow
	uint32_t midasEventNum = dataContainer.GetMidasEvent().GetSerialNumber();
	while(*m_TOFeventNumber < midasEventNum) {
		m_timingChain->GetEntry(++m_timingEntry);
	}
	
	///TODO get rid of assert
	assert(*m_TOFeventNumber < m_maxEntries);
	
	//If not required to skip, processes this event
	if( !skipMissingTimeEvents() ) return true;	
	
 	return *m_TOFeventNumber == midasEventNum;
	
}

void EventLoopBase::run() {
	
	if(m_midasFilenames.size() > 99) {
		std::cout << "Can't process more than 99 files. Exisitng" << std::endl;
		return;
	}
	
	//Reserve array of filenames
	std::array<char*, 100> argv;
	argv.fill(new char[200]);
	
	//Fill first entry with dummy - because ExecuteLoop expects the filename to be there
	strcpy(argv[0], "dummy");
	
	//Fill array with filenames from m_midasFilenames
	for(uint iFile = 0; iFile < m_midasFilenames.size(); ++iFile) {
		if(m_midasFilenames[iFile].size() > 200) {
			std::cout << "m_midasFilenames[" << iFile << "] = " << m_midasFilenames[iFile] << " > 200. Ignoring file" << std::endl;
			continue;
		}
		strcpy(argv[iFile+1], const_cast<char*>( m_midasFilenames[iFile].c_str() ) );
	}
	
	std::cout << "Executing loop with " << m_midasFilenames.size() << " files." << std::endl;
	
	//Eecute event loop
	if(m_midasFilenames.size() > 0) {
		ExecuteLoop(m_midasFilenames.size() + 1, argv.data());
	}
	//Release memory
// 	//TODO - find out why this gives a double free or corruption
// 	for (uint iElement = 0; iElement < m_midasFilenames.size() + 1; ++iElement) 
// 	{
// 		std::cout << iElement << std::endl;
// 		std::cout << "line = " << __LINE__ << std::endl;
// 		delete argv[iElement];
// 		std::cout << "line = " << __LINE__ << std::endl;
// 	}
}