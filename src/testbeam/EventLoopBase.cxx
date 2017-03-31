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
m_maxNumberOfFiles(0),
m_TOFeventNumber(std::make_shared<Long64_t>(-1)),
m_beta(std::make_shared<double>(0.0)),
m_betaError(std::make_shared<double>(0.0)),
m_timeAtCrystal_ns(std::make_shared<double>(0.0)),
m_timeAtCrystalError_ns(std::make_shared<double>(0.0)),
m_skipMissingTimeEvents(true),
m_skipPresentTimeEvents(false)
{
	//Disable root file output
	DisableRootOutput();
		
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
	if(m_maxEntries > 0) m_timingChain->GetEntry(m_timingEntry);
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
	assert(m_timingEntry < m_maxEntries);
	
	//If not required to skip, processes this event
	if( !getSkipMissingTimeEvents() ) return true;	
	
 	return *m_TOFeventNumber == midasEventNum;
	
}

void EventLoopBase::run(const std::string options) {
	
	if(m_midasFilenames.size() > 99) {
		std::cout << "Can't process more than 99 files. Exisitng" << std::endl;
		return;
	}
	
	//Reserve array of filenames
	constexpr int maxArgs = 100;	
	std::array<char*, maxArgs> argv;
	
	
	//Populate array. array.fill is not good because all elements will point to the same place.
	constexpr int maxCharPerFile = 200;
	for(uint iElement = 0; iElement < argv.size(); ++iElement)
		argv[iElement] = new char[maxCharPerFile];
	
	//Fill first entry with dummy - because ExecuteLoop expects the filename to be there
	strcpy(argv[0], "dummy");
	
	const bool isOptions = !options.empty();
	
	//Fill options
	if(isOptions)
		strcpy(argv[1],const_cast<char*>( options.c_str() ));
	
	//Max file to process
	const unsigned int maxFiles = ( m_maxNumberOfFiles != 0 and m_maxNumberOfFiles <= m_midasFilenames.size() ) ? m_maxNumberOfFiles : m_midasFilenames.size() ;
	
	//Warn if processing less files than found
	if(maxFiles < m_midasFilenames.size())
		std::cout << "\nWarning::: EventLoopBase::run: Processing " << maxFiles << " even though found " << m_midasFilenames.size() << " in given path\n" << std::endl;
	
	//Fill array with filenames from m_midasFilenames
	for(uint iFile = 0; iFile < maxFiles; ++iFile) {
		if(m_midasFilenames[iFile].size() > maxCharPerFile) {
			std::cout << "m_midasFilenames[" << iFile << "] = " << m_midasFilenames[iFile] << " > maxCharPerFile. Ignoring file" << std::endl;
			throw;
		}
		strcpy(argv[iFile+isOptions+1], const_cast<char*>( m_midasFilenames[iFile].c_str() ) );
	}
	
// 	for(auto s : argv) std::cout << s << std::endl;
	
	std::cout << "Executing loop with " << maxFiles << " files." << std::endl;
	
	//Eecute event loop
	if(maxFiles > 0) {
		ExecuteLoop(maxFiles + isOptions + 1, argv.data());
	}

	//Release memory
	for (uint iElement = 0; iElement < argv.size(); ++iElement) 
		delete argv[iElement];

}