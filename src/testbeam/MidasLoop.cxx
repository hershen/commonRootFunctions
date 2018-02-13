#include "testbeam/MidasLoop.h"

// std
#include <cassert>

// Mine
#include "fileFuncs.h"
#include "testbeam/testbeamFuncs.h"

using namespace myFuncs::testbeam;

MidasLoop::MidasLoop(const std::string midasFilesPath, const int runNum)
    : TRootanaEventLoop(), m_maxAmplitude(25.0), m_maxEntries(0), m_runNum(0), m_filesToSkip(0), m_maxNumberOfFiles(0),
      m_skipSignalEvents(false), m_skipNoiseEvents(true) {
  // Disable root file output
  DisableRootOutput();

  setRunNum(runNum);

  // Get midas files
  m_midasFilenames = getFilesRelatedToRun(midasFilesPath, getRunNum(), ".mid.gz");

  // Sort by filename so that event numbers are ordered
  std::sort(m_midasFilenames.begin(), m_midasFilenames.end());

  std::cout << "Info::: MidasLoop::MidasLoop: Found " << m_midasFilenames.size() << " .mid.gz files in " << midasFilesPath
            << std::endl;
}

bool MidasLoop::PreFilter(TDataContainer& dataContainer) {

  // Only process eventId = 1 - others are 15,100 - don't know what they are.
  if (dataContainer.GetMidasEvent().GetEventId() != 1)
    return false;

  return true;
}

void MidasLoop::run(const std::string options) {

  if (m_midasFilenames.size() > 99) {
    std::cout << "Can't process more than 99 files. Exisitng" << std::endl;
    return;
  }

  // Reserve array of filenames
  constexpr int maxArgs = 100;
  std::array<char*, maxArgs> argv;

  // Populate array. array.fill is not good because all elements will point to the same place.
  constexpr int maxCharPerFile = 200;
  for (uint iElement = 0; iElement < argv.size(); ++iElement)
    argv[iElement] = new char[maxCharPerFile];

  // Fill first entry with dummy - because ExecuteLoop expects the filename to be there
  strcpy(argv[0], "dummy");

  const bool isOptions = !options.empty();

  // Fill options
  if (isOptions)
    strcpy(argv[1], const_cast<char*>(options.c_str()));

  // Upper bound on number of files to process
  const int upperBoundNumFiles =
      (m_maxNumberOfFiles != 0 and m_maxNumberOfFiles <= m_midasFilenames.size()) ? m_maxNumberOfFiles : m_midasFilenames.size();

  const int numFilesToProcess = upperBoundNumFiles - m_filesToSkip;
  if (numFilesToProcess <= 0) {
    std::cout << "MidasLoop::run: Asked to process " << numFilesToProcess << ". Abroting!!\n";
    std::cout << "upperBoundNumFiles = " << upperBoundNumFiles << ", m_filesToSkip = " << m_filesToSkip << std::endl;
    return;
  }

  // Warn if processing less files than found
  if (numFilesToProcess < static_cast<int>(m_midasFilenames.size()))
    std::cout << "\nWarning::: MidasLoop::run: Asked to process " << numFilesToProcess << " files, even though found "
              << m_midasFilenames.size() << " in given path\n"
              << std::endl;

  std::cout << "\nINFO: Files to process :\n";
  // Fill array with filenames from m_midasFilenames
  int iFileIdx = isOptions + 1; // idx where to hold first filename in argv
  for (int iFile = m_filesToSkip; iFile < upperBoundNumFiles; ++iFile) {
    std::cout << m_midasFilenames[iFile] << "\n";
    if (m_midasFilenames[iFile].size() > maxCharPerFile) {
      std::cout << "m_midasFilenames[" << iFile << "] = " << m_midasFilenames[iFile] << " > maxCharPerFile. Ignoring file"
                << std::endl;
      throw;
    }
    strcpy(argv[iFileIdx], const_cast<char*>(m_midasFilenames[iFile].c_str()));
    ++iFileIdx;
  }

  std::cout << std::endl;

  // Eecute event loop
  std::cout << "Executing loop with " << numFilesToProcess << " files." << std::endl;

  // 	for(int i = 0 ; i < 10; ++i)
  // 		std::cout << "argv[" << i << "] = " << argv[i] << std::endl;

  ExecuteLoop(numFilesToProcess + isOptions + 1, argv.data());

  // Release memory
  for (uint iElement = 0; iElement < argv.size(); ++iElement)
    delete argv[iElement];
}
