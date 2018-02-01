#include "fileFuncs.h"

#include <iostream>
#include <limits>

// For getParams
#include <fstream>

// ROOT
#include "TChain.h"
#include "TFile.h"
#include "TObjArray.h"
#include "TROOT.h"
#include "TTree.h"
#include "dirent.h"
#include "stringFuncs.h"

// ToDo - Find why I need to do getEntries in order to use the tree
// ToDo - Change all TChain functions to NOT use pointers
// namespace fileFuncs
// {

//   enum {kBigNumber = 1234567890}; // Used for TChain default value

namespace myFuncs {

TChain *openChain(const std::string &treeName, const std::string &fileNamesString) {

  // create chain
  TChain *chain = new TChain(treeName.data());

  // Add files
  int numOfFiles = chain->Add(fileNamesString.data());

  // Check if files were added
  if (!numOfFiles) {
    std::cout << "fileFuncs::openChain: Error - No files found corresponding to string '" << fileNamesString.data()
              << "'. Returning 0." << std::endl;
    delete chain;
    return 0;
  }
  gROOT->cd();

  // Check if tree exists in the files
  if (!chain->GetListOfBranches()) {
    std::cout << "fileFuncs::openChain: Error - No TTree with name '" << treeName.data() << "' in files with string '"
              << fileNamesString.data() << "'. Returning 0." << std::endl;
    delete chain;
    return 0;
  }

  // Check if tree is empty
  if (!(chain->GetListOfBranches()->GetEntriesFast())) {
    std::cout << "fileFuncs::openChain: Error - TTree with name '" << treeName.data() << "' in files with string '"
              << fileNamesString.data() << "' has no branch (it's probably empty). Returning 0." << std::endl;
    //       delete branchArray;
    delete chain;
    return 0;
  }
  chain->SetBranchStatus("*", 0);

  return chain;
}

int setChainBranch(TChain *chain, const std::string &branchName, void *pointer) {
  if (!chain) {
    std::cout << "fileFuncs::setChainBranch: Error - empty TChain given. Aborting." << std::endl;
    return 100;
  }

  // ------------------------------------------------------------------
  // Set branch status (and all sub branches) to 1 so it can be accessed by getEntry(i)
  // ------------------------------------------------------------------
  chain->SetBranchStatus((branchName + "*").data(),
                         1); // TChain documentation says this should be set before the SetBranchAddress()

  // ------------------------------------------------------------------
  // Set branch address
  // ------------------------------------------------------------------
  int returnVal = 0;

  returnVal = chain->SetBranchAddress(branchName.data(), pointer) || returnVal;

  if (returnVal) {
    std::cout << "fileFuncs::setChainBranch: Error - Error while setting branch '" << branchName.data()
              << "' address - check branch name. Aborting" << std::endl;
    return returnVal;
  }

  return returnVal; // should be 0
}

int setChainBranch(TChain *chain, const std::vector<std::string> &branchNamesV, std::vector<void *> pointerV) {
  if (!chain) {
    std::cout << "fileFuncs::setChainBranch: Error - empty TChain given. Aborting." << std::endl;
    return 1000;
  }

  if (branchNamesV.size() != pointerV.size()) {
    std::cout << "fileFuncs::setChainBranch: Error - branchNamesV.size() = " << branchNamesV.size()
              << " != pointerV.size() = " << pointerV.size() << ". Aborting" << std::endl;
    return 1001;
  }

  if (branchNamesV.size() < 1) {
    std::cout << "fileFuncs::setChainBranch: Error - vector sizes different. Aborting" << std::endl;
    return 1002;
  }

  // ------------------------------------------------------------------
  // Set branch addresses and statuses
  // ------------------------------------------------------------------
  int returnVal = 0;
  for (unsigned int idx = 0; idx < branchNamesV.size(); ++idx)
    returnVal = setChainBranch(chain, branchNamesV[idx].data(), pointerV[idx]) || returnVal;

  if (returnVal) {
    std::cout << "fileFuncs::setChainBranch: Error - Error while setting branch address - check branch names. Aborting"
              << std::endl;
    return returnVal;
  }

  return returnVal; // shoud be 0
}

TChain *openChain_setBranch(const std::string &fileNameString, const std::string &treeName, const std::string &branchName,
                            void *pointer) {
  // create chain
  TChain *chain = openChain(treeName, fileNameString);
  if (!chain) {
    std::cout << "fileFuncs::openChain_setBranch: Could not open chain. Aborting" << std::endl;
    return 0;
  }

  // Set branch
  int returnVal = setChainBranch(chain, branchName, pointer);
  if (returnVal) {
    std::cout << "fileFuncs::openChain_setBranch: Could not set branch address. Recieved returnVal = " << returnVal
              << ". Aborting" << std::endl;
    delete chain;
    return 0;
  }

  return chain;
}

TChain *openChain_setBranch(const std::vector<std::string> &fileNameStringV, const std::string &treeName,
                            const std::string &branchName, void *pointer) {

  if (fileNameStringV.empty()) {
    std::cout << "fileFuncs::openChain_setBranch: fileNameStringV is empty. Aborting" << std::endl;
    return 0;
  }

  // ------------------------------------------------------------------
  // Create chain
  // ------------------------------------------------------------------
  TChain *chain = openChain_setBranch(fileNameStringV[0], treeName, branchName, pointer);

  if (!chain) {
    std::cout << "fileFuncs::openChain_setBranch: Error - Could not open chain. Aborting" << std::endl;
    return 0;
  }

  // ------------------------------------------------------------------
  // Add all other files to chain
  // ------------------------------------------------------------------
  for (size_t iFile = 1; iFile < fileNameStringV.size(); ++iFile) {
    int numFilesAdded = chain->Add(fileNameStringV[iFile].data());
    if (!numFilesAdded)
      std::cout << "fileFuncs::openChain_setBranch: Warning - fileNameStringV[" << iFile
                << "] = " << fileNameStringV[iFile].data() << " did not add any files to chain. Continuing." << std::endl;
  }

  return chain;
}

TChain *openChain_setBranch(const std::string &fileNameString, const std::string &treeName,
                            const std::vector<std::string> &branchNamesV, std::vector<void *> pointerV) {
  // create chain
  TChain *chain = openChain(treeName, fileNameString);
  if (!chain) {
    std::cout << "fileFuncs::openChain_setBranch: Could not open chain. Aborting" << std::endl;
    return 0;
  }

  // Set branch
  int returnVal = setChainBranch(chain, branchNamesV, pointerV);
  if (returnVal) {
    std::cout << "fileFuncs::openChain_setBranch: Could not set branch addresses. Recieved returnVal = " << returnVal
              << ". Aborting" << std::endl;
    delete chain;
    return 0;
  }

  return chain;
}

TChain *openChain_setBranch(const std::vector<std::string> &fileNameStringV, const std::string &treeName,
                            const std::vector<std::string> &branchNamesV, std::vector<void *> pointerV) {
  if (fileNameStringV.empty()) {
    std::cout << "fileFuncs::openChain_setBranch: fileNameStringV is empty. Aborting" << std::endl;
    return 0;
  }

  // ------------------------------------------------------------------
  // Create chain
  // ------------------------------------------------------------------
  TChain *chain = openChain_setBranch(fileNameStringV[0], treeName, branchNamesV, pointerV);

  if (!chain) {
    std::cout << "fileFuncs::openChain_setBranch: Error - Could not open chain. Aborting" << std::endl;
    return 0;
  }

  // ------------------------------------------------------------------
  // Add all other files to chain
  // ------------------------------------------------------------------
  for (size_t iFile = 1; iFile < fileNameStringV.size(); ++iFile) {
    int numFilesAdded = chain->Add(fileNameStringV[iFile].data());
    if (!numFilesAdded)
      std::cout << "fileFuncs::openChain_setBranch: Warning - fileNameStringV[" << iFile
                << "] = " << fileNameStringV[iFile].data() << " did not add any files to chain. Continuing." << std::endl;
  }

  return chain;
}

std::vector<std::string> getFilesEndingWith(const std::string &dirString, const std::string &ending) {
  //     #include "dirent.h"
  std::vector<std::string> rootFilesVec;

  DIR *dir;
  if ((dir = opendir(dirString.data())) != NULL) {
    /* print all the files and directories within directory */
    struct dirent *dirEntry;
    while ((dirEntry = readdir(dir)) != NULL) {
      std::string fileNameString(dirEntry->d_name);
      if (myFuncs::endsWith(fileNameString, ending))
        rootFilesVec.push_back(dirString + "/" + fileNameString);
    }
    closedir(dir);
  } else {
    /* could not open directory */
    std::cerr << "fileFuncs::getFilesEndingWith error: could not open directory " << dirString.data() << std::endl;
    return rootFilesVec;
  }

  return rootFilesVec;
}

TFile *openFile(const std::string &filename, const std::string &options) {
  TFile *file = new TFile(filename.data(), options.data());
  if (!file || file->IsZombie()) {
    std::cout << "openFile::Cannot open file  " << filename << ". Aborting" << std::endl;
    return nullptr;
  }
  return file;
}

std::ifstream openIfstream(const std::string &filename, const int numHeaders) {
  std::ifstream file(filename);

  if (!file.is_open()) {
    std::cout << "Could not open stream " << filename << ". Perhaps file doesn't exist???" << std::endl;
  }

  // Read header lines
  for (int lineNum = 0; lineNum < numHeaders; ++lineNum) {
    file.ignore(std::numeric_limits<std::streamsize>::max(), file.widen('\n'));
  }

  return file;
}

std::unordered_map<std::string, std::string> getParams(const std::string &filename, const int numHeaders) {
  std::ifstream file = myFuncs::openIfstream(filename, numHeaders);

  std::string paramName;
  std::string paramValue;

  std::unordered_map<std::string, std::string> params;
  while (file >> paramName >> paramValue) {
    params[paramName] = paramValue;
  }

  return params;
}

std::vector<std::string> readFile(const std::string &filename, const int numHeaders) {
  std::ifstream file = myFuncs::openIfstream(filename, numHeaders);

  std::string line;
  std::vector<std::string> outputs;
  std::unordered_map<std::string, std::string> params;
  while (file >> line) {
    outputs.push_back(line);
  }

  return outputs;
}

void copyTTreeFromFileToFile(const std::string &treename, const std::string &fromFilename, const std::string &toFilename) {
  TFile fromFile(fromFilename.c_str(), "READ");
  TTree *oldtree = (TTree *)fromFile.Get(treename.c_str());

  //Clone tree to new file
  TFile toFile(toFilename.c_str(), "UPDATE");
  oldtree->CloneTree();

  //Cleanup
  toFile.Write();
  fromFile.Close();
  toFile.Close();
}
} // namespace myFuncs
