#pragma once

// STL
#include <string>
#include <unordered_map>
#include <vector>

// Root

// ToDo - Find why I need to do getEntries in order to use the tree
// ToDo - Change all TChain functions to NOT use pointers

class TChain;
class TFile;
class TTree;

namespace myFuncs {
//--------------------------------------------------------------------------------------------
// openChain
//********************************************************************************************
// Create a new TChain pointer.
// The files corresponding to fileNamesString are added to the chain. fileNamesString can have wildcards
//<b>All</b> files are expected to contain a tree treeName
// NOTE - This also disables all branches (sets all branch statuses to 0). Each branch that needs to be read needs to be
// specifically enabled (set branch status to 1)
//--------------------------------------------------------------------------------------------
TChain* openChain(const std::string& treeName, const std::string& fileNamesString);

//--------------------------------------------------------------------------------------------
// setChainBranch
//********************************************************************************************
// Add branch branchName to chain with variable pointed at by pointer.
// Also set branch status (and sub branch) to 1 so that the branch can be read.
// There is no check of what pointer points to. The pointer should probably point to zero - this seems to work best.
//--------------------------------------------------------------------------------------------
int setChainBranch(TChain* chain, const std::string& branchName, void* pointer);

//--------------------------------------------------------------------------------------------
// setChainBranch
//********************************************************************************************
// Overloaded!
// Recieve a chain and vectors of branch names and pointers to variables.
// Set branch address for each branch name to the corresponding variable.
// Also set branch (and sub branch) status to 1 so that the branches can be read.
//--------------------------------------------------------------------------------------------
int setChainBranch(TChain* chain, const std::vector<std::string>& branchNamesV, std::vector<void*> pointerV);

//--------------------------------------------------------------------------------------------
// openChain_setBranch
//********************************************************************************************
// Create chain with tree treeName. Add file(s) corresponding to fileNameString(allows wildcards).
// Then set branch branchName with variable pointed at by pointer.
// There are no sanity checks on pointer.
// All other branche statuses in treeName are set to 0(disabled) - so they won't be read by GetEntries...
//--------------------------------------------------------------------------------------------
TChain* openChain_setBranch(const std::string& fileNameString, const std::string& treeName, const std::string& branchName,
                            void* pointer);

//--------------------------------------------------------------------------------------------
// openChain_setBranch - Overloaded
//********************************************************************************************
// Overloaded function - Create chain with tree treeName. Add file(s) corresponding to the elements in fileNameStringV(allows
// wildcards).  Then set branch branchName with variable pointed at by pointer.  There are no sanity checks on pointer.  All other
// branche statuses in treeName are set to 0(disabled) - so they won't be read by GetEntries...
//
// The first element in fileNameStringV has to have a tree treeName and branch branchName. There is no check(I think) that the
// other files have them...
//--------------------------------------------------------------------------------------------
TChain* openChain_setBranch(const std::vector<std::string>& fileNameStringV, const std::string& treeName,
                            const std::string& branchName, void* pointer);

//--------------------------------------------------------------------------------------------
// openChain_setBranch
//********************************************************************************************
// Overloaded!
// Create chain with tree treeName. Add file(s) corresponding to fileNameString(allows wildcards).
// Then set branches in branchNamesV with variables pointed at by pointers in pointerV.
// There are no sanity checks on the pointers.
// All other branch statuses in treeName are set to 0(disabled) - so they won't be read by GetEntries...
//--------------------------------------------------------------------------------------------
TChain* openChain_setBranch(const std::string& fileNameString, const std::string& treeName,
                            const std::vector<std::string>& branchNamesV, std::vector<void*> pointerV);

//--------------------------------------------------------------------------------------------
// openChain_setBranch - Overloaded
//********************************************************************************************
// Overloaded function - Create chain with tree treeName. Add file(s) corresponding to the elements in fileNameStringV(allows
// wildcards).  Then set branches in branchNamesV with variables pointed at by pointers in pointerV.  There are no sanity checks
// on pointer.  All other branche statuses in treeName are set to 0(disabled) - so they won't be read by GetEntries...
//
// The first element in fileNameStringV has to have a tree treeName and branch branchName. There is no check(I think) that the
// other files have them...
//--------------------------------------------------------------------------------------------
TChain* openChain_setBranch(const std::vector<std::string>& fileNameStringV, const std::string& treeName,
                            const std::vector<std::string>& branchNamesV, std::vector<void*> pointerV);

//--------------------------------------------------------------------------------------------
// getFilesEndingWith
//********************************************************************************************
// input:
// std::string dirString - directory to search for files
// std::string endsWith  - extention of filename
// output:
// std::vector<std::string> - list of full file names

// The function searches in dirString for files ending with endswith and returns all of them in a vector of std::strings.

/// Note - relies on dirent.h!!
//--------------------------------------------------------------------------------------------
std::vector<std::string> getFilesEndingWith(const std::string& dirString, const std::string& ending);

// Open a file and check that it's exists/ not zombie.
TFile* openFile(const std::string& filename, const std::string& options);
//   TFile getTFile(const std::string &fullFileName, const std::string &options);

std::ifstream openIfstream(const std::string& filename, const int numHeaders = 0);

//--------------------------------------------------------------------------------------------
// readParamFile
//********************************************************************************************
// Read paramter file.
// return string of pairs <paramName, paramValue>
// File is formated: paramName<whitespace>paramValue
std::unordered_map<std::string, std::string> getParams(const std::string& filename, const int numHeaders = 0);

//--------------------------------------------------------------------------------------------
// readFile
//********************************************************************************************
// Read a file.
// Treat each line as a string.
// Skip first numHeaders lines.
std::vector<std::string> readFile(const std::string& filename, const int numHeaders = 0);

void copyTTreeFromFileToFile(const std::string& treename, const std::string& fromPattern, const std::string& toFilename);

std::vector<std::string> myGlob(const std::string& pattern);
} // namespace myFuncs
