#include "eclCrystalDB.h"

#include "TChain.h"
#include "TMath.h"
#include "fileFuncs.h"

#include <iostream>

#define treeName "tree"

#define cellIdBranchName "cellId"
#define thetaIdBranchName "thetaId"
#define phiIdBranchName "phiId"
#define positionXbranchName "positionX"
#define positionYbranchName "positionY"
#define positionZbranchName "positionZ"
#define directionXbranchName "directionX"
#define directionYbranchName "directionY"
#define directionZbranchName "directionZ"

eclCrystalDB* eclCrystalDB::m_eclCrystalDB_Instance = nullptr;

eclCrystalDB* eclCrystalDB::instance() {
  if (!m_eclCrystalDB_Instance) // Only allow one instance of class to be generated.
    m_eclCrystalDB_Instance = new eclCrystalDB;
  return m_eclCrystalDB_Instance;
}

eclCrystalDB::eclCrystalDB(std::string fullFileName) : maxThetaId(-1) {
  int cellIdPointer = 0;
  int thetaIdPointer = 0;
  int phiIdPointer = 0;
  double positionXpointer = 0;
  double positionYpointer = 0;
  double positionZpointer = 0;
  double directionXpointer = 0;
  double directionYpointer = 0;
  double directionZpointer = 0;

  std::vector<std::string> branchNames = {cellIdBranchName,     thetaIdBranchName,    phiIdBranchName,
                                          positionXbranchName,  positionYbranchName,  positionZbranchName,
                                          directionXbranchName, directionYbranchName, directionZbranchName};

  std::vector<void*> pointers = {&cellIdPointer,    &thetaIdPointer,    &phiIdPointer,      &positionXpointer, &positionYpointer,
                                 &positionZpointer, &directionXpointer, &directionYpointer, &directionZpointer};

  TChain* inputChain = myFuncs::openChain_setBranch(fullFileName, std::string(treeName), branchNames, pointers);

  for (auto entry = 0; entry < inputChain->GetEntries(); ++entry) {
    inputChain->GetEntry(entry);
    crystalInfo tempCrystal{cellIdPointer,    thetaIdPointer,    phiIdPointer,      positionXpointer, positionYpointer,
                            positionZpointer, directionXpointer, directionYpointer, directionZpointer};
    crystalsDB.push_back(tempCrystal);
  }

  delete inputChain;
}

eclCrystalDB::~eclCrystalDB(void) { std::cout << "in destructor" << std::endl; }

void eclCrystalDB::getThetaPhi(const int& cellID, double& theta, double& phi) const {
  theta = eclCrystalDB::getTheta(cellID);
  phi = eclCrystalDB::getPhi(cellID);
}
//
double eclCrystalDB::getPhi(const int& cellID) const { return getPosition(cellID).Phi(); }

double eclCrystalDB::getTheta(const int& cellID) const { return getPosition(cellID).Theta(); }

void eclCrystalDB::getThetaIDPhiID(const int& cellID, int& thetaID, int& phiID) const {
  thetaID = eclCrystalDB::getThetaID(cellID);
  phiID = eclCrystalDB::getPhiID(cellID);
}

int eclCrystalDB::getThetaID(const int& cellID) const { return crystalsDB[cellID].thetaID; }

int eclCrystalDB::getPhiID(const int& cellID) const { return crystalsDB[cellID].phiID; }

TVector3 eclCrystalDB::getPosition(const int& cellID) const {
  double positionX = crystalsDB[cellID].positionX;
  double positionY = crystalsDB[cellID].positionY;
  double positionZ = crystalsDB[cellID].positionZ;
  return TVector3(positionX, positionY, positionZ);
}

TVector3 eclCrystalDB::getDirection(const int& cellID) const {
  double directionX = crystalsDB[cellID].directionX;
  double directionY = crystalsDB[cellID].directionY;
  double directionZ = crystalsDB[cellID].directionZ;
  return TVector3(directionX, directionY, directionZ);
}

double eclCrystalDB::getCylindricalR(const int& cellID) const {
  double x = getPositionX(cellID);
  double y = getPositionY(cellID);
  return TMath::Sqrt(std::pow(x, 2) + std::pow(y, 2));
}

int eclCrystalDB::getMaxThetaId() const {
  // If maximum has already been found, return it
  if (maxThetaId >= 0)
    return maxThetaId;

  // Else, look for it
  for (auto crystal : crystalsDB) {
    if (crystal.thetaID > maxThetaId)
      maxThetaId = crystal.thetaID;
  }

  return maxThetaId;
}

int eclCrystalDB::getRegion(const int& cellId) const {
  double theta = getTheta(cellId);
  if (theta < FEthetaMax)
    return 1; // forward endcap
  if (theta > BEthetaMin)
    return 3; // backward endcap
  else
    return 2; // barrel
}
