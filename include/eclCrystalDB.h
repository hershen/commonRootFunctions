#pragma once

#include "TVector3.h"
#include <string>
#include <vector>

// Cell ids start at 0?
// Units of distance in cm

class eclCrystalDB {

public:
  static eclCrystalDB *instance();

  inline int getNumCrystals() const { return crystalsDB.size(); }

  void getThetaPhi(const int &cellID, double &theta, double &phi) const;

  double getTheta(const int &cellID) const;

  double getPhi(const int &cellID) const;

  void getThetaIDPhiID(const int &cellID, int &thetaID, int &phiID) const;

  int getThetaID(const int &cellID) const;

  int getPhiID(const int &cellID) const;

  inline double getPositionX(const int &cellID) const { return crystalsDB[cellID].positionX; }

  inline double getPositionY(const int &cellID) const { return crystalsDB[cellID].positionY; }

  inline double getPositionZ(const int &cellID) const { return crystalsDB[cellID].positionZ; }

  double getCylindricalR(const int &cellID) const;

  inline double getCylindricalTheta(const int &cellID) const { return getPhi(cellID); }

  TVector3 getPosition(const int &cellID) const;

  TVector3 getDirection(const int &cellID) const;

  int getMaxThetaId() const;

  // 1 - forward endcap, 2 - Barrel, 3 - backward endcap
  int getRegion(const int &cellId) const;

private:
  eclCrystalDB(std::string fullFileName = "/home/hershen/PhD/ECLcalibration/crystalDB.root");

  ~eclCrystalDB();

  eclCrystalDB(eclCrystalDB const &); // copy constructor is private

  eclCrystalDB &operator=(eclCrystalDB const &); // assignment operator is private

  static eclCrystalDB *m_eclCrystalDB_Instance;

  struct crystalInfo {
    int cellID;
    int thetaID;
    int phiID;
    double positionX;
    double positionY;
    double positionZ;
    double directionX;
    double directionY;
    double directionZ;
  };

  std::vector<crystalInfo> crystalsDB;

  mutable int maxThetaId;

  static constexpr double FEthetaMax = 0.548033;
  static constexpr double BEthetaMin = 2.28115;
};
