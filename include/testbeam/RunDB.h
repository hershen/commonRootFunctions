#pragma once

#include <unordered_map>
#include <vector>

// Mine
#include "testbeam/constants.h"

namespace myFuncs {
namespace testbeam {

struct pedestalFitParams {
  const double p0;
  const double p1;
  const double p2;
};

//-----------------------------------------------------------
// Class representing one TB run and all associated parameters.
//-----------------------------------------------------------
class RunParams {
public:
  constexpr RunParams(const int runNum, const pedestalFitParams &fitParams_channel1, const pedestalFitParams &fitParams_channel15,
                      const pedestalFitParams &fitParams_summedChannels, const Crystal crystal, const double sourceDistance,
                      const double HV, const double nominalBeamMomentum, const double measuredBeamMomentum,
                      const double crystalFrontFaceToIncubatorSideWallDistance)
      : m_runNum(runNum), m_fitParams_channel1(fitParams_channel1), m_fitParams_channel15(fitParams_channel15),
        m_fitParams_summedChannels(fitParams_summedChannels), m_crystal(crystal), m_sourceDistance(sourceDistance), m_HV(HV),
        m_nominalBeamMomentum(nominalBeamMomentum), m_measuredBeamMomentum(measuredBeamMomentum),
        m_crystalFrontFaceToIncubatorSideWallDistance(crystalFrontFaceToIncubatorSideWallDistance) {}

  int getRunNum() const { return m_runNum; }

  Crystal getCrystal() const { return m_crystal; }

  // Negative means no source
  double getSourceDistance() const { return m_sourceDistance; }

  double getHV() const { return m_HV; }

  // Returns the nominal beam momentum - i.e. the momentum we thought the runs were at.
  //[MeV/c]. Negative means negative charged particles in beam
  double getNominalBeamMomentum() const { return m_nominalBeamMomentum; }

  // Returns the measured momentum from the analysis of the TOF data.
  //[MeV/c]. Negative means negative charged particles in beam
  double getMeasuredBeamMomentum() const { return m_measuredBeamMomentum; }

  // Returns the downstream TOF center to crystal center distance in meters.
  constexpr double getDownstream2crystalCenterDistance() const {
    return c_downstreamCenter2incubatorWall + c_incubatorWallSideWidth +
           m_crystalFrontFaceToIncubatorSideWallDistance + // different per run
           0.5 * c_crystalLength;
  }

  pedestalFitParams getPedestalFitParams(const int channel) const {
    if (channel == 1) {
      return m_fitParams_channel1;
    } else if (channel == 15) {
      return m_fitParams_channel15;

    } else if (channel == 16) {
      return m_fitParams_summedChannels;
    }

    return pedestalFitParams{};
  }

private:
  int m_runNum;
  const pedestalFitParams m_fitParams_channel1;       // Fit parameters for mean pedestal [ADC counts] as function of event number
  const pedestalFitParams m_fitParams_channel15;      // Fit parameters for mean pedestal [ADC counts] as function of event number
  const pedestalFitParams m_fitParams_summedChannels; // Fit parameters for mean pedestal [ADC counts] as function of event number
  Crystal m_crystal;
  double m_sourceDistance;       // Negative means no source.
  double m_HV;                   // [V]
  double m_nominalBeamMomentum;  // [MeV/c]. Negative means negative charged particles in beam.
  double m_measuredBeamMomentum; // [MeV/c]. Momentum measured after analyzing TOF data. Negative means negative charged particles
                                 // in beam.
  double m_crystalFrontFaceToIncubatorSideWallDistance; // meters
};

//-----------------------------------------------------------
// Singleton class to access DB
//-----------------------------------------------------------
class RunDB {
public:
  //-----------------------------------------------------------
  // Get instance to singleton
  //-----------------------------------------------------------
  static RunDB &instance() {
    static RunDB instance; // Guaranteed to be destroyed.
                           // Instantiated on first use.
    return instance;
  }

private:
  RunDB();

public:
  RunDB(RunDB const &) = delete;
  void operator=(RunDB const &) = delete;

  // access RunParams with run number runNum
  const RunParams &operator[](const int runNum) const;

private:
  // This holds the DB
  const std::unordered_map<int, const RunParams> m_DB;
};

} // namespace testbeam
} // namespace myFuncs
