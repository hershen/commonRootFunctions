#pragma once

#include <map>
#include <unordered_map>

namespace myFuncs {
namespace testbeam {

enum class Crystal { CsI_Tl_Belle, CsI_Tl_Babar, CsI_Ukrainian, CsI_Chinese };

constexpr double c_upstreamTOF2S0centersDistance = 0.44; // in meters #(average of 43.43 cm, 44.56 cm)

// old version - I don't know where the 43.1 is from...
// It came from a (probably incorrect calculation) done using elog 11 (appears in oneNote note).
// constexpr double c_upstreamTOF2S0distance = 0.437; //in meters #(average of 43.1 cm, 43.43 cm, 44.56 cm)

constexpr double c_downstreamTOF2upstreamTOFdistance = 3.092; // in meters
constexpr double c_downstreamTOF2S0distance = c_downstreamTOF2upstreamTOFdistance - c_upstreamTOF2S0centersDistance;

constexpr double c_incubatorWallSideWidth = 0.0525;         // meters
constexpr double c_downstreamCenter2incubatorWall = 0.0419; // meters

constexpr double c_downstreamFace2incubatorWall = 0.0034; // meters (based on pic 20150811_164604)

constexpr double c_upstreamLeadCollimatorFacesDistance = 0.02; // meters. Based on M11distances.pptx

constexpr double c_crystalLength = 0.3; // meters

// I.e., to align channel1 and channel 15 in time, channel 1 needs to be moved 1.088 ns to the right.
// Based on Timing resolution CR-RC4 + HPF
constexpr double c_channel1_channel15_timeDifference = 1.088; // ns.

using CrystalAdc2meV = std::unordered_map<int, double>;
const static CrystalAdc2meV c_CsI_Tl_Belle_Adc2MeV{{1, 0.7073}, {15, 0.7209}}; // [MeV / ADC counts]
const static std::map<Crystal, CrystalAdc2meV> c_crystal2_adc2mev{
    {Crystal::CsI_Tl_Belle, c_CsI_Tl_Belle_Adc2MeV}}; // Unordered map doesn't currently work with enums. Maybe with c++14

// maps of [beam momentum] = value.
// For mean and sigma.
// Values calculated from fitting a gaussian to around mean +- 2 sigma in beta distribtuion.
// Done in script findTOFcrystalTimesAndSigmas.C
const static std::unordered_map<int, double> muonBetaMean = {{-100, 0.72961}, {-120, 0.786136}, {-140, 0.828577}};
const static std::unordered_map<int, double> pionBetaMean = {{-100, 0.605955}, {-120, 0.680284}, {-140, 0.741776}};

const static std::unordered_map<int, double> electronBetaSigma = {{-100, 0.028034}, {-120, 0.02957}, {-140, 0.028172}};
const static std::unordered_map<int, double> muonBetaSigma = {{-100, 0.013659}, {-120, 0.015451}, {-140, 0.017359}};
const static std::unordered_map<int, double> pionBetaSigma = {{-100, 0.010849}, {-120, 0.011701}, {-140, 0.013979}};

const static std::unordered_map<int, double> c_max_min_maximumValue{{1, 37}, {15, 40}};
} // namespace testbeam
} // namespace myFuncs
