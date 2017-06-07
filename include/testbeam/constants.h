#pragma once

#include <unordered_map>

namespace myFuncs {
namespace testbeam{
	
constexpr double c_upstreamTOF2S0distance = 0.437; //in meters #(average of 43.1 cm, 43.43 cm, 44.56 cm)
constexpr double c_downstreamTOF2upstreamTOFdistance = 3.092; //in meters
constexpr double c_downstreamTOF2S0distance = c_downstreamTOF2upstreamTOFdistance - c_upstreamTOF2S0distance;

constexpr double c_incubatorWallSideWidth = 0.0525; // meters
constexpr double c_downstreamCenter2incubatorWall = 0.0419; // meters

constexpr double c_crystalLength = 0.3; // meters


//maps of [beam momentum] = value.
//For mean and sigma.
//Values calculated from fitting a gaussian to around mean +- 2 sigma in beta distribtuion.
//Done in script findTOFcrystalTimesAndSigmas.C
const std::unordered_map<int, double> muonBetaMean = { {-100, 0.72961},{-120, 0.786136}, {-140, 0.828577} }; 
const std::unordered_map<int, double> pionBetaMean = { {-100, 0.605955},{-120, 0.680284}, {-140, 0.741776} }; 

const std::unordered_map<int, double> electronBetaSigma = { {-100, 0.028034},{-120, 0.02957}, {-140, 0.028172} }; 
const std::unordered_map<int, double> muonBetaSigma = { {-100, 0.013659},{-120, 0.015451}, {-140, 0.017359} }; 
const std::unordered_map<int, double> pionBetaSigma = { {-100, 0.010849},{-120, 0.011701}, {-140, 0.013979} }; 


}
}