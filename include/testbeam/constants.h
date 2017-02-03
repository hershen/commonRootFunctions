#pragma once

namespace myFuncs {
namespace testbeam{
	
constexpr double c_upstreamTOF2S0distance = 0.437; //in meters #(average of 43.1 cm, 43.43 cm, 44.56 cm)
constexpr double c_downstreamTOF2upstreamTOFdistance = 3.092; //in meters
constexpr double c_downstreamTOF2S0distance = c_downstreamTOF2upstreamTOFdistance - c_upstreamTOF2S0distance;

constexpr double c_incubatorWallSideWidth = 0.0525; // meters
constexpr double c_downstreamCenter2incubatorWall = 0.0419; // meters

constexpr double c_crystalLength = 0.3; // meters

}
}