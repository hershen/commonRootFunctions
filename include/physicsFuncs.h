#pragma once

//Root
#include "TLorentzVector.h"
#include <cmath>
// #include "TMath.h"

namespace myFuncs
{

//gamma = 1/sqrt(1+beta^2)
inline double beta2gamma(const double beta) {
	return 1./std::sqrt(1.0 - beta * beta);
}

inline double gamma2beta(const double gamma) {
	return std::sqrt(gamma*gamma -1) / gamma;
}

//P = beta * gamma * mass
double beta_mass2momentum(const double beta, const double mass) {
	return beta * beta2gamma(beta) * mass;
}

}