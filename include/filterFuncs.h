#include <vector>

namespace myFuncs {
namespace DSP {

//----------------------------------------------------------
// filterCR_RC

// Do a digital CR-RC filter.

// Inputs:
//  xs - input values
//  tau - filter time constant
//  T - Sampling time

//----------------------------------------------------------
template <typename Type>
std::vector<double> filterCR_RC(const std::vector<Type> &xs, const double tau = 1., const double T = 0.);

//----------------------------------------------------------
// filterCR_RC2

// Do a digital CR-RC^2 filter.

// Inputs:
//  xs - input values
//  tau - filter time constant
//  T - Sampling time

//----------------------------------------------------------
template <typename Type>
std::vector<double> filterCR_RC2(const std::vector<Type> &xs, const double tau = 1., const double T = 0.);

//----------------------------------------------------------
// filterCR_RC4

// Do a digital CR-RC^4 filter.

// Inputs:
//  xs - input values
//  tau - filter time constant
//  T - Sampling time

//----------------------------------------------------------
template <typename Type>
std::vector<double> filterCR_RC4(const std::vector<Type> &xs, const double tau = 1., const double T = 0.);
} // namespace DSP

} // namespace myFuncs
