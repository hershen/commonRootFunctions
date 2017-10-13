#include <vector>

namespace myFuncs {
namespace DSP {

// T is sampling interval
std::vector<double> filterCR_RC(std::vector<double> xV, double tau = 1., double T = 0.);

template <typename Type>
std::vector<double> filterCR_RC2(std::vector<Type> xV, double tau = 1., double T = 0.);

template <typename Type>
std::vector<double> filterCR_RC4(std::vector<Type> xV, double tau = 1., double T = 0.);

} // namespace DSP

} // namespace myFuncs
