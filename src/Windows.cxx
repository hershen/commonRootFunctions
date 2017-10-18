#include "Windows.h"

// ROOT
#include "TMath.h"

namespace myFuncs {

WindowBase::WindowBase(const size_t windowLength) : m_windowLength(windowLength) {}

std::vector<double> WindowBase::getWindowValues() const {

  // create output vector
  std::vector<double> windowValues;
  windowValues.reserve(getWindowLength() + 1);

  // 			//Create window function
  // 			const T windowFunction(windowLength);

  for (size_t index = 0; index < windowValues.capacity(); ++index)
    windowValues.push_back(eval(index));

  return windowValues;
}

//-----------------------------------------------------------
// Hamming
//-----------------------------------------------------------

// Constructor
Hamming::Hamming(const size_t windowLength) : WindowBase(windowLength) {}

// Taken fromDiscrete time signal processing Oppenheim and Schafer, 3rd edition, page 536
double Hamming::eval(const long index) const {

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wsign-compare"
  if (index < 0 or index > getWindowLength())
    return 0.0;
#pragma GCC diagnostic pop

  return 0.54 - 0.46 * TMath::Cos(2.0 * TMath::Pi() * index / getWindowLength());
}

//-----------------------------------------------------------
// Hann
//-----------------------------------------------------------

// Constructor
Hann::Hann(const size_t windowLength) : WindowBase(windowLength) {}

// Taken fromDiscrete time signal processing Oppenheim and Schafer, 3rd edition, page 536
double Hann::eval(const long index) const {

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wsign-compare"
  if (index < 0 or index > getWindowLength())
    return 0.0;
#pragma GCC diagnostic pop

  return 0.5 - 0.5 * TMath::Cos(2.0 * TMath::Pi() * index / getWindowLength());
}
} // namespace myFuncs
