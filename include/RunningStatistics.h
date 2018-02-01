#pragma once

#include <cmath>
#include <iostream>

namespace myFuncs {

//-------------------------------------
// Class to calculate statistics on the run, without memory of previous values
// Ref - The art of computer programming, Knuth, Vol2, Second edition (1981), p. 216
//-------------------------------------
class RunningStatistics {
public:
  // Constructor
  RunningStatistics() = default;

  // Return num elements that were added
  inline size_t getNumElements() const { return m_numElements; }

  // Return sample mean
  inline double getSampleMean() const { return m_newM; }

  // Return sample variance
  inline double getSampleVariance() const { return getNumElements() > 1 ? (m_newS / (getNumElements() - 1)) : 0.0; }

  // Return sample standard deviation
  inline double getSampleStd() const { return std::sqrt(getSampleVariance()); }

  template <class T>
  void addElement(const T element) {
    ++m_numElements;
    m_oldM = m_newM;
    m_newM = m_oldM + (element - m_oldM) / m_numElements;

    m_oldS = m_newS;
    m_newS = m_oldS + (element - m_oldM) * (element - m_newM);
  }

private:
  size_t m_numElements = 0;
  double m_oldM;
  double m_oldS;
  double m_newM = 0.0;
  double m_newS = 0.0;
};

} // namespace myFuncs
