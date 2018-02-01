#pragma once

// STL
#include <iostream>

// ROOT
#include "TMatrixD.h"
#include "TPrincipal.h"

namespace myFuncs {

//-------------------------------------
// Class to calculate Covariance matrix on the run.
//-------------------------------------
class RunningCovMatrix : public TPrincipal {
public:
  // Delete empty Constructor
  RunningCovMatrix() = delete; // Empty constructor of TPrincipal messes things up.

  RunningCovMatrix(const int numVariables) : TPrincipal(numVariables, "N") {
    if (numVariables < 2) {
      throw std::invalid_argument("RunningCovMatrix::RunningCovMatrix(" + std::to_string(numVariables) +
                                  ") - numVariables needs > 1");
    }
  };

  //TPrincipal only accepts doubles
  void addVector(const std::vector<double> vector) {

    if (static_cast<int>(vector.size()) != GetCovarianceMatrix()->GetNcols()) {
      std::cerr << "vector.size() = " << vector.size() << " != matrix size  = " << GetCovarianceMatrix()->GetNcols() << "x"
                << GetCovarianceMatrix()->GetNrows() << ". Not adding vector!!!\n";
      return;
    }
    AddRow(vector.data());
  }

  // Return num elements that were added
  inline size_t getNumVariables() const { return fNumberOfVariables; }

  // Return num data points that were added.
  // Each vector is a datapoint
  inline size_t getNumDataPoints() const { return fNumberOfDataPoints; }

  inline TMatrixD getCovarianceMatrix() const { return *this->GetCovarianceMatrix(); }
};

} // namespace myFuncs
