#pragma once

// std
#include <iostream>

// Eigen
#include <Eigen/Dense>

// Root
#include "TMatrixD.h"
#include "TMatrixDSym.h"

namespace myFuncs {

// not tested
// Returns the contents of a TMatrix in an Eigen SelfadjointView.
// template <unsigned int UpLo = Eigen::Lower>
// Eigen::SelfAdjointView<Eigen::Map<Eigen::MatrixXd>, UpLo> TMatrixToEigenSymmetric(TMatrixD &matrix) {
//   typedef Eigen::Map<Eigen::MatrixXd> MapType;
//   MapType m2map(matrix.GetMatrixArray(), matrix.GetNrows(), matrix.GetNcols());
//   // Eigen::MatrixXd m2(matrix.GetNrows(), matrix.GetNcols());
//
//   // Matrix4d M = Map<Matrix<double,4,4,RowMajor> >(data);
//   return m2map.template selfadjointView<UpLo>();
// }

// Returns the contents of a TMatrix in an Eigen matrix.
//
Eigen::MatrixXd TMatrixToEigenMatrix(TMatrixD &matrix) {
  typedef Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> MapType;
  MapType m2map(matrix.GetMatrixArray(), matrix.GetNrows(), matrix.GetNcols());
  return m2map;
}
} // namespace myFuncs
