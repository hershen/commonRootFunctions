#pragma once

// std
#include <iostream>

// Eigen
#include <Eigen/Dense>

// Root
#include "TFile.h"
#include "TMatrixD.h"
#include "TMatrixDSym.h"
#include "TTree.h"

namespace myFuncs {

// Returns the contents of a TMatrix in an Eigen matrix.
Eigen::MatrixXd TMatrixToEigenMatrix(TMatrixD &matrix) {
  typedef Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> MapType;
  MapType m2map(matrix.GetMatrixArray(), matrix.GetNrows(), matrix.GetNcols());
  return m2map;
}

// Save Eigen matrix into TTree
template <typename EigenType>
void saveToFile(EigenType &matrix, const std::string &filename, const std::string &treeName = "tree",
                const std::string &options = "RECREATE") {
  TFile f(filename.c_str(), options.c_str());
  TTree tree(treeName.c_str(), "");
  const size_t numElements = matrix.rows() * matrix.cols();

  double *matrixData = matrix.data();
  int rows = matrix.rows();
  int cols = matrix.cols();
  tree.Branch("matrixArray", matrixData, ("matrixArray[" + std::to_string(numElements) + "]/D").c_str());
  tree.Branch("matrixRows", &rows);
  tree.Branch("matrixCols", &cols);
  tree.Fill();
  f.Write();
  f.Close();
}

// Load Eigen matrix from TTree
Eigen::MatrixXd loadFromFile(const std::string &filename, const std::string &treeName = "tree") {
  TFile f(filename.c_str(), "READ");
  TTree *tree = (TTree *)f.Get(treeName.c_str());
  int rows;
  int cols;
  tree->SetBranchAddress("matrixRows", &rows);
  tree->SetBranchAddress("matrixCols", &cols);

  // Read only cols and rows branches
  tree->SetBranchStatus("*", 0);
  tree->SetBranchStatus("matrixRows", 1);
  tree->SetBranchStatus("matrixCols", 1);
  tree->GetEntry(0);

  // Prealocate vector size so that data can be read into it.
  const int numElements = rows * cols;
  std::vector<double> matrixArray(numElements);

  // Read matrix array
  tree->SetBranchStatus("matrixArray", 1);
  tree->SetBranchAddress("matrixArray", matrixArray.data());

  tree->GetEntry(0);
  delete tree;
  f.Close();

  typedef Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>> MapType;
  MapType m2map(matrixArray.data(), rows, cols);
  return m2map;
}

template <class EigenType>
TMatrixD EigenToTMatrix(const EigenType &eigenMatrix) {
  return TMatrixD(eigenMatrix.rows(), eigenMatrix.cols(), eigenMatrix.data(), "F"); //"F" is for column wise
}
} // namespace myFuncs
