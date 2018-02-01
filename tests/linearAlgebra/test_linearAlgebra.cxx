#include "catch.hpp"

#include "linearAlgebraFuncs.h"

#include "TMatrixD.h"
#include "TRandom3.h"

using namespace myFuncs;
//
// // The matcher class
// template <typename T>
// struct MatrixEqualsMatcher : public Catch::MatcherBase<T> {
//
//   MatrixEqualsMatcher(T const &comparator, const double epsilon = std::numeric_limits<float>::epsilon() * 100)
//       : m_comparator(comparator), m_epsilon(epsilon) {}
//
//   virtual bool match(T const &matrix) const {
// if(matrix.GetNcolsreturn true; }
//   virtual std::string describe() const { return "Equals: i equal"; }
//   const T m_comparator;
//   const double m_epsilon;
// };
//
// template <typename T>
// MatrixEqualsMatcher<T> Equals(const T &comparator) {
//   return MatrixEqualsMatcher<T>(comparator);
// }
// // ...
//
// // Usage
// TEST_CASE("Integers are within a range") {
//   CHECK_THAT(TMatrixD(1,1), Equals(TMatrixD(1,1)));
//   CHECK_THAT(TMatrixD(10,10), Equals(TMatrixD(1,10)));
// }
//
// // TEST_CASE("Test symmertizeLowerLeftTriangle function", "[symmertizeLowerLeftTriangle]") {
// //   TRandom3 r(0);
// //   for (int i = 1; i < 20; ++i) {
// //     SECTION("Matrix with size " + std::to_string(i)) {
// //       TMatrixD matrix;
// //       symmertizeLowerLeftTriangle(matrix);
// //       CHECK(1 == 2);
// //     }
// //   }
// // }
//

template <class EigenMatrix>
void isEqual(const TMatrixD &tMatrix, const EigenMatrix &eigenMatrix) {
  CHECK(tMatrix.GetNrows() == eigenMatrix.rows());
  CHECK(tMatrix.GetNcols() == eigenMatrix.cols());
  for (int row = 0; row < tMatrix.GetNrows(); ++row) {
    for (int col = 0; col < tMatrix.GetNcols(); ++col) {
      CHECK(tMatrix(row, col) == Approx(eigenMatrix(row, col)));
    }
  }
}

// template <class EigenMatrix>
// bool isLowerEqual(const TMatrixD &tMatrix, const EigenMatrix &eigenMatrix) {
//   CHECK(tMatrix.GetNrows() == eigenMatrix.rows());
//   CHECK(tMatrix.GetNcols() == eigenMatrix.cols());
//   for (int row = 0; row < tMatrix.GetNrows(); ++row) {
//     for (int col = row; col < tMatrix.GetNcols(); ++col) {
//       if (tMatrix(row, col) == Approx(eigenMatrix(row, col))) {
//         return false;
//       }
//     }
//   }
//   return true;
// }

TEST_CASE("Test TMatrixToEigenMatrix function", "[TMatrixToEigenMatrix]") {
  SECTION("Check random matrices are equal") {
    TRandom3 r(0);
    for (int rows = 1; rows < 10; ++rows) {
      for (int cols = 1; cols < 10; ++cols) {
        SECTION("Matrix with size " + std::to_string(rows) + "x" + std::to_string(cols)) {
          TMatrixD tMatrix(rows, cols);
          double seed = r.Uniform(100);
          tMatrix.Randomize(1, 1e4, seed);
          const auto eigenMatrix = TMatrixToEigenMatrix(tMatrix);
          isEqual(tMatrix, eigenMatrix);
        }
      }
    }
  }
}

SCENARIO("Test saveToFile/loadFromFile function", "[saveToFile][loadFromFile]") {
  GIVEN("A matrix") {
    Eigen::MatrixXd mat = Eigen::MatrixXd::Random(2, 3);
    WHEN("trying to save") {
      const std::string filename = "/home/hershen/temp/testingEigen.root";
      myFuncs::saveToFile(mat, filename);
      THEN("Expect the matrix") {
        Eigen::MatrixXd loadedMatrix = myFuncs::loadFromFile(filename);
        CHECK(loadedMatrix.rows() == mat.rows());
        CHECK(loadedMatrix.cols() == mat.cols());
        for (int row = 0; row < mat.rows(); ++row) {
          for (int col = 0; col < mat.cols(); ++col) {
            CHECK(mat(row, col) == Approx(loadedMatrix(row, col)));
          }
        }
      }
    }

    WHEN("Saving a different matrix in a file with a specified tree name") {
      const std::string filename = "/home/hershen/temp/testingEigen2.root";
      mat = mat * 10;
      myFuncs::saveToFile(mat, filename, "myTree");
      THEN("Expect the matrix") {
        Eigen::MatrixXd loadedMatrix = myFuncs::loadFromFile(filename, "myTree");
        CHECK(loadedMatrix.rows() == mat.rows());
        CHECK(loadedMatrix.cols() == mat.cols());
        for (int row = 0; row < mat.rows(); ++row) {
          for (int col = 0; col < mat.cols(); ++col) {
            CHECK(mat(row, col) == Approx(loadedMatrix(row, col)));
          }
        }
      }
    }
  }
}

SCENARIO("Convert Eigen matrix to TMatrixD", "[EigenToTMatrix]") {
  GIVEN("A an Eigen Matrix") {
    Eigen::MatrixXd *eigenMatrix = new Eigen::MatrixXd(Eigen::MatrixXd::Random(8, 11));
    WHEN("Converting it to a TMatrixD") {
      TMatrixD tMatrix = myFuncs::EigenToTMatrix(*eigenMatrix);
      THEN("Matrices should be equal") { isEqual(tMatrix, *eigenMatrix); }
    }
  }
}
