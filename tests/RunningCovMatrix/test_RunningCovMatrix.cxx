#include "catch.hpp"
#include <iostream>

#include "RunningCovMatrix.h"
#include "mathFuncs.h"

#include "TMatrixD.h"

TEST_CASE("Test runningStatistics function with invalid size", "[runningStatistics]") {
  CHECK_THROWS_AS(myFuncs::RunningCovMatrix(1), std::invalid_argument);
}

TEST_CASE("Test runningStatistics function", "[runningStatistics]") {
  myFuncs::RunningCovMatrix runCovMat(4);

  SECTION("empty matrix") {
    CHECK(runCovMat.getNumVariables() == 4);
    CHECK(runCovMat.getNumDataPoints() == 0);
    CHECK(runCovMat.getCovarianceMatrix().GetNrows() == 4);
    CHECK(runCovMat.getCovarianceMatrix().GetNcols() == 4);
  }

  SECTION("Adding row to an empty RunningCovMatrix") {
    std::vector<double> vec4{1, 2, 3, 4};

    runCovMat.addVector(vec4);
    CHECK(runCovMat.getNumVariables() == 4);
    CHECK(runCovMat.getNumDataPoints() == 1);

    CHECK(runCovMat.getCovarianceMatrix().GetNrows() == 4);
    CHECK(runCovMat.getCovarianceMatrix().GetNcols() == 4);

    SECTION("Add another same size row to an empty RunningCovMatrix") {
      runCovMat.addVector(vec4);
      runCovMat.addVector(vec4);
      runCovMat.addVector(vec4);
      runCovMat.addVector(myFuncs::addToVector(vec4, 1.0));

      CHECK(runCovMat.getNumVariables() == 4);
      CHECK(runCovMat.getNumDataPoints() == 5);

      CHECK(runCovMat.getCovarianceMatrix().GetNrows() == 4);
      CHECK(runCovMat.getCovarianceMatrix().GetNcols() == 4);

      for (int i = 0; i < runCovMat.getCovarianceMatrix().GetNrows(); ++i) {
        for (int j = 0; j < runCovMat.getCovarianceMatrix().GetNcols(); ++j) {
          if (i >= j) {
            CHECK(runCovMat.getCovarianceMatrix()(i, j) == Approx(0.16));
          }
        }
      }
    }
  }
}

TEST_CASE("Test runningStatistics function with varrying numbers", "[runningStatistics]") {
  myFuncs::RunningCovMatrix runCovMat(4);
  runCovMat.addVector(std::vector<double>{1, 2, 3, 4});
  runCovMat.addVector(std::vector<double>{2, 4, 3, 5});
  runCovMat.addVector(std::vector<double>{8.1, 3.9, 4.7, 0.1});
  runCovMat.addVector(std::vector<double>{5, 6, 7, 8});
  runCovMat.addVector(std::vector<double>{8, 7, 1, 6});
  runCovMat.addVector(std::vector<double>{1.2, 8, 2.1, -3.2});
  runCovMat.addVector(std::vector<double>{9, 5.5, -1.5, 3.7});

  // Expectations from Excel's COVARIANCE.P function
  CHECK(runCovMat.getCovarianceMatrix()(0, 0) == Approx(10.56857142857140));
  CHECK(runCovMat.getCovarianceMatrix()(1, 0) == Approx(1.19));
  CHECK(runCovMat.getCovarianceMatrix()(2, 0) == Approx(-2.21142857142857));
  CHECK(runCovMat.getCovarianceMatrix()(3, 0) == Approx(2.37571428571429));
  CHECK(runCovMat.getCovarianceMatrix()(1, 1) == Approx(3.59714285714286));
  CHECK(runCovMat.getCovarianceMatrix()(2, 1) == Approx(-0.925714285714286));
  CHECK(runCovMat.getCovarianceMatrix()(3, 1) == Approx(-1.36857142857143));
  CHECK(runCovMat.getCovarianceMatrix()(2, 2) == Approx(6.21959183673469));
  CHECK(runCovMat.getCovarianceMatrix()(3, 2) == Approx(1.7330612244898));
  CHECK(runCovMat.getCovarianceMatrix()(3, 3) == Approx(12.1963265306122));
}
