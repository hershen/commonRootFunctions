#include "catch.hpp"
#include <iostream>

#include "RunningStatistics.h"

#include "TRandom3.h"

TEST_CASE("Test runningStatistics function", "[runningStatistics]") {

  SECTION("Empty statistics") {
    myFuncs::RunningStatistics runStat;
    CHECK(runStat.getNumElements() == 0);
    CHECK(runStat.getSampleMean() == 0);
    CHECK(runStat.getSampleStd() == 0);
    CHECK(runStat.getSampleVariance() == 0);
  }

  SECTION("Test getNumElements") {
    myFuncs::RunningStatistics runStat;
    CHECK(runStat.getNumElements() == 0);
    runStat.addElement(1);
    CHECK(runStat.getNumElements() == 1);
    runStat.addElement(-1);
    CHECK(runStat.getNumElements() == 2);
    runStat.addElement(1.1);
    CHECK(runStat.getNumElements() == 3);
  }

  SECTION("Test running function on small number of elements") {
    myFuncs::RunningStatistics runStat;

    runStat.addElement(1);
    CHECK(runStat.getSampleMean() == Approx(1));
    CHECK(runStat.getSampleStd() == Approx(0));
    CHECK(runStat.getSampleVariance() == Approx(0));
    runStat.addElement(2.0);
    CHECK(runStat.getSampleMean() == Approx(1.5));
    CHECK(runStat.getSampleStd() == Approx(0.707106781186548));
    CHECK(runStat.getSampleVariance() == Approx(0.5));
    runStat.addElement(3);
    CHECK(runStat.getSampleMean() == Approx(2));
    CHECK(runStat.getSampleStd() == Approx(1));
    CHECK(runStat.getSampleVariance() == Approx(1));
    runStat.addElement(4.0);
    CHECK(runStat.getSampleMean() == Approx(2.5));
    CHECK(runStat.getSampleStd() == Approx(1.29099444873581));
    CHECK(runStat.getSampleVariance() == Approx(1.66666666666667));
    runStat.addElement(5);
    CHECK(runStat.getSampleMean() == Approx(3));
    CHECK(runStat.getSampleStd() == Approx(1.58113883008419));
    CHECK(runStat.getSampleVariance() == Approx(2.5));
  }

  SECTION("Test running function on large sample from uniform distribution. These tests should pass ~99% of the times.") {
    myFuncs::RunningStatistics runStat;
    TRandom3 r(0);
    const long numElements = 1e4;
    const double offset = 1e9;
    for (long i = 0; i < numElements; ++i) {
      runStat.addElement(r.Uniform() + offset);
    }

    CHECK(runStat.getSampleMean() == Approx(offset + 0.5).margin(1. / std::sqrt(12) / sqrt(numElements) * 3.));

    // Standard deviation of variance
    // (https://stats.stackexchange.com/questions/105337/asymptotic-distribution-of-sample-variance-of-non-normal-sample/105338#105338)
    const double mu4 = 1. / 80;
    const double stdOfVariance = std::sqrt(1. / numElements * (mu4 - (numElements - 3) / (numElements - 1) / 144));
    CHECK(runStat.getSampleVariance() == Approx(1. / 12).margin(stdOfVariance * 3));
  }
}
