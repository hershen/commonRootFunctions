#include "catch.hpp"
#include <iostream>

#include "testbeam/testbeamFuncs.h"

using namespace myFuncs::testbeam;

SCENARIO("Test getTimes class", "[getTimes]") {
  WHEN("Calling the function with size 0 and dt 0") {
    THEN("Return empty vector") { CHECK(getTimes(0, 0) == std::vector<double>()); }
  }
  AND_WHEN("Calling the function with size 10 and dt 0") {
    THEN("Return a vector of zeros") { CHECK(getTimes(10, 0) == std::vector<double>(10, 0)); }
  }
  AND_WHEN("Calling the function with size 10 and dt 2") {
    THEN("Return vector with 2 increments") {
      std::vector<double> expected{0, 2, 4, 6, 8, 10, 12, 14, 16, 18};
      CHECK(getTimes(10, 2) == expected);
    }
  }
  AND_WHEN("Calling the function with size 10 and dt -2") {
    THEN("Return vector with -2 increments") {
      std::vector<double> expected{0, -2, -4, -6, -8, -10, -12, -14, -16, -18};
      CHECK(getTimes(10, -2) == expected);
    }
  }
}

SCENARIO("Test reductionFactorToEntriesToChop", "reductionFactorToEntriesToChop") {
  CHECK(reductionFactorToEntriesToChop(1) == 2800);
  CHECK(reductionFactorToEntriesToChop(20) == 140);
  CHECK(reductionFactorToEntriesToChop(25) == 112);
  CHECK(reductionFactorToEntriesToChop(100) == 28);
  CHECK(reductionFactorToEntriesToChop(-1) == 0);
}
