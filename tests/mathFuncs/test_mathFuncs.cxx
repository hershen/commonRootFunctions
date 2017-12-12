#include "catch.hpp"

#include "mathFuncs.h"

using Catch::Matchers::Equals;
using namespace Catch::Matchers;

TEST_CASE("Test roundKeepDigits function", "[roundKeepDigits]") {
  CHECK(myFuncs::roundKeepDigits(0.0, 3) == Approx(0.0).margin(1e-9));

  CHECK(myFuncs::roundKeepDigits(5.67890, -1) == Approx(5.67890).margin(1e-9));

  CHECK(myFuncs::roundKeepDigits(5.67890, 0) == Approx(6).margin(1e-9));
  CHECK(myFuncs::roundKeepDigits(5.67890, 1) == Approx(5.7).margin(1e-9));
  CHECK(myFuncs::roundKeepDigits(5.67890, 2) == Approx(5.68).margin(1e-9));
  CHECK(myFuncs::roundKeepDigits(5.67890, 3) == Approx(5.679).margin(1e-9));
  CHECK(myFuncs::roundKeepDigits(5.67890, 8) == Approx(5.67890).margin(1e-9));

  // return a number that has at max 2 significant digits.
  // If the 2 significant digits are xy, then if xy > 35, return a number with only x as significant digit
  // If xy <= 35, return a number with xy as significant digits.
  // y is always rounded using the smaller digits.
  // If error < 0 returns error
  // double round_35rule(const double error);
}
