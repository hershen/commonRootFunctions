#include "catch.hpp"

#include "mathFuncs.h"

using Catch::Matchers::Equals;
using namespace Catch::Matchers;

TEST_CASE("Test roundKeepDigits function", "[roundKeepDigits]") {
  CHECK(myFuncs::roundKeepDigits(0.0, 3) == Approx(0.0));

  CHECK(myFuncs::roundKeepDigits(5.67890, -1) == Approx(5.67890));
  CHECK(myFuncs::roundKeepDigits(5.67890, 0) == Approx(5.67890));

  CHECK(myFuncs::roundKeepDigits(5.67890, 1) == Approx(6));
  CHECK(myFuncs::roundKeepDigits(5.67890, 2) == Approx(5.7));
  CHECK(myFuncs::roundKeepDigits(5.67890, 3) == Approx(5.68));
  CHECK(myFuncs::roundKeepDigits(5.67890, 4) == Approx(5.679));
  CHECK(myFuncs::roundKeepDigits(5.67890, 5) == Approx(5.67890));

  CHECK(myFuncs::roundKeepDigits(1.111, 1) == Approx(1));

  CHECK(myFuncs::roundKeepDigits(567890, 1) == Approx(600000));
  CHECK(myFuncs::roundKeepDigits(567890, 2) == Approx(570000));
  CHECK(myFuncs::roundKeepDigits(567890, 3) == Approx(568000));
  CHECK(myFuncs::roundKeepDigits(567890, 4) == Approx(567900));
}

TEST_CASE("Test round_35rule function", "[round_35rule]") {
  CHECK(myFuncs::round_35rule(-0.01) == Approx(-0.01));

  CHECK(myFuncs::round_35rule(11000) == Approx(11000));
  CHECK(myFuncs::round_35rule(78912) == Approx(80000));
  CHECK(myFuncs::round_35rule(1.1234) == Approx(1.1));
  CHECK(myFuncs::round_35rule(3.1789) == Approx(3.2));
  CHECK(myFuncs::round_35rule(3.789) == Approx(4));
  CHECK(myFuncs::round_35rule(123.456789) == Approx(120));
  CHECK(myFuncs::round_35rule(123.1234) == Approx(120));
}

TEST_CASE("Test roundAccordingToError function", "[roundAccordingToError]") {
  // Assumes that error has only 2 significant digits
  CHECK(myFuncs::roundAccordingToError(1234, 0.0) == Approx(1234));
  CHECK(myFuncs::roundAccordingToError(-0.0543, 0.0) == Approx(-0.0543));

  CHECK(myFuncs::roundAccordingToError(5.67890, 1) == Approx(6));
  CHECK(myFuncs::roundAccordingToError(5.67890, 1.1) == Approx(5.7));

  CHECK(myFuncs::roundAccordingToError(5.67890, 11000) == Approx(0.0));
  CHECK(myFuncs::roundAccordingToError(67812.1234, 99000) == Approx(68000));
  CHECK(myFuncs::roundAccordingToError(666, 11000) == Approx(1000));
}

TEST_CASE("Test sampleMean function", "[sampleMean]") {

  SECTION("Empty vector") {
    std::vector<double> emptyVector;
    CHECK(myFuncs::sampleMean(emptyVector) == Approx(0));
    CHECK(myFuncs::sampleMean(emptyVector.begin(), emptyVector.end()) == Approx(0));
  }

  SECTION("Vector double, integer result") {
    std::vector<double> vector{1, 2, 3, 4, 5};
    CHECK(myFuncs::sampleMean(vector) == Approx(3));
    CHECK(myFuncs::sampleMean(vector.begin(), vector.end()) == Approx(3));
  }

  SECTION("Vector int, integer result") {
    std::vector<int> vector{1, 2, 3, 4, 5};
    CHECK(myFuncs::sampleMean(vector) == Approx(3));
    CHECK(myFuncs::sampleMean(vector.begin(), vector.end()) == Approx(3));
  }

  SECTION("Vector double, non integer result") {
    std::vector<double> vector{1, 2, 3, 4, 5, 6};
    CHECK(myFuncs::sampleMean(vector) == Approx(3.5));
    CHECK(myFuncs::sampleMean(vector.begin(), vector.end()) == Approx(3.5));
  }

  SECTION("Vector int, non integer result") {
    std::vector<int> vector{1, 2, 3, 4, 5, 6};
    CHECK(myFuncs::sampleMean(vector) == Approx(3.5));
    CHECK(myFuncs::sampleMean(vector.begin(), vector.end()) == Approx(3.5));
  }
}

TEST_CASE("Test linearInterpolate function", "[linearInterpolate]") {

  SECTION("Zero denominator") { CHECK_THROWS(myFuncs::linearInterpolate(1.0, 1.0, 2.0, 3.0, 1.5)); }
  SECTION("Double x, double y") { CHECK(myFuncs::linearInterpolate(1.0, 2.0, 2.0, 4.0, 1.5) == Approx(3)); }
  SECTION("Double x, int y") { CHECK(myFuncs::linearInterpolate(1.0, 2.0, 2, 4, 1.5) == Approx(3)); }
  SECTION("int x, double y") { CHECK(myFuncs::linearInterpolate(1, 2, 2.0, 4.0, 1.5) == Approx(3)); }
  SECTION("int x, int y") { CHECK(myFuncs::linearInterpolate(1, 2, 2, 4, 1.5) == Approx(3)); }
}
