#include "catch.hpp"
#include <iostream>

#include "MVectorTemplate.h"

#include "TF1.h"
#include "TFitResult.h"

#define cline std::cout << "line = " << __LINE__ << std::endl;

using namespace myFuncs;

SCENARIO("Constant vectors", "[MVectorTemplate]") {
  GIVEN("An constant vector") {
    std::vector<double> constant(1e3, 10.0);
    WHEN("Creating a template out of it") { CHECK_THROWS_AS(MVectorTemplate(constant, 1.0), std::invalid_argument); }
  }
  GIVEN("A constant vector with one different element") {
    const double dt = 2.0;
    std::vector<double> oneNonZero(1e3, 10.0);
    oneNonZero[500] = 300;
    WHEN("Creating a template out of it") {
      MVectorTemplate myTemplate(oneNonZero, dt);
      THEN("Expect throw when adding a constant vector") {
        CHECK_THROWS_AS(myTemplate.addVector(std::vector<double>(1e3, 10.0), 1.0, *myTemplate.getTF1()), std::invalid_argument);
      }
    }
  }
}

SCENARIO("Check template from single vector", "[MVectorTemplate]") {
  GIVEN("A constant vector with one different element") {
    const double dt = 2.0;
    std::vector<double> oneNonZero(1e3, 10.0);
    oneNonZero[500] = 300;
    WHEN("Creating a template out of it") {
      MVectorTemplate myTemplate(oneNonZero, dt);
      THEN("non zero element should be normalized") { CHECK(myTemplate.getTemplateValues()[500] == Approx(1)); }
      THEN("Num averaged functions should be 1") { CHECK(myTemplate.getNumAveragedFuncs() == 1); }
      THEN("TF1 is valid") { CHECK(myTemplate.getTF1()->IsValid() == true); }
      THEN("Peak index is correct") { CHECK(myTemplate.getPeakIdx() == 500); }
      THEN("x value of first entry is 0") { CHECK(myTemplate.getXvalueOfFirstTemplateEntry() == Approx(0)); }

      double minRange, maxRange;
      myTemplate.getTF1()->GetRange(minRange, maxRange);
      THEN("Function range begins at 0") { CHECK(minRange == Approx(0)); }
      const double expectedMaxRange = (oneNonZero.size() - 1) * dt;
      THEN("Function range ends at (size-1)*dt") { CHECK(maxRange == Approx(expectedMaxRange)); }

      AND_WHEN("Reseting template values") {
        myTemplate.resetTemplateValues();
        THEN("There are no template values") { CHECK(myTemplate.getTemplateValues().empty() == true); }
        THEN("Number of averaged functions == 0") { CHECK(myTemplate.getNumAveragedFuncs() == 0); }
        THEN("TF1 is invalid") { CHECK(myTemplate.getTF1()->IsValid() == false); }
      }
    }
  }
}

SCENARIO("Adding constant vector to a template", "[MVectorTemplate]") {
  GIVEN("A constant vector with one different element") {
    const double dt = 2.0;
    std::vector<double> oneNonZero(1e3, 10.0);
    oneNonZero[500] = 300;
    WHEN("Creating a template out of it") {
      MVectorTemplate myTemplate(oneNonZero, 1.0);
      AND_WHEN("Adding a consant vector to it") {
        std::vector<double> constant(1e3, 0.0);
        THEN("It throws") { CHECK_THROWS_AS(myTemplate.addVector(constant, 1.0, *myTemplate.getTF1()), std::invalid_argument); }
      }
    }
  }
}

SCENARIO("Check adding same shape vector", "[MVectorTemplate]") {
  GIVEN("A template from a constant vector with a few element different than the constant") {
    const double dt = 2.0;
    const double mainAmplitude = 300;
    std::vector<double> oneNonZero(1e3, 10.0);
    for (int i = 495; i <= 505; ++i) {
      oneNonZero[i] = mainAmplitude;
    }
    MVectorTemplate myTemplate(oneNonZero, dt);
    // myTemplate.setDebugLevel(16);
    WHEN("Adding a similar vector, offset in time and with a different pedestal and amplitude") {
      const double offsetPedestal = -5;
      std::vector<double> offsetVector(1e3, offsetPedestal);
      const int binsOffset = 5;
      const double offsetAmplitude = 100;
      for (int i = 495 + binsOffset; i <= 505 + binsOffset; ++i) {
        offsetVector[i] = offsetAmplitude;
      }
      THEN("Check that fit results and number of averaged functions are correct ") {
        for (int iAdd = 1; iAdd <= 50; ++iAdd) {
          auto fitResult = myTemplate.addVector(offsetVector, 10, *myTemplate.getTF1());
          CHECK(fitResult.Value(0) == Approx(offsetAmplitude - offsetPedestal));
          CHECK(fitResult.Value(1) == Approx(offsetPedestal));
          CHECK(fitResult.Value(2) == Approx(binsOffset * dt));
          CHECK(myTemplate.getNumAveragedFuncs() == iAdd + 1);
        }
        CHECK(myTemplate.getXvalueOfFirstTemplateEntry() == Approx(0));
        CHECK(myTemplate.getTemplateSize() == oneNonZero.size() - binsOffset);
      }
    }

    AND_WHEN("Adding a similar vector, offset in time to the other side") {
      const double offsetPedestal = -5;
      std::vector<double> offsetVector(1e3, offsetPedestal);
      const int binsOffset = -5;
      const double offsetAmplitude = 100;
      for (int i = 495 + binsOffset; i <= 505 + binsOffset; ++i) {
        offsetVector[i] = offsetAmplitude;
      }
      THEN("Check that fit results and number of averaged functions are correct ") {
        for (int iAdd = 1; iAdd <= 50; ++iAdd) {
          auto fitResult = myTemplate.addVector(offsetVector, 10, *myTemplate.getTF1());
          CHECK(fitResult.Value(0) == Approx(offsetAmplitude - offsetPedestal));
          CHECK(fitResult.Value(1) == Approx(offsetPedestal));
          CHECK(fitResult.Value(2) == Approx(binsOffset * dt));
          CHECK(myTemplate.getNumAveragedFuncs() == iAdd + 1);
        }
        CHECK(myTemplate.getXvalueOfFirstTemplateEntry() == Approx(-binsOffset * dt));
        CHECK(myTemplate.getTemplateSize() == oneNonZero.size() - std::abs(binsOffset));
      }
    }
  }
}
