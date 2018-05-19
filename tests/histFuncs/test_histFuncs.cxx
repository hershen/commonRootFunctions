#include "catch.hpp"

#include "histFuncs.h"
#include <iostream>

#include "TH1D.h"

using namespace myFuncs;
using Catch::Equals;

std::vector<double> oneToEleven{1.0,
                                std::pow(11, 0.1),
                                std::pow(11, 0.2),
                                std::pow(11, 0.3),
                                std::pow(11, 0.4),
                                std::pow(11, 0.5),
                                std::pow(11, 0.6),
                                std::pow(11, 0.7),
                                std::pow(11, 0.8),
                                std::pow(11, 0.9),
                                11.0};

std::vector<double> NightyOneToHundredOne{91.,
                                          91 * std::pow(101. / 91., 0.1),
                                          91 * std::pow(101. / 91., 0.2),
                                          91 * std::pow(101. / 91., 0.3),
                                          91 * std::pow(101. / 91., 0.4),
                                          91 * std::pow(101. / 91., 0.5),
                                          91 * std::pow(101. / 91., 0.6),
                                          91 * std::pow(101. / 91., 0.7),
                                          91 * std::pow(101. / 91., 0.8),
                                          91 * std::pow(101. / 91., 0.9),
                                          101.0};

SCENARIO("Check makeConstWidthOnLogScale function", "[makeConstWidthOnLogScale]") {
  GIVEN("A histogram with x axis between 1 and 11 with 10 bins") {
    TH1D hist("hist", "", 10, 1, 11);
    WHEN("converting axis to const log width") {
      myFuncs::makeConstWidthOnLogScale(hist.GetXaxis());
      THEN("It should equal what we expect") {
        for (int i = 0; i <= hist.GetNbinsX(); ++i) {
          CHECK((*hist.GetXaxis()->GetXbins())[i] == Approx(oneToEleven[i]));
        }
      }
    }
  }
  GIVEN("A histogram with x axis between 91 and 101 with 10 bins") {
    TH1D hist("hist", "", 10, 91, 101);
    WHEN("converting axis to const log width") {
      myFuncs::makeConstWidthOnLogScale(hist.GetXaxis());
      THEN("It should equal what we expect") {
        for (int i = 0; i <= hist.GetNbinsX(); ++i) {
          CHECK((*hist.GetXaxis()->GetXbins())[i] == Approx(NightyOneToHundredOne[i]));
        }
      }
    }
  }
}
