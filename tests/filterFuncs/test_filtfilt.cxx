#include "catch.hpp"

#include "filterFuncs.h"
#include "mathFuncs.h"

using Catch::Matchers::Equals;
using namespace Catch::Matchers;
using myFuncs::DSP::filtfilt;

const std::vector<double> inputs{1, 2, 3, 4, 5, 6, 5, 4, 3, 2, 1};

TEST_CASE("Test filtfilt function", "[filtfilt]") {

  SECTION("Empty denom") {
    const std::vector<double> nom;
    const std::vector<double> denom;

    CHECK_THAT(filtfilt(nom, denom, inputs), Equals(std::vector<double>()));
  }

  SECTION("Empty nom") {
    const std::vector<double> nom;
    const std::vector<double> denom{1};

    CHECK_THAT(filtfilt(nom, denom, inputs), Equals(std::vector<double>()));
  }

  // SECTION("Identity") {
  //   const std::vector<double> nom{1};
  //   const std::vector<double> denom{1};
  //
  //   CHECK_THAT(filtfilt(nom, denom, inputs), Equals(inputs));
  // }
  //
  // SECTION("Times 2") {
  //   const std::vector<double> nom{2};
  //   const std::vector<double> denom{1};
  //
  //   CHECK_THAT(filtfilt(nom, denom, inputs), Equals(myFuncs::scaleVector(inputs, 2.0)));
  // }
  //
  // SECTION("Times 0.5") {
  //   const std::vector<double> nom{1};
  //   const std::vector<double> denom{2};
  //
  //   CHECK_THAT(filtfilt(nom, denom, inputs), Equals(myFuncs::scaleVector(inputs, 0.5)));
  // }
  //
  // SECTION("Test1") {
  //   const std::vector<double> nom{2, 4};
  //   const std::vector<double> denom{5, 7};
  //   const std::vector<double> expected{0.400000000000000, 1.040000000000000, 1.344000000000000, 2.118400000000000,
  //                                      2.234240000000000, 3.272064000000001, 2.219110400000000, 2.493245440000000,
  //                                      0.909456384000000, 1.926761062400001, -0.697465487360000};
  //
  //   const auto filtered = filtfilt(nom, denom, inputs);
  //   REQUIRE(filtered.size() == expected.size());
  //   for (uint i = 0; i < inputs.size(); ++i)
  //     REQUIRE(std::abs(filtered[i] - expected[i]) < 1e-9);
  // }
  //
  // SECTION("Test2") {
  //   const std::vector<double> nom{0.5, 3.5, 5.8, 9.9, 10.33};
  //   const std::vector<double> denom{9.9, 10.33, 2.1, 3.33, 4.7};
  //   const std::vector<double> expected{0.050505050505051, 0.401846750331599, 1.014430551819594, 2.273623475335901,
  //                                      3.720973931816268, 4.604141062829779, 5.593554323003392, 5.514445867111669,
  //                                      5.456458878350476, 3.835196421697355, 2.062732046841669};
  //
  //   const auto filtered = filtfilt(nom, denom, inputs);
  //   REQUIRE(filtered.size() == expected.size());
  //   for (uint i = 0; i < inputs.size(); ++i)
  //     REQUIRE(std::abs(filtered[i] - expected[i]) < 1e-9);
  // }
  //
  // SECTION("Test3") {
  //   const std::vector<double> nom{0.5, 3.5, 5.8, 9.9, 10.33, 2, 3.1, 8.3, 9.7, 11.3, 1.7, 3, 2.1};
  //   const std::vector<double> denom{9.9, 10.33, 2.1, 3.33, 4.7, 5.1, 1.5, 2.5, 5.2, 3.6, 6.3, 1.3, 3.1};
  //   const std::vector<double> expected{0.050505050505051, 0.401846750331599, 1.014430551819594, 2.273623475335901,
  //                                      3.720973931816268, 4.780143511559500, 5.912414797970388, 6.618884007032738,
  //                                      7.128286168503705, 6.797498157036226, 5.307988818810417};
  //
  //   const auto filtered = filtfilt(nom, denom, inputs);
  //   REQUIRE(filtered.size() == expected.size());
  //   for (uint i = 0; i < inputs.size(); ++i)
  //     REQUIRE(std::abs(filtered[i] - expected[i]) < 1e-9);
  // }
  //
  // SECTION("CR-RC4") {
  //   // Create step function as input
  //   const int zeros = 300;
  //   std::vector<double> stepFunc(zeros, 0.0);
  //   stepFunc.insert(stepFunc.end(), zeros, 1.0);
  //
  //   std::vector<double> nom;
  //   std::vector<double> denom;
  //   std::tie(nom, denom) = myFuncs::DSP::getCR_RCnCoefficients(4, 500, 1.0 / 2.0);
  //   const auto filtered = filtfilt(nom, denom, stepFunc);
  //
  //   REQUIRE(filtered.size() == stepFunc.size());
  //
  //   // AnalyticalFunction
  //   auto analyticalFunc = myFuncs::analytical_RC_CRn(4, 500., 1.0, zeros * 2.0 - 2.0 / 2.0, stepFunc.size() * 2.0);
  //
  //   // Check filter
  //   for (uint i = 0; i < filtered.size(); ++i)
  //     REQUIRE(std::abs(filtered[i] - analyticalFunc.Eval(i * 2.0)) <= 1e-7); //Pretty low tollerance!
  //
  // }
}
