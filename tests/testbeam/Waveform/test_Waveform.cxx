#include "catch.hpp"
#include <iostream>

#include "testbeam/Waveform.h"

using namespace myFuncs::testbeam;

TEST_CASE("Test Waveform class", "[Waveform]") {
  SECTION("Throw - negative dt") {
    std::vector<int> vector;
    CHECK_THROWS_AS(Waveform<int>(vector, -1.5), std::invalid_argument);
  }
  SECTION("Throw - zero dt") {
    std::vector<int> vector;
    CHECK_THROWS_AS(Waveform<int>(vector, 0), std::invalid_argument);
  }
}
TEST_CASE("Test empty Waveform", "[Waveform]") {

  std::vector<int> vector;
  Waveform<int> waveform(vector, 1.5);
  CHECK(waveform.getTimes() == std::vector<double>());
  CHECK(waveform.getSamples() == std::vector<int>());
  CHECK(waveform.getDt() == Approx(1.5));

  const auto MaximumIdx_value = waveform.getMaximumIdx_value();
  CHECK(MaximumIdx_value.first == 0);
  CHECK(MaximumIdx_value.second == 0);

  CHECK(waveform.getMean(0,0) == Approx(0));
  CHECK(waveform.getMean() == Approx(0));
  CHECK(waveform.getStd(0,0) == Approx(0));
  CHECK(waveform.getStd() == Approx(0));

  Waveform<double> waveformDouble = waveform.transformToDouble();
  CHECK(waveformDouble.getSamples() == std::vector<double>());
  CHECK(waveformDouble.getDt() == Approx(1.5));

}
TEST_CASE("Test non empty Waveform", "[Waveform]") {

  SECTION("getTimes") {
    std::vector<int> vector(5);
    Waveform<int> waveform(vector, 1.5);
    std::vector<double> expected{0, 1.5, 3, 4.5, 6};
    CHECK(waveform.getTimes() == expected);
  }

  SECTION("Waveforms with content waveform") {
    std::vector<int> vector{1, 2, 3, 4, 5};
    Waveform<int> waveform(vector, 1.5);
    CHECK(waveform.getSamples() == vector);
    CHECK(waveform.getDt() == Approx(1.5));

    const auto MaximumIdx_value = waveform.getMaximumIdx_value();
    CHECK(MaximumIdx_value.first == 4);
    CHECK(MaximumIdx_value.second == 5);

    CHECK(waveform.getMean(1,3) == Approx(3));
    CHECK(waveform.getMean() == Approx(3.0));
    CHECK(waveform.getStd(1,3) == Approx(0.816496580927726));
    CHECK(waveform.getStd() == Approx(1.4142135623731));

    Waveform<double> waveformDouble = waveform.transformToDouble();
    CHECK(waveformDouble.getSamples() == std::vector<double>{1, 2, 3, 4, 5});
    CHECK(waveformDouble.getDt() == Approx(1.5));
  }
}
