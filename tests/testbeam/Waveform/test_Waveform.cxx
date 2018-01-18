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
  CHECK(waveform.getSamples() == std::vector<int>());
  CHECK(waveform.getDt() == Approx(1.5));

  const auto MaximumIdx_value = waveform.getMaximumIdx_value();
  CHECK(MaximumIdx_value.first == 0);
  CHECK(MaximumIdx_value.second == 0);

  CHECK(waveform.getMean(0, 0) == Approx(0));
  CHECK(waveform.getMean() == Approx(0));
  CHECK(waveform.getStd(0, 0) == Approx(0));
  CHECK(waveform.getStd() == Approx(0));

  // Transform to double
  Waveform<double> waveformDouble = waveform.transformToDouble();
  CHECK(waveformDouble.getSamples() == std::vector<double>());
  CHECK(waveformDouble.getDt() == Approx(1.5));

  // Moving average
  const auto averageEach_nWaveform = waveform.averageEach_n(2);
  CHECK(averageEach_nWaveform.getSamples() == std::vector<double>());
  CHECK(averageEach_nWaveform.getDt() == Approx(3.0));
}
TEST_CASE("Test non empty Waveform", "[Waveform]") {

  SECTION("Waveforms with content waveform") {
    std::vector<int> vector{1, 2, 3, 4, 5};
    Waveform<int> waveform(vector, 1.5);
    CHECK(waveform.getSamples() == vector);
    CHECK(waveform.getDt() == Approx(1.5));

    const auto MaximumIdx_value = waveform.getMaximumIdx_value();
    CHECK(MaximumIdx_value.first == 4);
    CHECK(MaximumIdx_value.second == 5);

    CHECK(waveform.getMean(1, 3) == Approx(3));
    CHECK(waveform.getMean() == Approx(3.0));
    CHECK(waveform.getStd(1, 3) == Approx(0.816496580927726));
    CHECK(waveform.getStd() == Approx(1.4142135623731));

    Waveform<double> waveformDouble = waveform.transformToDouble();
    CHECK(waveformDouble.getSamples() == std::vector<double>{1, 2, 3, 4, 5});
    CHECK(waveformDouble.getDt() == Approx(1.5));

    SECTION("average every n") {
      const auto averageEach_nWaveform = waveform.averageEach_n(2);
      CHECK(averageEach_nWaveform.getSamples() == std::vector<double>{1.5, 3.5});
      CHECK(averageEach_nWaveform.getDt() == Approx(3.0));

      // Moving average
      const auto averageEach_nWaveform1 = averageEach_nWaveform.averageEach_n(1);
      CHECK(averageEach_nWaveform1.getSamples() == std::vector<double>{1.5, 3.5});
      CHECK(averageEach_nWaveform1.getDt() == Approx(3.0));

      // Moving average
      const auto averageEach_nWaveform2 = averageEach_nWaveform1.averageEach_n(2);
      CHECK(averageEach_nWaveform2.getSamples() == std::vector<double>{2.5});
      CHECK(averageEach_nWaveform2.getDt() == Approx(6.0));
    }

    SECTION("time shift") {
      CHECK(waveform.timeShift(0.5 * 1.5).getSamples() == std::vector<double>{1, 1.5, 2.5, 3.5, 4.5});
      CHECK(waveform.timeShift(-0.5 * 1.5).getSamples() == std::vector<double>{1.5, 2.5, 3.5, 4.5, 5});
      CHECK(waveform.timeShift(1.5).getSamples() == std::vector<double>{1, 1, 2, 3, 4});
      CHECK(waveform.timeShift(-1.5).getSamples() == std::vector<double>{2, 3, 4, 5, 5});

      //Shift by +- 3 bins
      CHECK(waveform.timeShift(1.5*3.0).getSamples() == std::vector<double>{1, 1, 1, 1, 2});
      CHECK(waveform.timeShift(-1.5*3.0).getSamples() == std::vector<double>{4, 5, 5, 5, 5});
    }
  }
}

TEST_CASE("Test operator+", "[Waveform]") {
  SECTION("Different dts"){
    Waveform<int> dt1(std::vector<int>(), 1.0);
    Waveform<double> dt2(std::vector<double>(), 2.0);

    CHECK_THROWS_AS(dt1+dt2, std::invalid_argument);
  }
  SECTION("Empty waveforms"){
    Waveform<int> emptywaveformInt(std::vector<int>(), 1.0);
    Waveform<double> emptywaveformDouble(std::vector<double>(), 1.0);

    CHECK((emptywaveformInt+emptywaveformDouble).getSamples() == std::vector<double>());
  }
  SECTION("1 empty 1 not"){
    Waveform<int> emptywaveformInt(std::vector<int>(), 1.0);
    Waveform<double> emptywaveformDouble(std::vector<double>{1,2,3,4}, 1.0);

    CHECK_THROWS_AS((emptywaveformInt+emptywaveformDouble), std::invalid_argument);
  }
  SECTION("both non empty"){
    Waveform<int> emptywaveformInt(std::vector<int>{1,2,3,4}, 1.0);
    Waveform<double> emptywaveformDouble(std::vector<double>{1.5,2.5,3.5,4.5}, 1.0);

    CHECK((emptywaveformInt+emptywaveformDouble).getSamples() == std::vector<double>{2.5,4.5,6.5,8.5});
  }

}
