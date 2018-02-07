#include "catch.hpp"
#include <iostream>

#include "testbeam/Waveform.h"

using namespace myFuncs::testbeam;

TEST_CASE("Test Waveform class", "[Waveform]") {
  SECTION("Throw - negative dt") {
    std::vector<int> vector;
    CHECK_THROWS_AS(Waveform(vector, -1.5), std::invalid_argument);
  }
  SECTION("Throw - zero dt") {
    std::vector<int> vector;
    CHECK_THROWS_AS(Waveform(vector, 0), std::invalid_argument);
  }
}
TEST_CASE("Test empty Waveform", "[Waveform]") {

  std::vector<int> vector;
  Waveform waveform(vector, 1.5);
  CHECK(waveform.getSamples() == std::vector<double>());
  CHECK(waveform.getDt() == Approx(1.5));

  const auto MaximumIdx_value = waveform.getMaximumIdx_value();
  CHECK(MaximumIdx_value.first == 0);
  CHECK(MaximumIdx_value.second == 0);

  CHECK(waveform.getMean(0, 0) == Approx(0));
  CHECK(waveform.getMean() == Approx(0));
  CHECK(waveform.getStd(0, 0) == Approx(0));
  CHECK(waveform.getStd() == Approx(0));

  // Moving average
  auto averageEach_nWaveform = waveform;
  averageEach_nWaveform.averageEach_n(2);
  CHECK(averageEach_nWaveform.getSamples() == std::vector<double>());
  CHECK(averageEach_nWaveform.getDt() == Approx(3.0));
}
TEST_CASE("Test non empty Waveform", "[Waveform]") {

  SECTION("Waveforms with content waveform") {
    std::vector<int> vectorInt{1, 2, 3, 4, 5};
    std::vector<double> vectorDouble{1, 2, 3, 4, 5};
    Waveform waveform(vectorInt, 1.5);
    CHECK(waveform.getSamples() == vectorDouble);
    CHECK(waveform.getDt() == Approx(1.5));

    const auto MaximumIdx_value = waveform.getMaximumIdx_value();
    CHECK(MaximumIdx_value.first == 4);
    CHECK(MaximumIdx_value.second == 5);

    CHECK(waveform.getMean(1, 3) == Approx(3));
    CHECK(waveform.getMean() == Approx(3.0));
    CHECK(waveform.getStd(1, 3) == Approx(0.816496580927726));
    CHECK(waveform.getStd() == Approx(1.4142135623731));

    CHECK(waveform.getMaxMinusPedestal(1) == Approx(4));
    CHECK(waveform.getMaxMinusPedestal(2) == Approx(3.5));
    CHECK(waveform.getMaxMinusPedestal(5) == Approx(2));

    SECTION("average every n") {
      auto averageEach_nWaveform = waveform;
      averageEach_nWaveform.averageEach_n(2);
      CHECK(averageEach_nWaveform.getSamples() == std::vector<double>{1.5, 3.5});
      CHECK(averageEach_nWaveform.getDt() == Approx(3.0));

      // Moving average
      averageEach_nWaveform.averageEach_n(1);
      CHECK(averageEach_nWaveform.getSamples() == std::vector<double>{1.5, 3.5});
      CHECK(averageEach_nWaveform.getDt() == Approx(3.0));

      // Moving average
      auto averageEach_nWaveform2 = averageEach_nWaveform;
      averageEach_nWaveform2.averageEach_n(2);
      CHECK(averageEach_nWaveform2.getSamples() == std::vector<double>{2.5});
      CHECK(averageEach_nWaveform2.getDt() == Approx(6.0));
    }
  }
}

TEST_CASE("timeShift", "[Waveform]") {
  std::vector<double> vectorDouble{1, 2, 3, 4, 5};
  Waveform waveform(vectorDouble, 1.5);
  SECTION("time shift + 0.5 bins") {
    waveform.timeShift(0.5 * 1.5);
    CHECK(waveform.getSamples() == std::vector<double>{1, 1.5, 2.5, 3.5, 4.5});
  }
  SECTION("time shift - 0.5 bins") {
    waveform.timeShift(-0.5 * 1.5);
    CHECK(waveform.getSamples() == std::vector<double>{1.5, 2.5, 3.5, 4.5, 5});
  }
  SECTION("time shift +1 bins") {
    waveform.timeShift(1.5);
    CHECK(waveform.getSamples() == std::vector<double>{1, 1, 2, 3, 4});
  }
  SECTION("time shift -1 bins") {
    waveform.timeShift(-1.5);
    CHECK(waveform.getSamples() == std::vector<double>{2, 3, 4, 5, 5});
  }
  SECTION("time shift +3 bins") {
    waveform.timeShift(1.5 * 3.0);
    CHECK(waveform.getSamples() == std::vector<double>{1, 1, 1, 1, 2});
  }
  SECTION("time shift -3 bins") {
    waveform.timeShift(-1.5 * 3.0);
    CHECK(waveform.getSamples() == std::vector<double>{4, 5, 5, 5, 5});
  }
}

// TEST_CASE("Test operator+", "[Waveform]") {
//   SECTION("Different dts"){
//     Waveform<int> dt1(std::vector<int>(), 1.0);
//     Waveform<double> dt2(std::vector<double>(), 2.0);
//
//     CHECK_THROWS_AS(dt1+dt2, std::invalid_argument);
//   }
//   SECTION("Empty waveforms"){
//     Waveform<int> emptywaveformInt(std::vector<int>(), 1.0);
//     Waveform<double> emptywaveformDouble(std::vector<double>(), 1.0);
//
//     CHECK((emptywaveformInt+emptywaveformDouble).getSamples() == std::vector<double>());
//   }
//   SECTION("1 empty 1 not"){
//     Waveform<int> emptywaveformInt(std::vector<int>(), 1.0);
//     Waveform<double> emptywaveformDouble(std::vector<double>{1,2,3,4}, 1.0);
//
//     CHECK_THROWS_AS((emptywaveformInt+emptywaveformDouble), std::invalid_argument);
//   }
//   SECTION("both non empty"){
//     Waveform<int> emptywaveformInt(std::vector<int>{1,2,3,4}, 1.0);
//     Waveform<double> emptywaveformDouble(std::vector<double>{1.5,2.5,3.5,4.5}, 1.0);
//
//     CHECK((emptywaveformInt+emptywaveformDouble).getSamples() == std::vector<double>{2.5,4.5,6.5,8.5});
//   }
//
// }
