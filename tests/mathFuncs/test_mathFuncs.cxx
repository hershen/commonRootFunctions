#include "catch.hpp"

#include "mathFuncs.h"

#include "TGraph.h"

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

template <class T>
void testshiftVector_empty(const std::string &typeName) {
  SECTION("Empty vector " + typeName) {
    std::vector<T> emptyVector;
    CHECK(myFuncs::shiftVector(emptyVector, 1).size() == 0);
    CHECK(myFuncs::shiftVector(emptyVector, 1.5).size() == 0);
    CHECK(myFuncs::shiftVector(emptyVector, 1.5, 2.5).size() == 0);
  }
}
template <class T>
void testshiftVector_linear(const std::string &typeName) {

  SECTION("Linear vector " + typeName) {
    std::vector<T> inputs{1, 2, 3, 4, 5, 6};

    //---------------------------------------
    // Integer bins
    //---------------------------------------
    SECTION("Shift by integer bins") {
      SECTION("Shift by +1 bin") {
        const std::vector<double> expected{1, 1, 2, 3, 4, 5};
        const auto shifted = myFuncs::shiftVector(inputs, 1);
        CHECK(expected.size() == shifted.size());
        for (uint i = 0; i < inputs.size(); ++i)
          CHECK(shifted[i] == Approx(expected[i]));
      }
      SECTION("Shift by +3 bin") {
        const std::vector<double> expected{1, 1, 1, 1, 2, 3};
        const auto shifted = myFuncs::shiftVector(inputs, 3);
        CHECK(expected.size() == shifted.size());
        for (uint i = 0; i < inputs.size(); ++i)
          CHECK(shifted[i] == Approx(expected[i]));
      }
      SECTION("Shift by +5 bin") {
        const std::vector<double> expected{1, 1, 1, 1, 1, 1};
        const auto shifted = myFuncs::shiftVector(inputs, 5);
        CHECK(expected.size() == shifted.size());
        for (uint i = 0; i < inputs.size(); ++i)
          CHECK(shifted[i] == Approx(expected[i]));
      }
      SECTION("Shift by +10 bin") {
        const std::vector<double> expected{1, 1, 1, 1, 1, 1};
        const auto shifted = myFuncs::shiftVector(inputs, 10);
        CHECK(expected.size() == shifted.size());
        for (uint i = 0; i < inputs.size(); ++i)
          CHECK(shifted[i] == Approx(expected[i]));
      }
      SECTION("Shift by -1 bin") {
        const std::vector<double> expected{2, 3, 4, 5, 6, 6};
        const auto shifted = myFuncs::shiftVector(inputs, -1);
        CHECK(expected.size() == shifted.size());
        for (uint i = 0; i < inputs.size(); ++i)
          CHECK(shifted[i] == Approx(expected[i]));
      }
      SECTION("Shift by -3 bin") {
        const std::vector<double> expected{4, 5, 6, 6, 6, 6};
        const auto shifted = myFuncs::shiftVector(inputs, -3);
        CHECK(expected.size() == shifted.size());
        for (uint i = 0; i < inputs.size(); ++i)
          CHECK(shifted[i] == Approx(expected[i]));
      }
      SECTION("Shift by -5 bin") {
        const std::vector<double> expected{6, 6, 6, 6, 6, 6};
        const auto shifted = myFuncs::shiftVector(inputs, -5);
        CHECK(expected.size() == shifted.size());
        for (uint i = 0; i < inputs.size(); ++i)
          CHECK(shifted[i] == Approx(expected[i]));
      }
      SECTION("Shift by -10 bin") {
        const std::vector<double> expected{6, 6, 6, 6, 6, 6};
        const auto shifted = myFuncs::shiftVector(inputs, -10);
        CHECK(expected.size() == shifted.size());
        for (uint i = 0; i < inputs.size(); ++i)
          CHECK(shifted[i] == Approx(expected[i]));
      }
    }
    //---------------------------------------
    // Fractional bins
    //---------------------------------------
    SECTION("Shift by fractional bin") {
      SECTION("Shift by +0.5 bin") {
        const std::vector<double> expected{1, 1.5, 2.5, 3.5, 4.5, 5.5};
        const auto shifted = myFuncs::shiftVector(inputs, 0.5);
        CHECK(expected.size() == shifted.size());
        for (uint i = 0; i < inputs.size(); ++i)
          CHECK(shifted[i] == Approx(expected[i]));
      }
      SECTION("Shift by +1.5 bin") {
        const std::vector<double> expected{1, 1, 1.5, 2.5, 3.5, 4.5};
        const auto shifted = myFuncs::shiftVector(inputs, 1.5);
        CHECK(expected.size() == shifted.size());
        for (uint i = 0; i < inputs.size(); ++i)
          CHECK(shifted[i] == Approx(expected[i]));
      }
      SECTION("Shift by -0.5 bin") {
        const std::vector<double> expected{1.5, 2.5, 3.5, 4.5, 5.5, 6};
        const auto shifted = myFuncs::shiftVector(inputs, -0.5);
        CHECK(expected.size() == shifted.size());
        for (uint i = 0; i < inputs.size(); ++i)
          CHECK(shifted[i] == Approx(expected[i]));
      }
      SECTION("Shift by -1.5 bin") {
        const std::vector<double> expected{2.5, 3.5, 4.5, 5.5, 6, 6};
        const auto shifted = myFuncs::shiftVector(inputs, -1.5);
        CHECK(expected.size() == shifted.size());
        for (uint i = 0; i < inputs.size(); ++i)
          CHECK(shifted[i] == Approx(expected[i]));
      }
    }
  }
}

template <class T>
void testshiftVector_parabola(const std::string &typeName) {

  SECTION("Parabola vector " + typeName) {
    std::vector<T> inputs{9, 4, 1, 0, 1, 4, 9};

    //---------------------------------------
    // Integer bins
    //---------------------------------------
    SECTION("Shift by integer bins") {
      SECTION("Shift by +1 bin") {
        const std::vector<double> expected{9, 9, 4, 1, 0, 1, 4};
        const auto shifted = myFuncs::shiftVector(inputs, 1);
        CHECK(expected.size() == shifted.size());
        for (uint i = 0; i < inputs.size(); ++i)
          CHECK(shifted[i] == Approx(expected[i]));
      }
      SECTION("Shift by +3 bin") {
        const std::vector<double> expected{9, 9, 9, 9, 4, 1, 0};
        const auto shifted = myFuncs::shiftVector(inputs, 3);
        CHECK(expected.size() == shifted.size());
        for (uint i = 0; i < inputs.size(); ++i)
          CHECK(shifted[i] == Approx(expected[i]));
      }
      SECTION("Shift by -1 bin") {
        const std::vector<double> expected{4, 1, 0, 1, 4, 9, 9};
        const auto shifted = myFuncs::shiftVector(inputs, -1);
        CHECK(expected.size() == shifted.size());
        for (uint i = 0; i < inputs.size(); ++i)
          CHECK(shifted[i] == Approx(expected[i]));
      }
      SECTION("Shift by -3 bin") {
        const std::vector<double> expected{0, 1, 4, 9, 9, 9, 9};
        const auto shifted = myFuncs::shiftVector(inputs, -3);
        CHECK(expected.size() == shifted.size());
        for (uint i = 0; i < inputs.size(); ++i)
          CHECK(shifted[i] == Approx(expected[i]));
      }
    }
    //---------------------------------------
    // Fractional bins
    //---------------------------------------
    SECTION("Shift by fractional bin") {
      SECTION("Shift by +0.5 bin") {
        const std::vector<double> expected{9, 6.5, 2.5, 0.5, 0.5, 2.5, 6.5};
        const auto shifted = myFuncs::shiftVector(inputs, 0.5);
        CHECK(expected.size() == shifted.size());
        for (uint i = 0; i < inputs.size(); ++i)
          CHECK(shifted[i] == Approx(expected[i]));
      }
      SECTION("Shift by +1.5 bin") {
        const std::vector<double> expected{9, 9, 6.5, 2.5, 0.5, 0.5, 2.5};
        const auto shifted = myFuncs::shiftVector(inputs, 1.5);
        CHECK(expected.size() == shifted.size());
        for (uint i = 0; i < inputs.size(); ++i)
          CHECK(shifted[i] == Approx(expected[i]));
      }
      SECTION("Shift by -0.5 bin") {
        const std::vector<double> expected{6.5, 2.5, 0.5, 0.5, 2.5, 6.5, 9};
        const auto shifted = myFuncs::shiftVector(inputs, -0.5);
        CHECK(expected.size() == shifted.size());
        for (uint i = 0; i < inputs.size(); ++i)
          CHECK(shifted[i] == Approx(expected[i]));
      }
      SECTION("Shift by -1.5 bin") {
        const std::vector<double> expected{2.5, 0.5, 0.5, 2.5, 6.5, 9, 9};
        const auto shifted = myFuncs::shiftVector(inputs, -1.5);
        CHECK(expected.size() == shifted.size());
        for (uint i = 0; i < inputs.size(); ++i)
          CHECK(shifted[i] == Approx(expected[i]));
      }
    }
  }
}

TEST_CASE("Test shiftVector function", "[shiftVector]") {
  testshiftVector_empty<int>("int");
  testshiftVector_empty<double>("double");
  testshiftVector_linear<int>("int");
  testshiftVector_linear<double>("double");
}

TEST_CASE("Test parabola_xMax function", "[parabola_xMax]") {
  SECTION("Using raw numbers") {
    CHECK(myFuncs::parabola_xMax(0, 0) == Approx(0));
    CHECK(myFuncs::parabola_xMax(0, 1) == Approx(0));
    CHECK(myFuncs::parabola_xMax(10, 2) == Approx(-2.5));
    CHECK(myFuncs::parabola_xMax(10.0, 2.0) == Approx(-2.5));
  }

  SECTION("Using TFitResult") {
    // Assuming y = 1 + 10x + 2x^2
    std::vector<double> xs{0, 1, 2, 3, 4};
    std::vector<double> ys{1, 13, 29, 49, 73};
    TGraph graph(5, xs.data(), ys.data());
    auto fitResult = graph.Fit("pol2", "SQ");
    CHECK(myFuncs::parabola_xMax(fitResult) == Approx(-2.5));
  }
}

TEST_CASE("Test parabola_maxValue function", "[parabola_maxValue]") {
  SECTION("Using raw numbers") {
    CHECK(myFuncs::parabola_maxValue(0, 0, 0) == Approx(0));

    CHECK(myFuncs::parabola_maxValue(1, 10, 2) == Approx(-11.5));
    CHECK(myFuncs::parabola_maxValue(1, 2, 3) == Approx(2. / 3.));
    CHECK(myFuncs::parabola_maxValue(-0.5, 0.5, 1.5) == Approx(-0.54166666666666666667)); // Not best number to check
  }

  SECTION("Using TFitResult") {
    // Assuming y = 1 + 10x + 2x^2
    std::vector<double> xs{0, 1, 2, 3, 4};
    std::vector<double> ys{1, 13, 29, 49, 73};
    TGraph graph(5, xs.data(), ys.data());
    auto fitResult = graph.Fit("pol2", "SQ");
    CHECK(myFuncs::parabola_maxValue(fitResult, -2.5) == Approx(-11.5));
    CHECK(myFuncs::parabola_maxValue(fitResult) == Approx(-11.5));
  }
}

TEST_CASE("Test linear_crossValue function", "[linear_crossValue]") {
  SECTION("Using raw numbers") {
    CHECK(myFuncs::linear_crossValue(0, 0, 0) == Approx(0));
    CHECK(myFuncs::linear_crossValue(0, 1, 0) == Approx(0));
    CHECK(myFuncs::linear_crossValue(0, 1, 2) == Approx(2));
    CHECK(myFuncs::linear_crossValue(10, 2, 3) == Approx(-3.5));
  }

  SECTION("Using TFitResult") {
    // Assuming y = 10 + 2x
    std::vector<double> xs{0, 1, 2, 3, 4};
    std::vector<double> ys{10, 12, 14, 16, 18};
    TGraph graph(5, xs.data(), ys.data());
    auto fitResult = graph.Fit("pol1", "SQ");

    CHECK(myFuncs::linear_crossValue(fitResult, 10) == Approx(0.0).margin(1e-10));
    CHECK(myFuncs::linear_crossValue(fitResult, 0) == Approx(-5));
  }
}

TEST_CASE("Test sumElementByElement function", "[sumElementByElement]") {
  SECTION("Empty vectors") {
    std::vector<int> vectorInt;
    std::vector<double> vectorDouble;
    CHECK(myFuncs::sumElementByElement(vectorInt, vectorInt).size() == 0);
    CHECK(myFuncs::sumElementByElement(vectorDouble, vectorInt).size() == 0);
    CHECK(myFuncs::sumElementByElement(vectorInt, vectorDouble).size() == 0);
    CHECK(myFuncs::sumElementByElement(vectorDouble, vectorDouble).size() == 0);
  }

  SECTION("Non empty vectors") {
    std::vector<int> vectorInt{1, 2, 3, 4};
    std::vector<double> vectorDouble{1.5, 1.5, -1.5, 3.1};

    SECTION("int + int") {
      const std::vector<int> expected{2, 4, 6, 8};
      const auto result = myFuncs::sumElementByElement(vectorInt, vectorInt);
      CHECK(expected.size() == result.size());
      for (uint i = 0; i < result.size(); ++i)
        CHECK(result[i] == Approx(expected[i]));
    }
    SECTION("int + double") {
      const std::vector<double> expected{2.5, 3.5, 1.5, 7.1};
      const auto result = myFuncs::sumElementByElement(vectorInt, vectorDouble);
      CHECK(expected.size() == result.size());
      for (uint i = 0; i < result.size(); ++i)
        CHECK(result[i] == Approx(expected[i]));
    }
    SECTION("double + double") {
      const std::vector<double> expected{3, 3, -3, 6.2};
      const auto result = myFuncs::sumElementByElement(vectorDouble, vectorDouble);
      CHECK(expected.size() == result.size());
      for (uint i = 0; i < result.size(); ++i)
        CHECK(result[i] == Approx(expected[i]));
    }
  }
}

TEST_CASE("Test findMaxEvery_n function", "[findMaxEvery_n]") {
  SECTION("Empty vectors") {
    std::vector<int> vectorInt;
    std::vector<double> vectorDouble;
    CHECK(myFuncs::findMaxEvery_n(vectorInt.begin(), vectorInt.end(), 0) == vectorInt.begin());
    CHECK(myFuncs::findMaxEvery_n(vectorDouble.begin(), vectorDouble.end(), 0) == vectorDouble.begin());
    CHECK(myFuncs::findMaxEvery_n(vectorInt.begin(), vectorInt.end(), 3) == vectorInt.begin());
    CHECK(myFuncs::findMaxEvery_n(vectorDouble.begin(), vectorDouble.end(), 3) == vectorDouble.begin());
    CHECK(myFuncs::findMaxEvery_n(vectorInt, 0) == vectorInt.begin());
    CHECK(myFuncs::findMaxEvery_n(vectorDouble, 0) == vectorDouble.begin());
    CHECK(myFuncs::findMaxEvery_n(vectorInt, 3) == vectorInt.begin());
    CHECK(myFuncs::findMaxEvery_n(vectorDouble, 3) == vectorDouble.begin());
  }

  SECTION("Non empty vectors") {
    std::vector<double> vector{1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16,
                               0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0,  0,  0,  0,  0};

    SECTION("Find every 1") {
      CHECK(*myFuncs::findMaxEvery_n(vector.begin(), vector.end(), 1) == 16);
      CHECK(*myFuncs::findMaxEvery_n(vector, 1) == 16);
    }
    SECTION("Find every 2") {
      CHECK(*myFuncs::findMaxEvery_n(vector.begin(), vector.end(), 2) == 15);
      CHECK(*myFuncs::findMaxEvery_n(vector, 2) == 15);
    }
    SECTION("Find every 3") {
      CHECK(*myFuncs::findMaxEvery_n(vector.begin(), vector.end(), 3) == 16);
      CHECK(*myFuncs::findMaxEvery_n(vector, 3) == 16);
    }
    SECTION("Find every 4") {
      CHECK(*myFuncs::findMaxEvery_n(vector.begin(), vector.end(), 4) == 13);
      CHECK(*myFuncs::findMaxEvery_n(vector, 4) == 13);
    }
    SECTION("Find every 5") {
      CHECK(*myFuncs::findMaxEvery_n(vector.begin(), vector.end(), 5) == 16);
      CHECK(*myFuncs::findMaxEvery_n(vector, 5) == 16);
    }
    SECTION("Find every 6") {
      CHECK(*myFuncs::findMaxEvery_n(vector.begin(), vector.end(), 6) == 13);
      CHECK(*myFuncs::findMaxEvery_n(vector, 6) == 13);
    }
    SECTION("Find every 7") {
      CHECK(*myFuncs::findMaxEvery_n(vector.begin(), vector.end(), 7) == 15);
      CHECK(*myFuncs::findMaxEvery_n(vector, 7) == 15);
    }
    SECTION("Find every 8") {
      CHECK(*myFuncs::findMaxEvery_n(vector.begin(), vector.end(), 8) == 9);
      CHECK(*myFuncs::findMaxEvery_n(vector, 8) == 9);
    }
    SECTION("Find every 9") {
      CHECK(*myFuncs::findMaxEvery_n(vector.begin(), vector.end(), 9) == 10);
      CHECK(*myFuncs::findMaxEvery_n(vector, 9) == 10);
    }
  }
}

TEST_CASE("Test findMaxEvery_n_ThenBetween function", "[findMaxEvery_n_ThenBetween]") {
  SECTION("Empty vectors") {
    std::vector<int> vectorInt;
    std::vector<double> vectorDouble;
    CHECK(myFuncs::findMaxEvery_n_ThenBetween(vectorInt.begin(), vectorInt.end(), 0) == vectorInt.begin());
    CHECK(myFuncs::findMaxEvery_n_ThenBetween(vectorDouble.begin(), vectorDouble.end(), 0) == vectorDouble.begin());
    CHECK(myFuncs::findMaxEvery_n_ThenBetween(vectorInt.begin(), vectorInt.end(), 1) == vectorInt.begin());
    CHECK(myFuncs::findMaxEvery_n_ThenBetween(vectorDouble.begin(), vectorDouble.end(), 1) == vectorDouble.begin());
    CHECK(myFuncs::findMaxEvery_n_ThenBetween(vectorInt.begin(), vectorInt.end(), 3) == vectorInt.begin());
    CHECK(myFuncs::findMaxEvery_n_ThenBetween(vectorDouble.begin(), vectorDouble.end(), 3) == vectorDouble.begin());
  }

  SECTION("Non empty vectors") {
    std::vector<double> vector{0,  1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15,
                               16, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0,  0,  0,  0,  0};

    SECTION("Find every 0") {
      CHECK(myFuncs::findMaxEvery_n_ThenBetween(vector.begin(), vector.end(), 0) == vector.begin());
      CHECK(myFuncs::findMaxEvery_n_ThenBetween(vector, 0) == vector.begin());
    }
    SECTION("Find every 1") {
      CHECK(*myFuncs::findMaxEvery_n_ThenBetween(vector.begin(), vector.end(), 1) == 16);
      CHECK(*myFuncs::findMaxEvery_n_ThenBetween(vector, 1) == 16);
    }
    SECTION("Find every 10") {
      CHECK(*myFuncs::findMaxEvery_n_ThenBetween(vector.begin(), vector.end(), 10) == 16);
      CHECK(*myFuncs::findMaxEvery_n_ThenBetween(vector, 10) == 16);
    }
    SECTION("Find every 20") {
      CHECK(*myFuncs::findMaxEvery_n_ThenBetween(vector.begin(), vector.end(), 20) == 16);
      CHECK(*myFuncs::findMaxEvery_n_ThenBetween(vector, 20) == 16);
    }
    SECTION("Out of range") {
      SECTION("Find every 50") {
        CHECK(*myFuncs::findMaxEvery_n_ThenBetween(vector.begin(), vector.end(), 50) == 16);
        CHECK(*myFuncs::findMaxEvery_n_ThenBetween(vector, 50) == 16);
      }
    }
  }
}
