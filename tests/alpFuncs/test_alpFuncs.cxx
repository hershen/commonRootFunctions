#include "catch.hpp"

#include "alpFuncs.h"

using namespace myFuncs;

SCENARIO("Check inEMC", "[inEMC]") {
  GIVEN("inEMC") {
    WHEN("Checking barrel section") {
      CHECK(inEMC(0, 0, 0) == false);
      CHECK(inEMC(100, 0, 0) == true);
      CHECK(inEMC(0, 90, 170) == false);
      CHECK(inEMC(0, 95, 170) == true);
      CHECK(inEMC(0, 90, -170) == false);
      CHECK(inEMC(0, 95, -170) == true);
    }
    WHEN("Checking endcap section") {
      CHECK(inEMC(55, 0, 190) == false);
      CHECK(inEMC(69, 0, 190) == true);
      CHECK(inEMC(100, 90, 170) == true);
      CHECK(inEMC(56, 0, 195.35) == true);
    }
  }
}
