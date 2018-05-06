#include "catch.hpp"

#include "LHEreader.h"
#include <iostream>

using namespace myFuncs;
using Catch::Equals;

inline bool operator==(const Particle& lhs, const Particle& rhs) {
  return lhs.pdg == rhs.pdg and lhs.px == Approx(rhs.px )  and lhs.py == Approx(rhs.py ) and
         lhs.pz == Approx(rhs.pz )  and lhs.E == Approx(rhs.E )  and lhs.mass == Approx(rhs.mass );
}
inline bool operator!=(const Particle& lhs, const Particle& rhs) {
   return !(lhs == rhs); }

std::vector<Particle> firstEvent{ {-11,  +0.0000000000e+00, +0.0000000000e+00, +5.2896999753e+00, 5.2897000000e+00, 5.1100000000e-04},
{11, -0.0000000000e+00, -0.0000000000e+00, -5.2896999753e+00, 5.2897000000e+00, 5.1100000000e-04},
{9000006, +3.8579933588e+00, -3.8086585203e-01, +1.3593358906e+00, 6.4712415784e+00, 5.0000001949e+00},
{22, -3.8579933588e+00, +3.8086585203e-01, -1.3593358906e+00, 4.1081584216e+00, 0.0000000000e+00},
{22, +2.6531798569e-01, +1.7471644229e+00, -5.8991031695e-01, 1.8630543030e+00, 0.0000000000e+00},
{22, +3.5926753731e+00, -2.1280302749e+00, +1.9492462075e+00, 4.6081872754e+00, 0.0000000000e+00}};

std::vector<Particle> secondEvent{ {-11,  +0.0000000000e+00, +0.0000000000e+00, +5.2896999753e+00, 5.2897000000e+00, 5.1100000000e-04},
{11, -0.0000000000e+00, -0.0000000000e+00, -5.2896999753e+00, 5.2897000000e+00, 5.1100000000e-04},
{9000006, +9.1155157064e-01, -2.3802561235e+00, +3.2218659963e+00, 6.4712415345e+00, 5.0000001020e+00},
{22, -9.1155157064e-01, +2.3802561235e+00, -3.2218659963e+00, 4.1081584655e+00, 0.0000000000e+00},
{22, -9.4807115451e-02, -2.8754016736e+00, +4.1264363460e+00, 5.0303479096e+00, 0.0000000000e+00},
{22, +1.0063586861e+00, +4.9514555009e-01, -9.0457034968e-01, 1.4408936249e+00, 0.0000000000e+00}};

SCENARIO("Read lhe file", "[lheFuncs]") {
  GIVEN("a bad filename") {
    const std::string filename = "bad.lhe";
    WHEN("Reading a file") {
      THEN("It should throw") { CHECK_THROWS(LHEreader(filename)); }
    }
  }
  GIVEN("an lhe file") {
    const std::string filename = "/home/hershen/PhD/Root/commonRootFunctions/tests/lheFuncs/sampleLhefile.lhe";
    WHEN("Reading a file") {
      LHEreader lheReader(filename);
      THEN("") {
        WHEN("Reading the first event") {
          auto particles = lheReader.getEventParticles();
          THEN("particles should match first event") {
            REQUIRE_THAT( particles, Equals(firstEvent) );
           }

           WHEN("Reading the second event") {
             particles = lheReader.getEventParticles();
             THEN("particles should match second event") {
               REQUIRE_THAT( particles, Equals(secondEvent) );
              }
           }
        }

      }
    }

  }
}
