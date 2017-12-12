#include "catch.hpp"

#include "filterFuncs.h"
#include "mathFuncs.h"

using Catch::Matchers::Equals;
using namespace Catch::Matchers;
using myFuncs::DSP::filtfilt;

const std::vector<double> inputs{1, 2, 3, 4, 5, 6, 5, 4, 3, 2, 1};

TEST_CASE("Test filtfilt function", "[filtfilt]") {

  SECTION("Empty denom and nom") {
    const std::vector<double> nom;
    const std::vector<double> denom;
    REQUIRE(filtfilt(nom, denom, inputs).size() == 0);
  }

  SECTION("Empty denom") {
    const std::vector<double> nom{1};
    const std::vector<double> denom;
    REQUIRE(filtfilt(nom, denom, inputs).size() == 0);
  }

  SECTION("Empty nom") {
    const std::vector<double> nom;
    const std::vector<double> denom{1};
    REQUIRE(filtfilt(nom, denom, inputs).size() == 0);
  }

  SECTION("Identity") {
    const std::vector<double> nom{1};
    const std::vector<double> denom{1};

    CHECK_THAT(filtfilt(nom, denom, inputs), Equals(inputs));
  }

  SECTION("Times 4") {
    const std::vector<double> nom{2};
    const std::vector<double> denom{1};

    CHECK_THAT(filtfilt(nom, denom, inputs), Equals(myFuncs::scaleVector(inputs, 4.0)));
  }

  SECTION("Times 0.5") {
    const std::vector<double> nom{1};
    const std::vector<double> denom{2};

    CHECK_THAT(filtfilt(nom, denom, inputs), Equals(myFuncs::scaleVector(inputs, 0.25)));
  }

  SECTION("Test1") {
    const std::vector<double> nom{1.1, 0.01};
    const std::vector<double> denom{1.3, 0.7};
    const std::vector<double> expected{0.26498722288220061000, 0.66605629631788910000, 0.84731747385076850000,
                                       1.36597437072092930000, 1.29617249519140070000, 2.29868932211725330000,
                                       1.29663770017605850000, 1.36485992361482670000, 0.84952204325215963000,
                                       0.66188946727560383000, 0.27276472745023994000};

    const auto filtered = filtfilt(nom, denom, inputs);
    REQUIRE(filtered.size() == expected.size());
    for (uint i = 0; i < inputs.size(); ++i) {
      REQUIRE(std::abs(filtered[i] - expected[i]) < 0.17);
    }
  }

  SECTION("Test2") {
    const std::vector<double> nom{1.1, 0.9};
    const std::vector<double> denom{0.9, 1.1, 0.7, 0.001};
    const std::vector<double> expected{
        0.66532406028368873000, -0.84242373439698515000, 4.01456302133739480000, 1.11685760359906980000,
        1.60549364253181140000, 6.07097873944464710000,  1.23630521598657460000, 1.69015482249061220000,
        3.51488221941513770000, -0.69952138765133287000, 1.02478765211951290000,
    };

    const auto filtered = filtfilt(nom, denom, inputs);
    REQUIRE(filtered.size() == expected.size());
    for (uint i = 0; i < inputs.size(); ++i)
      REQUIRE(std::abs(filtered[i] - expected[i]) < 0.12);
  }
}
