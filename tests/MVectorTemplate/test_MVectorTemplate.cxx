#include "catch.hpp"
#include <iostream>

#include "MVectorTemplate.h"

#include "TCanvas.h"
#include "TF1.h"
#include "TH1D.h"
#include "TRandom3.h"
#include "TStopwatch.h"
#include "TStyle.h"
#include "mathFuncs.h"


#include "histFuncs.h"
#define cline std::cout << "line = " << __LINE__ << std::endl;

using namespace myFuncs;

SCENARIO("Empty vector", "[MVectorTemplate]") {
  GIVEN("An empty vector") {
		std::vector<double> empty(1e3, 0.0);
		WHEN("Adding many of these") {
			MVectorTemplate myTemplate(empty, 1.0);
			TCanvas* canvas = new TCanvas("canvas", "", 0, 0 1200, 900);
			myTemplate.getTF1()->Draw();
			waitForDoubleclick();
			// myTemplate.addVector(empty, 1e-2, *myTemplate.getTF1());
		}

	}
}
