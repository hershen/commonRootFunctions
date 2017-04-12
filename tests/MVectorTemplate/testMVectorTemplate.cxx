#include <iostream>

#include "MVectorTemplate.h"

#include "TF1.h"
#include "TCanvas.h"
#include "TRandom3.h"
#include "TStopwatch.h"
#include "TH1D.h"
#include "mathFuncs.h"
#include "TStyle.h"
#include "TRandom3.h"

#define cline std::cout << "line = " << __LINE__ << std::endl;


const double pedestalOffset = 8;
const int indexOffset = 0;
const double ampScale = 0.5;

std::vector<double> v1;
std::vector<double> v2;
std::vector<double> times;

bool testVectorTemplateNoXshift() {
	bool pass = true;
	v1 = {-0.5,0.5,-0.5,0.5,-0.5,0.5,-0.5,0.5,-0.5,0.5,-0.5,0.5,-0.5,0.5,-0.5,0.5,-0.5,0.5,-0.5,0.5,-0.5,0.5,-0.5,0.5,-0.5,0.5,-0.5,0.5,-0.5,0.5,2,3,4,5,6,7,8,8,8,7,6,5,4,3,2,1,0,-1};
	
	for (size_t i = 0; i < v1.size(); ++i) {		
		v2.push_back(ampScale * v1[i] + pedestalOffset);
	}
	
	
	//Makes v2 inconsistent with v1.
	v2[37] += 1;
	
	//===========================================
	//Average v1 and v2 once
	myFuncs::MVectorTemplate vecTemplate(v1,2.0);
	vecTemplate.addVector(v2, 0.1);
	
	auto templateValues = vecTemplate.getTemplateValues();
	
	std::vector<double> expectedAverage {-0.0556693, 0.0556693, -0.0556693, 0.0556693, -0.0556693, 0.0556693, -0.0556693, 0.0556693, -0.0556693, 0.0556693, -0.0556693, 0.0556693, -0.0556693, 0.0556693, -0.0556693, 0.0556693, -0.0556693, 0.0556693, -0.0556693, 0.0556693, -0.0556693, 0.0556693, -0.0556693, 0.0556693, -0.0556693, 0.0556693, -0.0556693, 0.0556693, -0.0556693, 0.0556693, 0.222677, 0.334016, 0.445355, 0.556693, 0.668032, 0.77937, 0.890709, 1, 0.890709, 0.77937, 0.668032, 0.556693, 0.445355, 0.334016, 0.222677, 0.111339, 0};
	
	for(size_t idx = 0; idx < templateValues.size(); ++idx ) {
		if( std::abs(templateValues[idx] - expectedAverage[idx]) > 1e-5) {
			std::cout << "average of v1 and v2[" << idx << "] = " << templateValues[idx] << ", expecting " << expectedAverage[idx] << std::endl;
			pass = false;
		}
		
	}
	//===========================================
	//Average v1 and v2 many time and assume output will be very similar to v2.
	for(size_t iTime = 0; iTime < 1000; ++iTime) {
		vecTemplate.addVector(v2, 0.1);
	}
	
	//Make vectorTemplate from v2, in order to compare to average
	myFuncs::MVectorTemplate v2Template(v2,2.0);
	auto v2TemplateVaues = v2Template.getTemplateValues();
	auto averageVaues = vecTemplate.getTemplateValues();
	
	for(size_t idx = 0; idx < averageVaues.size(); ++idx ) {
		if( std::abs(averageVaues[idx] - v2TemplateVaues[idx]) > 6e-4) {
			std::cout << "average of v1 and v2 many times[" << idx << "] = " << averageVaues[idx] << ", expecting " << v2TemplateVaues[idx] << std::endl;
			pass = false;
		}
		
	}
	
	
	return pass;
}



void testMVectorTemplate()
{
  bool pass = testVectorTemplateNoXshift();
	if(!pass) std::cout << "testVectorTemplateNoXshift() FAILED!!!" << std::endl;
	else std::cout << "testVectorTemplateNoXshift() SUCEEDED!!!" << std::endl;
}
  

int main(/*int argc, char *argv[]*/)
{

  testMVectorTemplate();
  return 0;
}