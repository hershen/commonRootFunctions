#include "mathFuncs.h"
#include "TF1.h"
#include <iostream>

bool testcalcResiduals() 
{
	bool pass = true;
	TF1 quadFunc("quadFunc","x*x",-10,10);
	std::vector<double> xValues = {-3, -1, 0, 1, 5};

	//------------------------------------------
	//All residuals should be zero
	//------------------------------------------
	std::vector<double> goodYvalues = {9, 1, 0, 1, 25};
	TF1 a = myFuncs::CFD(quadFunc, 1, 0.5);
	std::vector<double> goodResiduals = myFuncs::calcResiduals<double,double>(xValues, goodYvalues, 0.5, quadFunc);
	
	for(auto goodResidual : goodResiduals) 
		if( std::abs(goodResidual) > 1e-6)
		{
			std::cout << "testcalcResiduals: expecting residual = 0, obtained " << goodResidual << std::endl;
			pass = false;
		}
	
	
	//------------------------------------------
	//Some residuals are not zero. Single std
	//------------------------------------------
	std::vector<double> badYvalues = {8, 1, 3, -4, -1, 0, 0, 0};  //Larger on purpose
	std::vector<double> expectedResiduals = {2, 0, -6, 10, 52};
	
	std::vector<double> badResiduals = myFuncs::calcResiduals(xValues, badYvalues, 0.5, quadFunc);
	for(int iResid = 0; iResid < badResiduals.size(); ++iResid) {
		
		if( std::abs(badResiduals[iResid] - expectedResiduals[iResid]) > 1e-6)
		{
			std::cout << "testcalcResiduals: expecting residual = " << expectedResiduals[iResid] << ", obtained " << badResiduals[iResid] << std::endl;
			pass = false;
		}
	}
	
	
	//------------------------------------------
	//Some residuals are not zero. Vector std 
	//------------------------------------------
	std::vector<double> stds = {1, 0.5, 2, 3, 0.01, 0, 0, 0};  //Larger on purpose
	expectedResiduals = std::vector<double>{1, 0, -3./2., 5./3., 2600.};
	
	badResiduals = myFuncs::calcResiduals(xValues, badYvalues, stds, quadFunc);
	for(int iResid = 0; iResid < badResiduals.size(); ++iResid) {
		if( std::abs(badResiduals[iResid] - expectedResiduals[iResid]) > 1e-6)
		{
			std::cout << "testcalcResiduals: func val = " << quadFunc.Eval(xValues[iResid]) << ", value = " << badYvalues[iResid] << ", expecting residual = " << expectedResiduals[iResid] << ", obtained " << badResiduals[iResid] << std::endl;
			pass = false;
		}
	}
	
	
	try{
		myFuncs::calcResiduals(xValues, std::vector<double>{1}, 1, quadFunc);
	}
	catch (const std::invalid_argument&) {
		std::cout << "in catch" << std::endl;
	}
	
	try{
		myFuncs::calcResiduals(xValues, xValues, std::vector<double>{1}, quadFunc);
	}
	catch (const std::invalid_argument&) {
		std::cout << "in catch" << std::endl;
	}
	
	return pass;
}

void testMathFuncs()
{
	auto result = testcalcResiduals();
	std::cout << "testcalcResiduals test result = " << result << std::endl;
}

int main()
{
	testMathFuncs();
}