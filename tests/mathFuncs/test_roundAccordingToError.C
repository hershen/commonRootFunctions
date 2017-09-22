#include <boost/lexical_cast.hpp>

#include "mathFuncs.h"


const std::vector<double> xToRounds {1234, 1234.567898888, 1234.8, 8888.8888, 1.0, 1.1, 5.0, 5.5, 0.01, 0.11111111, 0.8, 0.08888888, 1234.1111111, 1234.01111, 1.234e-9};

void try_roundAccordingToError(const double xToRound, const double error, const double answer, const double epsilon = 1e-9) {
	
	const double functionReturns = myFuncs::roundAccordingToError(xToRound, error);
	
	if( std::abs(functionReturns - answer) > epsilon )
		std::cout << "WRONG : ";
	else 
		std::cout << "RIGHT : ";
	
	std::cout << "roundAccordingToError(" << xToRound << ", " << error << ") = " << functionReturns << ". Should be " << answer << std::endl;
	
}

void test_roundAccordingToError() {
	const double y = 0.603;
	
	for(auto xToRound : xToRounds)
		std::cout << std::setprecision(9) << "roundAccordingToError(" << xToRound << ", " << y << ") = " << myFuncs::roundAccordingToError(xToRound, y) << std::endl;

	std::cout << std::endl << std::endl;
	
	try_roundAccordingToError(888.888, 10, 890);
	try_roundAccordingToError(888.888, 33, 889);
	try_roundAccordingToError(888.888, 1, 889);
	try_roundAccordingToError(888.888, 0.8, 888.9);
	try_roundAccordingToError(888.888, 0.088, 888.888);
	
	//negative number
	try_roundAccordingToError(-888.888, 10, -890);
	try_roundAccordingToError(-888.888, 33, -889);
	try_roundAccordingToError(-888.888, 1, -889);
	try_roundAccordingToError(-888.888, 0.8, -888.9);
	try_roundAccordingToError(-888.888, 0.088, -888.888);
	
	//Providing wrong answer
	try_roundAccordingToError(888.888, 0.088, -30);
}

void lex(const double x) {
	std::cout << boost::lexical_cast<std::string>(x) << std::endl;
}
