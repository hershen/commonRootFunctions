#include "mathFuncs.h"

std::vector<double> numbers {1.0, 1.113, 1.87, 9.97, 1.000113, 1.000813, 5.15, 3.05, 3.11,  8.0, 0.000113, 0.000813, 0.813, 0.113, 0.000000000000117, 0.000000000813, 101, 350, 355, 367, 100000000, 39300000, 366};

void testRound_35rule() {
	for(const auto number : numbers) {
		std::cout << "testRound_35rule(" << number << ") = " << myFuncs::round_35rule(number) << std::endl;
		
	}
	
}