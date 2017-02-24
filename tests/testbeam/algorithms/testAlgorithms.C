#include "algorithms.h"

bool testFindFirstBiggerLastSmaller() 
{
	bool pass = true;
	
	//----------------------------------------
	//Test mixed vector
	//----------------------------------------
	std::vector<double> vector = {1.1, 2.2, 3.3, 4.1, 3.3, -10.0, 5.2, 6.1, 2, 7, 8, 9};
	auto result = myFuncs::findFirstBiggerLastSmaller (vector.begin(), vector.end(), 4);
	
	if( std::distance(vector.begin(), result.first) != 3 ) {
		std::cout << "std::distance(vector.begin(), result.first) = " << std::distance(vector.begin(), result.first) << " != 3" << std::endl;
		pass = false;
	}
	if( std::distance(vector.begin(), result.second) != 8 ) {
		std::cout << "std::distance(vector.begin(), result.second) = " << std::distance(vector.begin(), result.second) << " != 8" << std::endl;
		pass = false;
	}
	
	
	
	//----------------------------------------
	//Test all smaller than threshold
	//----------------------------------------
	vector = {0, 0, 0, 0, 0};
	
	result = myFuncs::findFirstBiggerLastSmaller (vector.begin(), vector.end(), 4);
	
	if( std::distance(vector.begin(), result.first) != vector.size() ) pass = false;
	if( std::distance(vector.begin(), result.second) != vector.size() ) pass = false;
	
	if( !pass ) std::cout << "testFindFirstBiggerLastSmaller failed on all smaller vector" << std::endl;
	
	//----------------------------------------
	//Test all larger than threshold
	//----------------------------------------
	vector = {10, 10, 10, 10, 10};
	
	result = myFuncs::findFirstBiggerLastSmaller (vector.begin(), vector.end(), 4);
	
	if( std::distance(vector.begin(), result.first) != 0 ) pass = false;
	if( std::distance(vector.begin(), result.second) != 0 ) pass = false;
	
	if( !pass ) std::cout << "testFindFirstBiggerLastSmaller failed on all larger vector" << std::endl;
	
	return pass;
}

void testAlgorithms()
{
	auto result = testFindFirstBiggerLastSmaller();
	std::cout << "testFindFirstBiggerLastSmaller test result = " << result << std::endl;
}

int main()
{
	testFindFirstBiggerLastSmaller();
}