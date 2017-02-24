#pragma once

#include <algorithm>
#include <iterator>
#include <functional>

namespace myFuncs {

	//----------------------------------------------------------
	//Inputs:
	//  first - forward itirator point to first element to check.
	//  last - forward itirator pointing to last element to check.
	//  threshold - threshold to compare to
	
	//Returns a pair of itirators, the first points to the first element that's greater than threshold. The second points to the last element that's smaller than threshold.
	
	//If all elements are smaller than T, both returned itirators will point to last.
	//If all elements are larger than T, both returned itirators will point to first.
	//----------------------------------------------------------
	template <class FwItir, class T >
	std::pair<FwItir,FwItir> findFirstBiggerLastSmaller(FwItir first, FwItir last, T threshold) {
		
		//Find first element greather than threshold
		//Note - can't user lower/upper_bound because it assumes the container is sorted.
		const auto firstGreater = std::find_if(first, last, [threshold](typename std::iterator_traits<FwItir>::value_type value){return value > threshold;});//first element greater than threshold
		
		//Find last element smaller than threshold - starting from back.
		//Note - the in firstGreater + 1 is because reverse_iterator points to one before the base iterator
		const auto lastSmaller = std::find_if(std::reverse_iterator<FwItir>(last), std::reverse_iterator<FwItir>(firstGreater + 1), [threshold](typename std::iterator_traits<FwItir>::value_type value){return value < threshold;} );
				
		return std::make_pair(firstGreater, (lastSmaller + 1).base());
		
	}	

} //myFuncs