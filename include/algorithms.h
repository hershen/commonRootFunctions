#pragma once

#include <algorithm>
#include <functional>
#include <iterator>

namespace myFuncs {

//----------------------------------------------------------
// Works on unsorted containers (in contrast to upper/lower bound!!
// Inputs:
//  first - forward itirator point to first element to check.
//  last - forward itirator pointing to last element to check.
//  threshold - threshold to compare to

// Returns an itirator to the first element that's greater than threshold.
// If all elements are smaller than T, returned itirator will point to last.
// If all elements are larger than T, returned itirator will point to first.
//----------------------------------------------------------
template <class FwItir, class T>
FwItir findFirstBigger(FwItir first, FwItir last, T threshold) {
  return std::find_if(first, last, [threshold](typename std::iterator_traits<FwItir>::value_type value) {
    return value > threshold;
  }); //first element greater than threshold
}

//----------------------------------------------------------
// Works on unsorted containers (in contrast to upper/lower bound!!
// Inputs:
//  first - forward itirator point to first element to check.
//  last - forward itirator pointing to last element to check.
//  threshold - threshold to compare to

// find_if is supposed to work the range [first, last).
// But it seems that this works in the range [first, last]

// Returns an itirator to the last element that's smaller than threshold.
// If all elements are larger than T, returned itirator will point to last.
//----------------------------------------------------------
template <class FwItir, class T>
FwItir findLastSmaller(FwItir first, FwItir last, T threshold) {

  const auto rbegin = std::reverse_iterator<FwItir>(last);
  const auto rend = std::reverse_iterator<FwItir>(first);

  // 		std::cout << "rbegin = "<< *rbegin << ", rend = " << *rend << std::endl;
  const auto reverseLastSmaller = std::find_if(
      rbegin, rend, [threshold](typename std::iterator_traits<FwItir>::value_type value) { return value < threshold; });

  // If nothing found return last (this is because we're going to do return ( reverseLastSmaller + 1).base(); and base moves one
  // after the reverse iterator)
  if (reverseLastSmaller == rend)
    return last;

  // The +1 is becuase base moves one after the reverse iterator so we're going one before, before doing base
  return (reverseLastSmaller + 1).base();
}

//----------------------------------------------------------
// Works on unsorted containers (in contrast to upper/lower bound!!
// Inputs:
//  first - forward itirator point to first element to check.
//  last - forward itirator pointing to last element to check.
//  threshold - threshold to compare to

// Returns a pair of itirators, the first points to the first element that's greater than threshold. The second points to the last
// element that's smaller than threshold.

// If all elements are smaller than T, both returned itirators will point to last.
// If all elements are larger than T, both returned itirators will point to first.
//----------------------------------------------------------
template <class FwItir, class T>
std::pair<FwItir, FwItir> findFirstBiggerLastSmaller(FwItir first, FwItir last, T threshold) {

  // Find first element greather than threshold
  // Note - can't user lower/upper_bound because it assumes the container is sorted.
  const auto firstGreater = findFirstBigger(first, last, threshold); // first element greater than threshold

  return std::make_pair(firstGreater, findLastSmaller(firstGreater, last, threshold));
}
} // namespace myFuncs
