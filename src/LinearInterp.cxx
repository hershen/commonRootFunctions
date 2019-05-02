#include "LinearInterp.h"

#include <algorithm>
#include <stdexcept>

double LinearInterp::interpolate(const double x) const {
  
  if(xVals.empty()) {
    throw std::invalid_argument("X is empty");
  }
  
  if(x<xVals.front()) {
    return yVals.front();
  } else if (x>xVals.back()) {
    return yVals.back();
  }    
  
  const auto x1_it = std::upper_bound(xVals.begin(), xVals.end(), x);
  const uint idx2 = std::distance(xVals.begin(), x1_it);
  
  const auto x1 = xVals[idx2-1];
  const auto x2 = xVals[idx2];
  const auto y1 = yVals[idx2-1];
  const auto y2 = yVals[idx2];
  
  return y1 + (y2-y1)/(x2-x1)*(x-x1);
}
