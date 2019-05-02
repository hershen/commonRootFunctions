#include <vector>

class LinearInterp {
public:
  LinearInterp(const std::vector<double>& xVals, const std::vector<double>& yVals):
  xVals(xVals), 
  yVals(yVals) 
  {};
  
  double interpolate(const double x) const;
  
  std::vector<double> getXvals() const {return xVals;}
  std::vector<double> getYvals() const {return yVals;}
  
private:
  std::vector<double> xVals;
  std::vector<double> yVals;
};
