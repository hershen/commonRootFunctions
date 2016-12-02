#include <vector>

namespace filterFuncs
{
  
  //T is sampling interval
  std::vector<double> DfilterCR_RC ( std::vector<double> xV,  double tau = 1.,  double T = 0.);
  
  template <typename Type>
  std::vector<double> DfilterCR_RC2( std::vector<Type> xV,  double tau = 1.,  double T = 0.);
  
  template <typename Type>
  std::vector<double> DfilterCR_RC4( std::vector<Type> xV,  double tau = 1.,  double T = 0.);
  
  
  
  
  

}
