#include <string>
#include <vector>

namespace myFuncs
{
  struct detectorRegion
  {
    std::string name; 
    double min;
    double max;  
  };

  const std::vector<detectorRegion> ECLdetectorRegions = { {"FWD",0.296706, 0.548033}, // = 17..31.4 deg
			        {"BRL",0.561996, 2.24624},  // = 32.2..128.7 deg
			        {"BWD",2.28115, 2.61799}  // = 130.7..150 deg
			      };
  

  std::string getECLdetectorRegion(const double theta)
  {
    for(const auto& ECLdetectorRegion : ECLdetectorRegions)
      if( theta > ECLdetectorRegion.min && theta < ECLdetectorRegion.max ) return ECLdetectorRegion.name;
    
    return "other";
    
  }
}