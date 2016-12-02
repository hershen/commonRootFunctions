#include <map>

namespace myFuncs
{
  
  //--------------------------------------------------------------------------------------------
  //increaseMapValue
  //********************************************************************************************
  //Increase map value related to key by deltaValue.
  //If the map doesn't contain the provided key, a new one will be added with value deltaValue (i.e. assuming value starts at 0).
  //--------------------------------------------------------------------------------------------
  template <typename keyType, typename valueType>
  void increaseMapValue(std::map<keyType,valueType>& map, const keyType key, const valueType deltaValue )
  {
    auto returnPair = map.emplace(key,deltaValue);
    
    //returnPair.second is bool (true if key was added).
    if( !returnPair.second ) returnPair.first->second += deltaValue;  
  }
  
  template <typename keyType, typename valueType>
  std::vector<std::pair<keyType, valueType> > sortMapByValue(const std::map<keyType, valueType>& map)
  {
    std::vector<std::pair<keyType, valueType>> pairs;
    for (auto pair : map) pairs.push_back(pair);
    
    std::sort(pairs.begin(), pairs.end(), [=](std::pair<keyType, valueType>& a, std::pair<keyType, valueType>& b)
    {
	return a.second < b.second;
    }
    );
    
    return pairs;
  }
}