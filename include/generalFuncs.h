#include <unordered_map>

namespace myFuncs {

//--------------------------------------------------------------------------------------------
// increaseMapValue
//********************************************************************************************
// Increase map value related to key by deltaValue.
// If the map doesn't contain the provided key, a new one will be added with value deltaValue (i.e. assuming value starts at 0).
//--------------------------------------------------------------------------------------------
template <typename keyType, typename valueType>
void increaseMapValue(std::unordered_map<keyType, valueType> &map, const keyType key, const valueType deltaValue) {
  auto returnPair = map.emplace(key, deltaValue);

  // returnPair.second is bool (true if key was added).
  if (!returnPair.second)
    returnPair.first->second += deltaValue;
}

template <typename keyType, typename valueType>
std::vector<std::pair<keyType, valueType>> sortMapByValue(const std::unordered_map<keyType, valueType> &map) {
  std::vector<std::pair<keyType, valueType>> pairs;
  for (auto pair : map)
    pairs.push_back(pair);

  std::sort(pairs.begin(), pairs.end(),
            [=](std::pair<keyType, valueType> &a, std::pair<keyType, valueType> &b) { return a.second < b.second; });

  return pairs;
}

template <class T>
std::vector<T> multiplyElementByElement(const std::vector<T> &vec) {
  return vec;
}

template <class T1, class T2>
std::vector<decltype(T1() * T2())> multiplyElementByElement(const std::vector<T1> &vec1, const std::vector<T2> &vec2) {
  std::vector<decltype(T1() * T2())> returnVec;
  returnVec.reserve(vec1.size());
  std::transform(vec1.begin(), vec1.end(),
                 vec2.begin(),                  // Second vector
                 std::back_inserter(returnVec), // Insert result to returnVec
                 std::multiplies<void>());
  return returnVec;
}

// I don't know how to extract the type of T1*T2*Ts[0]*Ts[1]...
// won't work if (for example) vector<int>, vector<int>, vector<double>.
template <class T1, class T2, class... Ts>
std::vector<decltype(T1() * T2())> multiplyElementByElement(const std::vector<T1> &first, const std::vector<T2> &second, Ts &... others) {
  return multiplyElementByElement(multiplyElementByElement(first, second), others...);
}

} // namespace myFuncs
