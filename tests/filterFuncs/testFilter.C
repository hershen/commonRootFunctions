#include "filterFuncs.h"

void testFilter() {
  const std::vector<double> nom = {4, 1, 0.5, 3, 5};
  const std::vector<double> denom = {2.0, 1.5};
  const std::vector<double> inputs = {1.0, 2.0, 3.0, 4.0};

  // const double sum = std::inner_product(inputs.begin(), inputs.end(), inputs.rbegin(), 0.0);
  // std::cout << sum << std::endl;

  // return;
  const auto output = myFuncs::DSP::filter(nom, denom, inputs);

  std::for_each(output.begin(), output.end(), [](const auto &x) { std::cout << x << "\n"; });
}
