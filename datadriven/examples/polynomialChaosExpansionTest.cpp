#include <math.h>

#include <cmath>
#include <iomanip>
#include <iostream>
#include <sgpp/datadriven/tools/PolynomialChaosExpansion.hpp>
#include <sgpp/globaldef.hpp>
#include <sgpp_base.hpp>

#include "sgpp/base/datatypes/DataVector.hpp"
// function to be integrated
double e(const sgpp::base::DataVector& vec) {
  return std::pow(vec.get(0), 2) + std::pow(vec.get(1), 3);
}
int main() {
  sgpp::datadriven::PolynomialChaosExpansion ee = sgpp::datadriven::PolynomialChaosExpansion(
      e, 4, std::vector<std::string>{"legendre", "legendre"},
      std::vector<std::pair<double, double>>{std::pair<double, double>{-1, 1},
                                             std::pair<double, double>{-1, 1}},
      0, 0);
  std::cout << "integral: " << ee.monteCarloQuad(e, 1000000) << std::endl;
  std::cout << "coefficients: " << std::endl;

  ee.calculateCoefficients();
  auto stuff = ee.getCoefficients();
  for (auto entry : stuff) {
    std::cout << entry << " ";
  }
  std::cout << std::endl;
  std::cout << "evaluation at (2,2): " << ee.evalExpansion(sgpp::base::DataVector{2, 2})
            << std::endl;
  return 0;
}
