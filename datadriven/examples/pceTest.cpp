#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/datadriven/tools/PolynomialChaosExpansion.hpp>
#include <sgpp/globaldef.hpp>
#include <sgpp_base.hpp>
#include <string>
#include <utility>
#include <vector>
// functions to be integrated
double e(const sgpp::base::DataVector& vec) { return vec[0] * vec[0] * vec[0] - vec[1] * vec[1]; }
double f(const sgpp::base::DataVector& vec) { return 1.0; }
double g(const sgpp::base::DataVector& vec) {
  return 1 + (std::sin(vec[0]) + std::cos(vec[1])) / std::exp(vec[1]);
}
int main() {
  std::fstream of;
  std::string path;
  std::cout << "enter path: " << '\n';
  std::cin >> path;
  of.open("plot/" + path + "(pce).txt", std::ios::out | std::ios::trunc);
  of << std::fixed;
  of << std::setprecision(9);
  for (int i = 1; i < 20; ++i) {
    sgpp::datadriven::PolynomialChaosExpansion ee =
        sgpp::datadriven::PolynomialChaosExpansion(g, i,
                                                   std::vector<sgpp::datadriven::distributionType>{
                                                       sgpp::datadriven::distributionType::Uniform,
                                                       sgpp::datadriven::distributionType::Uniform},
                                                   std::vector<std::pair<double, double>>{
                                                       std::pair<double, double>{-1, 1},
                                                       std::pair<double, double>{-1, 1},
                                                   },
                                                   0, 0);
    auto re = ee.calculateCoefficients();
    for (auto entry : re) {
      of << entry << ',';
    }
  }
  of.close();
}
