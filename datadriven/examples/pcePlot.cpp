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
  sgpp::datadriven::PolynomialChaosExpansion ee = sgpp::datadriven::PolynomialChaosExpansion(
      e, 5,
      std::vector<sgpp::datadriven::distributionType>{sgpp::datadriven::distributionType::Uniform,
                                                      sgpp::datadriven::distributionType::Uniform},
      std::vector<std::pair<double, double>>{
          std::pair<double, double>{-1, 1},
          std::pair<double, double>{-1, 1},
      },
      0, 0);
  std::fstream of;
  std::string path;
  std::cout << "enter path: " << '\n';
  std::cin >> path;
  of.open("plot/" + path + ".txt", std::ios::out | std::ios::trunc);
  of << std::fixed;
  of << std::setprecision(9);
  int iters = 8;
  int dim = 2;
  ee.printGrid(dim, 5, "plot/" + path + "(base-level5).txt");
  for (int i = 2; i <= iters; ++i) {
    auto re = ee.adaptiveQuadrature(e, dim, i, 10);
    of << re << ',';
    std::cout << i << '\n';
  }
  ee.printAdaptiveGrid(e, dim, 5, 10, "plot/" + path + "(e-level5).txt");
  of << '\n';
  for (int i = 2; i <= iters; ++i) {
    auto re = ee.adaptiveQuadrature(f, dim, i, 10);
    of << re << ',';
    std::cout << i << '\n';
  }
  ee.printAdaptiveGrid(f, dim, 5, 10, "plot/" + path + "(f-level5).txt");
  of << '\n';
  for (int i = 2; i <= iters; ++i) {
    auto re = ee.adaptiveQuadrature(g, dim, i, 10);
    of << re << ',';
    std::cout << i << '\n';
  }
  ee.printAdaptiveGrid(g, dim, 5, 10, "plot/" + path + "(g-level5).txt");
  of << '\n';
  for (int i = 2; i <= iters; ++i) {
    auto re = ee.adaptiveQuadratureL2(e, dim, i, 10);
    of << re << ',';
    std::cout << i << '\n';
  }
  of << '\n';
  for (int i = 2; i <= iters; ++i) {
    auto re = ee.adaptiveQuadratureL2(f, dim, i, 10);
    of << re << ',';
    std::cout << i << '\n';
  }
  of << '\n';
  for (int i = 2; i <= iters; ++i) {
    auto re = ee.adaptiveQuadratureL2(g, dim, i, 10);
    of << re << ',';
    std::cout << i << '\n';
  }
  of << '\n';
  of.close();
}
