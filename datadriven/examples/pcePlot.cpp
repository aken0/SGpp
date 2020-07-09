#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
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
  of.open("plot_pce/" + path + ".txt", std::ios::out | std::ios::trunc);
  of << std::fixed;
  of << std::setprecision(std::numeric_limits<long double>::digits10 + 1);
  int points = 15000;
  int dim = 2;
  ee.printGrid(dim, 1, "plot_pce/" + path + "(base-level5).txt");
  std::cout << "base" << '\n';
  for (int i = 2; i <= points; i *= 2) {
    of << i << ',';
  }
  of << '\n';
  for (int i = 2; i <= points; i *= 2) {
    auto re = ee.adaptiveQuadrature(e, dim, i);
    of << re << ',';
  }
  ee.printAdaptiveGrid(e, dim, 400, "plot_pce/" + path + "(e-level5).txt");
  std::cout << "e" << '\n';
  of << '\n';
  for (int i = 2; i <= points; i *= 2) {
    auto re = ee.adaptiveQuadrature(f, dim, i);
    of << re << ',';
  }
  ee.printAdaptiveGrid(f, dim, 400, "plot_pce/" + path + "(f-level5).txt");
  std::cout << "f" << '\n';
  of << '\n';
  for (int i = 2; i <= points; i *= 2) {
    auto re = ee.adaptiveQuadrature(g, dim, i);
    of << re << ',';
  }
  ee.printAdaptiveGrid(g, dim, 400, "plot_pce/" + path + "(g-level5).txt");
  std::cout << "g" << '\n';
  of << '\n';
  for (int i = 2; i <= points; i *= 2) {
    auto re = ee.adaptiveQuadratureL2(e, dim, i);
    of << re << ',';
  }
  of << '\n';
  std::cout << "e" << '\n';
  for (int i = 2; i <= points; i *= 2) {
    auto re = ee.adaptiveQuadratureL2(f, dim, i);
    of << re << ',';
  }
  of << '\n';
  std::cout << "f" << '\n';
  for (int i = 2; i <= points; i *= 2) {
    auto re = ee.adaptiveQuadratureL2(g, dim, i);
    of << re << ',';
  }
  of << '\n';
  std::cout << "g" << '\n';
  of.close();
}
