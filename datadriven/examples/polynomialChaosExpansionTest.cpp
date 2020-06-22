#include <math.h>

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sgpp/datadriven/tools/PolynomialChaosExpansion.hpp>
#include <sgpp/globaldef.hpp>
#include <sgpp_base.hpp>
#include <utility>
#include <vector>

#include "sgpp/base/datatypes/DataVector.hpp"
// function to be integrated
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
  // std::cout << "integral(monteCarloQuad): " << ee.monteCarloQuad(e, 1000000) << std::endl;
  // std::cout << "integral(sparseGridQuadrature): " << ee.sparseGridQuadrature(e, 1, 2) << '\n';
  // std::cout << "integral(adaptiveQuadrature): " << ee.adaptiveQuadrature(e, 1, 2, 5) << '\n';
  std::fstream of;
  std::string path;
  std::cout << "enter path: " << '\n';
  std::cin >> path;
  of.open("data/" + path + ".txt", std::ios::out | std::ios::trunc);
  of << std::fixed;
  of << std::setprecision(9);
  int iters = 10;
  int dim = 2;
  for (int i = 2; i <= iters; ++i) {
    auto re = ee.sparseGridQuadrature(e, dim, i);
    of << re << ',';
    std::cout << i << '\n';
  }
  of << '\n';
  for (int i = 2; i <= iters; ++i) {
    auto re = ee.sparseGridQuadrature(f, dim, i);
    of << re << ',';
    std::cout << i << '\n';
  }
  of << '\n';
  for (int i = 2; i <= iters; ++i) {
    auto re = ee.sparseGridQuadrature(g, dim, i);
    of << re << ',';
    std::cout << i << '\n';
  }
  of << '\n';
  /*
  for (int i = 2; i < 15; ++i) {
    of << ee.adaptiveQuadrature(e, 2, i, 5) << ',';
    std::cout << i << '\n';
  }*/
  of.close();
  /*
  std::cout << "coefficients: " << std::endl;
    ee.calculateCoefficients();
    auto stuff = ee.getCoefficients();
    for (auto entry : stuff) {
      std::cout << entry << " ";
    }
    std::cout << std::endl;
    std::cout << "evaluation at (2,2): " << ee.evalExpansion(sgpp::base::DataVector{2, 2})
              << std::endl;
              */
  return 0;
}
