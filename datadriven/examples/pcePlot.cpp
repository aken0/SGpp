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
double e(const sgpp::base::DataVector& vec) {
  return sin(vec[0]) + 3 * vec[0] * pow(sin(vec[1]), 3) +
         5 * exp(-100 * (pow(vec[0] - 0.1, 2) + pow(vec[1], 2)));
}
double f(const sgpp::base::DataVector& vec) { return 1.0; }
double g(const sgpp::base::DataVector& vec) {
  // return 1 + (std::sin(vec[0]) + std::cos(vec[1])) / std::exp(vec[1]);
  return sin(vec[0]) + 3 * vec[0] * pow(sin(vec[1]), 3);
  // return std::exp(-200 * (std::pow(vec[0] - 0.8, 2) + std::pow(vec[1] - 0.8, 2)));
}
int main() {
  sgpp::base::DistributionsVector dists;
  auto dist1 = std::make_shared<sgpp::base::DistributionNormal>(0, .1);
  auto dist2 = std::make_shared<sgpp::base::DistributionNormal>(0, .1);
  dists.push_back(dist1);
  dists.push_back(dist2);
  sgpp::datadriven::PolynomialChaosExpansion ee =
      sgpp::datadriven::PolynomialChaosExpansion(e, 5, dists);
  std::fstream of;
  std::string path;
  std::cout << "enter path: " << '\n';
  std::cin >> path;
  of.open("plot_pce/" + path + ".txt", std::ios::out | std::ios::trunc);
  of << std::fixed;
  of << std::setprecision(std::numeric_limits<double>::digits10 + 1);
  int points = 20;
  int dim = 2;
  ee.printGrid(dim, 1, "plot_pce/" + path + "(base-level5).txt");
  std::cout << "base" << '\n';
  for (int i = 50; i <= points; i *= 1.5) {
    of << i << ',';
  }
  of << '\n';
  for (int i = 50; i <= points; i *= 1.5) {
    auto re = ee.adaptiveQuadratureWeighted(e, dim, i, 100);
    of << re << ',';
  }
  std::cout << "e" << '\n';
  of << '\n';
  for (int i = 50; i <= points; i *= 1.5) {
    auto re = ee.adaptiveQuadratureWeighted(f, dim, i, 100);
    of << re << ',';
  }

  ee.printAdaptiveGrid(e, dim, 500, "plot_pce/" + path + "(e-level5).txt");
  // ee.printGrid(dim, 200, "plot_pce/" + path + "(e-level5).txt");
  ee.printAdaptiveGrid(f, dim, 500, "plot_pce/" + path + "(f-level5).txt");
  std::cout << "f" << '\n';
  of << '\n';
  for (int i = 50; i <= points; i *= 1.5) {
    auto re = ee.adaptiveQuadratureWeighted(g, dim, i, 100);
    of << re << ',';
  }
  ee.printAdaptiveGrid(g, dim, 500, "plot_pce/" + path + "(g-level5).txt");
  std::cout << "g" << '\n';
  of << '\n';
  for (int i = 50; i <= points; i *= 1.5) {
    auto re = ee.adaptiveQuadratureL2(e, dim, i);
    of << re << ',';
  }
  of << '\n';
  std::cout << "e" << '\n';
  for (int i = 50; i <= points; i *= 1.5) {
    auto re = ee.adaptiveQuadratureL2(f, dim, i);
    of << re << ',';
  }
  of << '\n';
  std::cout << "f" << '\n';
  for (int i = 50; i <= points; i *= 1.5) {
    auto re = ee.adaptiveQuadratureL2(g, dim, i);
    of << re << ',';
  }
  // sparseGridQuadrature
  of << '\n';
  std::cout << "g" << '\n';
  for (int i = 50; i <= points; i *= 1.5) {
    auto re = ee.sparseGridQuadrature(e, dim, i, 100);
    of << re << ',';
  }
  std::cout << "e2" << '\n';
  of << '\n';
  for (int i = 50; i <= points; i *= 1.5) {
    auto re = ee.sparseGridQuadrature(f, dim, i, 100);
    of << re << ',';
  }
  std::cout << "f2" << '\n';
  of << '\n';
  for (int i = 50; i <= points; i *= 1.5) {
    auto re = ee.sparseGridQuadrature(g, dim, i, 100);
    of << re << ',';
  }
  std::cout << "g2" << '\n';
  of << '\n';
  for (int i = 50; i <= points; i *= 1.5) {
    auto re = ee.sparseGridQuadratureL2(e, dim, i);
    of << re << ',';
  }
  of << '\n';
  std::cout << "e2" << '\n';
  for (int i = 50; i <= points; i *= 1.5) {
    auto re = ee.sparseGridQuadratureL2(f, dim, i);
    of << re << ',';
  }
  of << '\n';
  std::cout << "f2" << '\n';
  for (int i = 50; i <= points; i *= 1.5) {
    auto re = ee.sparseGridQuadratureL2(g, dim, i);
    of << re << ',';
  }
  // Monte-Carlo
  of << '\n';
  std::cout << "g2" << '\n';
  for (int i = 50; i <= points; i *= 1.5) {
    sgpp::base::DataVector result(500);
    for (int j = 0; j < 500; ++j) {
      result[j] = ee.monteCarloQuad(e, i);
    }
    double re = result.sum() / static_cast<double>(result.size());
    of << re << ',';
  }
  of << '\n';
  std::cout << "emc" << '\n';
  for (int i = 50; i <= points; i *= 1.5) {
    sgpp::base::DataVector result(500);
    for (int j = 0; j < 500; ++j) {
      result[j] = ee.monteCarloQuad(f, i);
    }
    double re = result.sum() / static_cast<double>(result.size());
    of << re << ',';
  }
  of << '\n';
  std::cout << "fmc" << '\n';
  for (int i = 50; i <= points; i *= 1.5) {
    sgpp::base::DataVector result(500);
    for (int j = 0; j < 500; ++j) {
      result[j] = ee.monteCarloQuad(g, i);
    }
    double re = result.sum() / static_cast<double>(result.size());
    of << re << ',';
  }
  of << '\n';
  std::cout << "gmc" << '\n';

  for (auto method : {"sparseGrid", "adaptiveWeighted"}) {
    for (auto function : {e, f, g}) {
      for (auto order : {1, 3, 10}) {
        for (int i = 50; i <= points; i *= 1.5) {
          auto ee = sgpp::datadriven::PolynomialChaosExpansion(function, order, dists);
          ee.calculateCoefficients(i, method);
          auto re = ee.getL2Error(i, method);
          of << re << ',';
          ee.clearCoefficients();
        }
        std::cout << "pce" << '\n';
        of << '\n';
      }
    }
  }
  of << '\n';
  std::cout << "done" << '\n';
  for (int i = 50; i <= points; i *= 1.5) {
    std::unique_ptr<sgpp::base::Grid> grid(sgpp::base::Grid::createLinearBoundaryGrid(dim));
    sgpp::base::GridStorage& gridStorage = grid->getStorage();
    int j = 1;
    while (gridStorage.getSize() < i) {
      gridStorage.clear();
      grid->getGenerator().regular(j);
      ++j;
    }
    of << gridStorage.getSize() << ',';
  }
  of.close();
}
