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
  return (vec[0] * vec[0] * vec[0] * std::exp(-std::pow(vec[0], 2) / 2) -
          vec[1] * vec[1] * std::exp(-std::pow(vec[1], 2) / 2) +
          vec[2] * std::exp(-std::pow(vec[2], 2) / 2));
}
double f(const sgpp::base::DataVector& vec) { return 1.0; }
int main() {
  sgpp::base::DistributionsVector dists;
  auto dist1 = std::make_shared<sgpp::base::DistributionNormal>(0, 1);
  auto dist2 = std::make_shared<sgpp::base::DistributionNormal>(0, 1);
  auto dist3 = std::make_shared<sgpp::base::DistributionNormal>(0, 1);
  dists.push_back(dist1);
  dists.push_back(dist2);
  dists.push_back(dist3);
  sgpp::datadriven::PolynomialChaosExpansion ee =
      sgpp::datadriven::PolynomialChaosExpansion(e, 5, dists);
  std::fstream of;
  std::string path;
  std::cout << "enter path: " << '\n';
  std::cin >> path;
  of.open("plot_pce/" + path + ".txt", std::ios::out | std::ios::trunc);
  of << std::fixed;
  of << std::setprecision(std::numeric_limits<double>::digits10 + 1);
  int points = 500;
  int dim = 3;
  std::cout << "base" << '\n';
  for (int i = 20; i <= points; i *= 1.5) {
    of << i << ',';
  }
  of << '\n';
  for (int i = 20; i <= points; i *= 1.5) {
    auto re = ee.adaptiveQuadrature(e, dim, i);
    of << re << ',';
  }
  std::cout << "e" << '\n';
  of << '\n';
  for (int i = 20; i <= points; i *= 1.5) {
    auto re = ee.adaptiveQuadrature(f, dim, i);
    of << re << ',';
  }
  std::cout << "f" << '\n';
  of << '\n';
  for (int i = 20; i <= points; i *= 1.5) {
    auto re = ee.adaptiveQuadratureL2(e, dim, i);
    of << re << ',';
  }
  of << '\n';
  std::cout << "e" << '\n';
  for (int i = 20; i <= points; i *= 1.5) {
    auto re = ee.adaptiveQuadratureL2(f, dim, i);
    of << re << ',';
  }
  of << '\n';
  std::cout << "f" << '\n';
  // sparseGridQuadrature
  for (int i = 20; i <= points; i *= 1.5) {
    auto re = ee.sparseGridQuadrature(e, dim, i);
    of << re << ',';
  }
  std::cout << "e2" << '\n';
  of << '\n';
  for (int i = 20; i <= points; i *= 1.5) {
    auto re = ee.sparseGridQuadrature(f, dim, i);
    of << re << ',';
  }
  of << '\n';
  std::cout << "f2" << '\n';
  for (int i = 20; i <= points; i *= 1.5) {
    auto re = ee.sparseGridQuadratureL2(e, dim, i);
    of << re << ',';
  }
  of << '\n';
  std::cout << "e2" << '\n';
  for (int i = 20; i <= points; i *= 1.5) {
    auto re = ee.sparseGridQuadratureL2(f, dim, i);
    of << re << ',';
  }
  of << '\n';
  std::cout << "f2" << '\n';

  for (auto method : {"sparseGrid", "adaptiveGrid"}) {
    for (auto function : {e, f}) {
      for (auto order : {1, 3, 5}) {
        for (int i = 20; i <= points; i *= 1.5) {
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
  for (int i = 20; i <= points; i *= 1.5) {
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
