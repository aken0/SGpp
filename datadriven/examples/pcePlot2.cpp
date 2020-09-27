#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/datadriven/tools/PolynomialChaosExpansion.hpp>
#include <sgpp/globaldef.hpp>
#include <sgpp_base.hpp>
#include <sgpp_optimization.hpp>
#include <string>
#include <utility>
#include <vector>
// functions to be integrated
double e(const sgpp::base::DataVector& vec) {
  return 1 - ((4 * vec[1]) / (5 * 225 * vec[0])) -
         ((vec[2] * vec[2]) / (25 * 225 * vec[0] * vec[0]));
}
int main() {
  sgpp::base::DistributionsVector dists;

  auto dist1 = std::make_shared<sgpp::base::DistributionLogNormal>(5, .5);
  auto dist2 = std::make_shared<sgpp::base::DistributionNormal>(2000, 400);
  auto dist3 = std::make_shared<sgpp::base::DistributionNormal>(500, 100);
  sgpp::base::DataVector lb{exp(.5), -1600, -400};
  sgpp::base::DataVector ub{exp(9.5), 5600, 1400};
  dists.push_back(dist1);
  dists.push_back(dist2);
  dists.push_back(dist3);

  class Functe : public sgpp::base::ScalarFunction {
   public:
    explicit Functe(int dim) : sgpp::base::ScalarFunction(dim) {}
    double eval(const sgpp::base::DataVector& vec) { return e(vec); }
    virtual void clone(std::unique_ptr<sgpp::base::ScalarFunction>& clone) const {
      clone = std::unique_ptr<sgpp::base::ScalarFunction>(new Functe(*this));
    }
  };

  sgpp::datadriven::PolynomialChaosExpansion ee =
      sgpp::datadriven::PolynomialChaosExpansion(e, 5, dists);
  std::fstream of;
  std::string path;
  std::cout << "enter path: " << '\n';
  std::cin >> path;
  of.open("plot_pce/" + path + ".txt", std::ios::out | std::ios::trunc);
  of << std::fixed;
  of << std::setprecision(std::numeric_limits<double>::digits10 + 1);
  int points = 400;
  int dim = 3;
  std::vector<int> gridPoints;

  for (int i = 50; i <= points; i *= 2) {
    std::unique_ptr<sgpp::base::Grid> grid(sgpp::base::Grid::createLinearBoundaryGrid(dim));
    sgpp::base::GridStorage& gridStorage = grid->getStorage();
    int j = 1;
    while (gridStorage.getSize() < i) {
      gridStorage.clear();
      grid->getGenerator().regular(j);
      ++j;
    }
    gridPoints.push_back(gridStorage.getSize());
  }

  for (int i : gridPoints) {
    of << i << ',';
  }
  of << '\n';

  for (int i : gridPoints) {
    auto re = ee.adaptiveQuadratureWeighted(e, dim, i, 100);
    of << re << ',';
  }
  of << '\n';

  for (int i : gridPoints) {
    auto re = ee.adaptiveQuadratureL2(e, dim, i);
    of << re << ',';
  }
  of << '\n';

  // sparseGridQuadrature
  for (int i : gridPoints) {
    auto re = ee.sparseGridQuadrature(e, dim, i, 100);
    of << re << ',';
  }
  of << '\n';

  for (int i : gridPoints) {
    auto re = ee.sparseGridQuadratureL2(e, dim, i);
    of << re << ',';
  }
  of << '\n';

  for (auto method : {"sparseGrid", "adaptiveWeighted"}) {
    for (auto function : {e}) {
      for (auto order : {1, 3, 5}) {
        for (int i : gridPoints) {
          auto ee = sgpp::datadriven::PolynomialChaosExpansion(function, order, dists);
          ee.calculateCoefficients(i, method);
          auto re = ee.getL2Error(i, method);
          of << re << ',';
          ee.clearCoefficients();
        }
        of << '\n';
      }
    }
  }

  for (int i : gridPoints) {
    auto e1 = std::make_shared<Functe>(dim);
    sgpp::optimization::SplineResponseSurface surface(e1, lb, ub,
                                                      sgpp::base::GridType::NakBsplineExtended);
    surface.regularByPoints(i);
    auto gen = [&dists, &surface]() {
      sgpp::base::DataVector randvec = dists.sample();
      return std::pow(e(randvec) - surface.eval(randvec), 2);
    };
    size_t num = 100000;
    sgpp::base::DataVector results(num);
    std::generate(results.begin(), results.end(), gen);
    of << (results.sum() / static_cast<double>(num)) << ',';
  }
  of << '\n';

  for (int i : gridPoints) {
    auto e1 = std::make_shared<Functe>(dim);
    sgpp::optimization::SplineResponseSurface surface(e1, lb, ub,
                                                      sgpp::base::GridType::NakBsplineExtended);
    surface.surplusAdaptive(i, 1);
    auto gen = [&dists, &surface]() {
      sgpp::base::DataVector randvec = dists.sample();
      return std::pow(e(randvec) - surface.eval(randvec), 2);
    };
    size_t num = 100000;
    sgpp::base::DataVector results(num);
    std::generate(results.begin(), results.end(), gen);
    of << (results.sum() / static_cast<double>(num)) << ',';
  }
  of << '\n';

  for (auto method : {"sparseGrid", "adaptiveWeighted"}) {
    for (auto function : {e}) {
      for (auto order : {1, 3, 5}) {
        for (int i : gridPoints) {
          auto ee = sgpp::datadriven::PolynomialChaosExpansion(function, order, dists);
          ee.calculateCoefficients(i, method);
          of << ee.getMean(i, method) << ',';
          ee.clearCoefficients();
        }
        of << '\n';
      }
    }
  }

  for (int i : gridPoints) {
    auto e1 = std::make_shared<Functe>(dim);
    sgpp::optimization::SplineResponseSurface surface(e1, lb, ub,
                                                      sgpp::base::GridType::NakBsplineExtended);
    surface.regularByPoints(i);
    of << surface.getMean(dists, 100) << ',';
  }
  of << '\n';

  for (int i : gridPoints) {
    auto e1 = std::make_shared<Functe>(dim);
    sgpp::optimization::SplineResponseSurface surface(e1, lb, ub,
                                                      sgpp::base::GridType::NakBsplineExtended);
    surface.surplusAdaptive(i, 1);
    of << surface.getMean(dists, 100) << ',';
  }
  of << '\n';

  for (auto method : {"sparseGrid", "adaptiveWeighted"}) {
    for (auto function : {e}) {
      for (auto order : {1, 3, 5}) {
        for (int i : gridPoints) {
          auto ee = sgpp::datadriven::PolynomialChaosExpansion(function, order, dists);
          ee.calculateCoefficients(i, method);
          of << ee.getVariance(i, method) << ',';
          ee.clearCoefficients();
        }
        of << '\n';
      }
    }
  }

  for (int i : gridPoints) {
    auto e1 = std::make_shared<Functe>(dim);
    sgpp::optimization::SplineResponseSurface surface(e1, lb, ub,
                                                      sgpp::base::GridType::NakBsplineExtended);
    surface.regularByPoints(i);
    of << surface.getVariance(dists, 100)[0] << ',';
  }
  of << '\n';

  for (int i : gridPoints) {
    auto e1 = std::make_shared<Functe>(dim);
    sgpp::optimization::SplineResponseSurface surface(e1, lb, ub,
                                                      sgpp::base::GridType::NakBsplineExtended);
    surface.surplusAdaptive(i, 1);
    of << surface.getVariance(dists, 100)[0] << ',';
  }
  of << '\n';
  of.close();
}
