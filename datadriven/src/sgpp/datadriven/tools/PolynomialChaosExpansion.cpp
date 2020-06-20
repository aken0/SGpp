// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <functional>
#include <map>
#include <numeric>
#include <random>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/datadriven/tools/PolynomialChaosExpansion.hpp>
#include <sgpp/globaldef.hpp>
#include <sgpp_base.hpp>
#include <string>
#include <utility>
#include <vector>
namespace sgpp {
namespace datadriven {
PolynomialChaosExpansion::PolynomialChaosExpansion(std::function<double(const base::DataVector&)> f,
                                                   int order, std::vector<distributionType> types,
                                                   std::vector<std::pair<double, double>> ranges,
                                                   double alpha, double beta)
    : func(f), order(order), types(types), ranges(ranges), alpha(alpha), beta(beta) {
  // hermite,jacobi,legendre,laguerre,genlaguerre
  this->weights = std::vector<std::function<double(double)>>{
      {[](double x) { return std::exp(-std::pow(x, 2) / 2); }},
      {[this](double x) -> double {
        return std::pow((1 - x), this->alpha) * std::pow((1 + x), this->beta);
      }},
      {[](double x) { return 1.0; }},
      {[](double x) { return std::exp(-x); }},
      {[this](double x) -> double { return std::pow(x, this->alpha) * std::exp(-x); }}};
  this->denoms = std::vector<std::function<double(double)>>{
      {[](double j) { return std::sqrt(2 * M_PI) * std::tgamma(j + 1); }},
      {[this](double j) {
        return (
            (std::pow(2, this->alpha + this->beta + 1) / (2 * j + this->alpha + this->beta + 1)) *
            ((std::tgamma(j + this->alpha + 1) * std::tgamma(j + this->beta + 1)) /
             (std::tgamma(j + this->alpha + this->beta + 1) * std::tgamma(j + 1))));
      }},
      {[](double j) { return 2 / ((2 * j) + 1); }},
      {[](double j) { return 1.0; }},
      {[this](double j) { return std::tgamma(j + this->alpha + 1) / std::tgamma(j + 1); }}};

  this->evals = std::vector<std::function<double(double, double)>>{
      {[this](int n, double x) { return evalHermite(n, x); }},
      {[this](int n, double x) { return evalJacobi(n, x); }},
      {[this](int n, double x) { return evalLegendre(n, x); }},
      {[this](int n, double x) { return evalLaguerre(n, x); }},
      {[this](int n, double x) { return evalGenLaguerre(n, x); }}};
}
PolynomialChaosExpansion::~PolynomialChaosExpansion() {}
double PolynomialChaosExpansion::evalLegendre(int n, double x) {
  if (n == 0) {
    return 1.0;
  } else if (n == 1) {
    return x;
  } else {
    double next = 0.0;
    // n-1
    double last = 1.0;
    // n
    double curr = x;
    for (double i = 1.0; i < n; ++i) {
      next = ((2.0 * i + 1.0) / (i + 1.0)) * x * curr - (i / (i + 1.0)) * last;
      last = curr;
      curr = next;
    }
    return next;
  }
}
double PolynomialChaosExpansion::evalHermite(int n, double x) {
  if (n == 0) {
    return 1.0;
  } else if (n == 1) {
    return x;
  } else {
    double next = 0.0;
    double last = 1.0;
    double curr = x;
    for (int i = 1.0; i < n; ++i) {
      next = x * curr - i * last;
      last = curr;
      curr = next;
    }
    return next;
  }
}
double PolynomialChaosExpansion::evalLaguerre(int n, double x) {
  if (n == 0) {
    return 1.0;
  } else if (n == 1) {
    return 1 - x;
  } else {
    double next = 0.0;
    double last = 1.0;
    double curr = 1 - x;
    for (double i = 1.0; i < n; ++i) {
      next = ((2.0 * i + 1.0 - x) * curr - i * last) / (i + 1);
      last = curr;
      curr = next;
    }
    return next;
  }
  return 0;
}

double PolynomialChaosExpansion::evalJacobi(int n, double x) {
  if (n == 0) {
    return 1.0;
  } else if (n == 1) {
    return (alpha + 1) + (alpha + beta + 2) * ((x - 1) / 2);
  } else {
    double next = 0.0;
    double last = 1.0;
    double curr = (alpha + 1) + (alpha + beta + 2) * ((x - 1) / 2);
    for (double i = 2.0; i <= n; ++i) {
      double q1 =
          ((2 * i + alpha + beta - 1) * ((2 * i + alpha + beta) * (2 * i + alpha + beta - 2) * x +
                                         std::pow(alpha, 2) - std::pow(beta, 2))) /
          (2 * i * (i + alpha + beta) * (2 * i + alpha + beta - 2));
      double q2 = (2 * (i + alpha - 1) * (i + beta - 1) * (2 * i + alpha + beta)) /
                  (2 * i * (i + alpha + beta) * (2 * i + alpha + beta - 2));
      next = q1 * curr - q2 * last;
      last = curr;
      curr = next;
    }
    return next;
  }
  return 0;
}

double PolynomialChaosExpansion::evalGenLaguerre(int n, double x) {
  if (n == 0) {
    return 1.0;
  } else if (n == 1) {
    return 1 + alpha - x;
  } else {
    double next = 0.0;
    // n-1
    double last = 1.0;
    // n
    double curr = 1 + alpha - x;
    for (double i = 1.0; i < n; ++i) {
      next = ((2.0 * i + 1.0 + alpha - x) * curr - (i + alpha) * last) / (i + 1);
      last = curr;
      curr = next;
    }
    return next;
  }
  return 0;
}
double PolynomialChaosExpansion::monteCarloQuad(
    std::function<double(const base::DataVector&)> funct, size_t n) {
  std::vector<std::uniform_real_distribution<double>> dists(ranges.size());
  double prod = 1;
  for (auto pair : ranges) {
    prod *= (pair.second - pair.first);
  }
  double factor = prod * (1.0 / static_cast<double>(n));
  std::random_device dev;
  std::mt19937_64 mersenne{dev()};
  for (std::vector<std::pair<double, double>>::size_type i = 0; i < ranges.size(); ++i) {
    dists[i] = std::uniform_real_distribution<double>{ranges[i].first, ranges[i].second};
  }
  auto gen = [&funct, &dists, &mersenne]() {
    std::vector<double> randvec(dists.size());
    for (std::vector<std::uniform_real_distribution<double>>::size_type j = 0; j < dists.size();
         ++j) {
      randvec[j] = dists[j](mersenne);
    }
    return funct(base::DataVector(randvec));
  };
  base::DataVector results(n);
  std::generate(results.begin(), results.end(), gen);
  return factor * results.sum();
}
// wip
double PolynomialChaosExpansion::sparseGridQuadrature(
    std::function<double(const base::DataVector&)> funct, int dim, int level) {
  auto numfunc = [&funct](const base::DataVector& input,
                          std::vector<std::pair<double, double>>& ranges) {
    base::DataVector temp(input.size());
    for (int i = 0; i < input.size(); ++i) {
      temp[i] = input[i] * (ranges[i].second - ranges[i].first) + ranges[i].first;
    }
    return funct(temp);
  };
  std::unique_ptr<sgpp::base::Grid> grid(sgpp::base::Grid::createLinearBoundaryGrid(dim));
  sgpp::base::GridStorage& gridStorage = grid->getStorage();
  // std::cout << "dimensionality:        " << gridStorage.getDimension() << std::endl;

  grid->getGenerator().regular(level);
  std::cout << "number of grid points: " << gridStorage.getSize() << std::endl;

  /**
   * Calculate the surplus vector alpha for the interpolant of \f$
   * f(x)\f$.  Since the function can be evaluated at any
   * point. Hence. we simply evaluate it at the coordinates of the
   * grid points to obtain the nodal values. Then we use
   * hierarchization to obtain the surplus value.
   *
   */
  sgpp::base::DataVector evals(gridStorage.getSize());
  base::DataVector p(dim);
  for (size_t i = 0; i < gridStorage.getSize(); i++) {
    sgpp::base::GridPoint& gp = gridStorage.getPoint(i);
    for (int j = 0; j < dim; ++j) {
      p[j] = gp.getStandardCoordinate(j);
    }
    evals[i] = numfunc(p, ranges);
  }
  std::cout << "Hierarchizing...\n\n";
  sgpp::base::DataVector coeffs(evals.getSize());
  sgpp::base::HierarchisationSLE hierSLE(*grid);
  sgpp::base::sle_solver::Eigen sleSolver;

  // solve linear system
  if (!sleSolver.solve(hierSLE, evals, coeffs)) {
    std::cout << "Solving failed, exiting.\n";
    return 1;
  }
  /*std::unique_ptr<base::OperationHierarchisation>(
      sgpp::op_factory::createOperationHierarchisation(*grid))
      ->doHierarchisation(alpha);
      */

  // direct quadrature
  std::unique_ptr<sgpp::base::OperationQuadrature> opQ(
      sgpp::op_factory::createOperationQuadrature(*grid));
  double res = opQ->doQuadrature(coeffs);
  double prod = 1;
  for (auto pair : ranges) {
    prod *= (pair.second - pair.first);
  }
  return res * prod;
}
// work in progress
double PolynomialChaosExpansion::adaptiveQuadrature(
    std::function<double(const base::DataVector&)> funct, int dim, int level, int steps) {
  auto numfunc = [&funct](const base::DataVector& input,
                          std::vector<std::pair<double, double>>& ranges) {
    base::DataVector temp(input.size());
    for (int i = 0; i < input.size(); ++i) {
      temp[i] = (input[i]) * (ranges[i].second - ranges[i].first) + ranges[i].first;
    }
    return funct(temp);
  };

  std::unique_ptr<sgpp::base::Grid> grid(sgpp::base::Grid::createLinearGrid(dim));
  sgpp::base::GridStorage& gridStorage = grid->getStorage();

  grid->getGenerator().regular(level);

  /**
   * Create coefficient vector with size corresponding to the grid size.
   * Initially, all the values are set to zero.
   */
  base::DataVector alpha(gridStorage.getSize());

  /**
   * Create a vector for storing (possibly expensive) function evaluations at
   * each gridpoint.
   */
  base::DataVector funEvals(gridStorage.getSize());
  for (size_t i = 0; i < gridStorage.getSize(); i++) {
    base::GridPoint& gp = gridStorage.getPoint(i);
    base::DataVector vec(dim);
    for (int j = 0; j < dim; ++j) {
      vec[j] = gp.getStandardCoordinate(j);
    }
    funEvals[i] = numfunc(vec, ranges);
  }

  /**
   * create a vector for storing newly added points by their sequence id.
   */
  std::vector<size_t> addedPoints;

  /**
   * Refine adaptively 5 times.
   */
  for (int step = 0; step < steps; step++) {
    /**
     * Refine a single grid point each time.
     * The SurplusRefinementFunctor chooses the grid point with the highest absolute surplus.
     * Refining the point means, that all children of this point (if not already present) are
     * added to the grid. Also all missing parents are added (recursively).
     * All new points are appended to the addedPoints vector.
     */
    base::SurplusRefinementFunctor functor(alpha, 1);
    grid->getGenerator().refine(functor, &addedPoints);

    /**
     * Extend alpha and funEval vector (new entries uninitialized). Note that right now, the size
     * of both vectors
     * matches number of gridpoints again, but the values of the new points are set to zero.
     */
    alpha.resize(gridStorage.getSize());
    funEvals.resize(gridStorage.getSize());

    /**
     * Evaluate the function f at the newly created gridpoints and set the
     * corresponding entries in the funEval vector with these values.
     */
    for (size_t i = 0; i < addedPoints.size(); i++) {
      size_t seq = addedPoints[i];
      base::GridPoint& gp = gridStorage.getPoint(seq);
      base::DataVector vec(dim);
      for (int j = 0; j < dim; ++j) {
        vec[j] = gp.getStandardCoordinate(j);
      }
      funEvals[seq] = numfunc(vec, ranges);
    }

    /**
     * Reset the alpha vector to function evals to prepare for hierarchisation.
     */
    alpha.copyFrom(funEvals);

    /**
     * Each time, we have to hierarchize the grid again, because in the previous interation,
     * new grid points have been added.
     */
    sgpp::op_factory::createOperationHierarchisation(*grid)->doHierarchisation(alpha);

    /**
     * Clear the addedPoints vector for the next iteration.
     */
    addedPoints.clear();

    // std::cout << "refinement step " << step + 1 << ", new grid size: " << alpha.getSize() <<
    // std::endl;
  }
  // for (size_t i = 0; i < alpha.getSize(); i++) {
  //  std::cout << alpha[i] << std::endl;
  //}

  /**
   * Calculate the surplus vector alpha for the interpolant of \f$
   * f(x)\f$.  Since the function can be evaluated at any
   * point. Hence. we simply evaluate it at the coordinates of the
   * grid points to obtain the nodal values. Then we use
   * hierarchization to obtain the surplus value.
   *
   */
  /**
    sgpp::base::DataVector alpha(gridStorage.getSize());
    double p[2];

    for (size_t i = 0; i < gridStorage.getSize(); i++) {
      sgpp::base::GridPoint& gp = gridStorage.getPoint(i);
      p[0] = gp.getStandardCoordinate(0);
      p[1] = gp.getStandardCoordinate(1);
      alpha[i] = numfunc(2, p, nullptr);
    }

    std::unique_ptr<base::OperationHierarchisation>(
        sgpp::op_factory::createOperationHierarchisation(*grid))
        ->doHierarchisation(alpha);
  */

  double prod = 1;
  for (auto pair : ranges) {
    prod *= (pair.second - pair.first);
  }
  // direct quadrature
  std::unique_ptr<sgpp::base::OperationQuadrature> opQ(
      sgpp::op_factory::createOperationQuadrature(*grid));
  double res = opQ->doQuadrature(alpha);
  return res * prod;
}

std::vector<std::vector<int>> PolynomialChaosExpansion::multiIndex(int dimension, int order) {
  std::vector<std::vector<int>> index;
  std::vector<int> curr(dimension);
  while (true) {
    index.push_back(curr);
    curr[0]++;
    for (int i = 0; i < dimension; ++i) {
      if (curr[dimension - 1] > order) {
        break;
      }
      if (curr[i] > order) {
        curr[i] -= order + 1;
        curr[i + 1]++;
      }
    }
    if (curr[dimension - 1] > order) {
      break;
    }
  }
  index.erase(std::remove_if(index.begin(), index.end(),
                             [order](std::vector<int> ee) {
                               return std::accumulate(ee.begin(), ee.end(), 0) > order;
                             }),
              index.end());

  return index;
}
base::DataVector PolynomialChaosExpansion::calculateCoefficients() {
  auto index = PolynomialChaosExpansion::multiIndex(static_cast<int>(types.size()), order);
  base::DataVector result(index.size());
  // calculate aj for each entry in the multiindex
  for (std::vector<std::vector<double>>::size_type j = 0; j < index.size(); ++j) {
    // lambdas for composite function to be integrated
    auto numfunc = [this, &index, j](const base::DataVector& vec) {
      double prd = 1;
      for (std::vector<double>::size_type i = 0; i < vec.getSize(); ++i) {
        prd *= evals[static_cast<int>(types[i])](index[j][i], vec[i]) *
               weights[static_cast<int>(types[i])](vec[i]);
      }
      return prd;
    };
    auto intfunc = [this, &numfunc](const base::DataVector& vec) {
      return numfunc(vec) * func(vec);
    };
    // integrate above function using mc
    // double num = monteCarloQuad(intfunc, 10000000);

    // integrate using sparsegrids
    // double num = sparseGridQuadrature(intfunc, static_cast<int>(types.size()), 15);

    // integrate using adaptive sparsegrids
    double num = adaptiveQuadrature(intfunc, static_cast<int>(types.size()), 8, 10);

    // calculate denominator
    double denom = 1.0;
    for (std::vector<distributionType>::size_type i = 0; i < types.size(); ++i) {
      denom *= denoms[static_cast<int>(types[i])](index[j][i]);
    }
    double aj = num / denom;
    result[j] = aj;
  }
  this->coefficients = result;
  return result;
}
base::DataVector PolynomialChaosExpansion::getCoefficients() { return coefficients; }

double PolynomialChaosExpansion::evalExpansion(const base::DataVector& xi) {
  if (coefficients.empty()) {
    calculateCoefficients();
  }
  auto index = PolynomialChaosExpansion::multiIndex(static_cast<int>(types.size()), order);
  double sum = 0.0;
  for (std::vector<std::vector<int>>::size_type j = 0; j < index.size(); ++j) {
    double prod = 1.0;
    for (std::vector<double>::size_type i = 0; i < index[j].size(); ++i) {
      prod *= evals[static_cast<int>(types[i])](index[j][i], xi[i]);
    }
    sum += prod * coefficients[j];
  }
  return sum;
}
// sample target and compare to pce eval
double PolynomialChaosExpansion::getL2Error() { return 0.0; }
}  // namespace datadriven
}  // namespace sgpp
