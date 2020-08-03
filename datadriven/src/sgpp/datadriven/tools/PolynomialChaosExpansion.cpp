// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <functional>
#include <iostream>
#include <map>
#include <memory>
#include <numeric>
#include <random>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/tools/GridPrinter.hpp>
#include <sgpp/datadriven/tools/PolynomialChaosExpansion.hpp>
#include <sgpp/globaldef.hpp>
#include <sgpp_base.hpp>
#include <string>
#include <utility>
#include <vector>

#include "sgpp/base/tools/DistributionNormal.hpp"
namespace sgpp {
namespace datadriven {
PolynomialChaosExpansion::PolynomialChaosExpansion(std::function<double(const base::DataVector&)> f,
                                                   int order, std::vector<distributionType> types,
                                                   std::vector<std::pair<double, double>> ranges,
                                                   double alpha, double beta)
    : func(f), order(order), types(types), ranges(ranges), alpha(alpha), beta(beta) {
  // order:hermite,jacobi,legendre,laguerre,genlaguerre
  this->weights = std::vector<std::function<double(double)>>{
      {[](double x) { return std::exp(-std::pow(x, 2) / 2); }},
      {[this](double x) -> double {
        return std::pow((1.0 - x), this->alpha) * std::pow((1.0 + x), this->beta);
      }},
      {[](double x) { return 1.0; }},
      {[](double x) { return std::exp(-x); }},
      {[this](double x) -> double { return std::pow(x, this->alpha) * std::exp(-x); }}};
  this->denoms = std::vector<std::function<double(double)>>{
      {[](double j) { return std::sqrt(2.0 * M_PI) * std::tgamma(j + 1.0); }},
      {[this](double j) {
        return ((std::pow(2.0, this->alpha + this->beta + 1.0) /
                 (2.0 * j + this->alpha + this->beta + 1.0)) *
                ((std::tgamma(j + this->alpha + 1.0) * std::tgamma(j + this->beta + 1.0)) /
                 (std::tgamma(j + this->alpha + this->beta + 1.0) * std::tgamma(j + 1.0))));
      }},
      {[](double j) { return 2.0 / ((2.0 * j) + 1.0); }},
      {[](double j) { return 1.0; }},
      {[this](double j) { return std::tgamma(j + this->alpha + 1.0) / std::tgamma(j + 1.0); }}};

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
    double last = 1.0;
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
    return 1.0 - x;
  } else {
    double next = 0.0;
    double last = 1.0;
    double curr = 1.0 - x;
    for (double i = 1.0; i < n; ++i) {
      next = ((2.0 * i + 1.0 - x) * curr - i * last) / (i + 1.0);
      last = curr;
      curr = next;
    }
    return next;
  }
}

double PolynomialChaosExpansion::evalJacobi(int n, double x) {
  if (n == 0) {
    return 1.0;
  } else if (n == 1) {
    return (alpha + 1.0) + (alpha + beta + 2.0) * ((x - 1.0) / 2.0);
  } else {
    double next = 0.0;
    double last = 1.0;
    double curr = (alpha + 1.0) + (alpha + beta + 2.0) * ((x - 1.0) / 2.0);
    for (double i = 2.0; i <= n; ++i) {
      double q1 = ((2.0 * i + alpha + beta - 1.0) *
                   ((2.0 * i + alpha + beta) * (2.0 * i + alpha + beta - 2.0) * x +
                    std::pow(alpha, 2) - std::pow(beta, 2))) /
                  (2.0 * i * (i + alpha + beta) * (2.0 * i + alpha + beta - 2.0));
      double q2 = (2.0 * (i + alpha - 1.0) * (i + beta - 1.0) * (2.0 * i + alpha + beta)) /
                  (2.0 * i * (i + alpha + beta) * (2.0 * i + alpha + beta - 2.0));
      next = q1 * curr - q2 * last;
      last = curr;
      curr = next;
    }
    return next;
  }
}

double PolynomialChaosExpansion::evalGenLaguerre(int n, double x) {
  if (n == 0) {
    return 1.0;
  } else if (n == 1) {
    return 1.0 + alpha - x;
  } else {
    double next = 0.0;
    double last = 1.0;
    double curr = 1.0 + alpha - x;
    for (double i = 1.0; i < n; ++i) {
      next = ((2.0 * i + 1.0 + alpha - x) * curr - (i + alpha) * last) / (i + 1.0);
      last = curr;
      curr = next;
    }
    return next;
  }
}

double PolynomialChaosExpansion::monteCarloQuad(
    const std::function<double(const base::DataVector&)>& funct, const size_t& n) {
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
    base::DataVector randvec(dists.size());
    for (std::vector<std::uniform_real_distribution<double>>::size_type i = 0; i < dists.size();
         ++i) {
      randvec[i] = dists[i](mersenne);
    }
    return funct(base::DataVector(randvec));
  };
  base::DataVector results(n);
  std::generate(results.begin(), results.end(), gen);
  return factor * results.sum();
}

void PolynomialChaosExpansion::printGrid(int dim, int n, std::string tFilename) {
  std::unique_ptr<sgpp::base::Grid> grid(sgpp::base::Grid::createNakBsplineBoundaryGrid(dim, 3));
  sgpp::base::GridStorage& gridStorage = grid->getStorage();
  int i = 0;
  while (gridStorage.getSize() < n) {
    gridStorage.clear();
    grid->getGenerator().regular(i);
    ++i;
  }

  std::ofstream fileout;
  fileout.open(tFilename.c_str());
  for (size_t i = 0; i < grid->getSize(); ++i) {
    auto coords = gridStorage.getCoordinates(gridStorage.getPoint(i));
    for (size_t j = 0; j < grid->getDimension(); ++j) {
      fileout << coords.get(j);
      if (j < grid->getDimension() - 1) {
        fileout << ',';
      }
    }
    fileout << '\n';
  }
}

void PolynomialChaosExpansion::printAdaptiveGrid(
    std::function<double(const base::DataVector&)> funct, int dim, size_t n,
    std::string tFilename) {
  auto numfunc = [&funct](const base::DataVector& input,
                          std::vector<std::pair<double, double>>& ranges) {
    base::DataVector temp(input.size());
    for (base::DataVector::size_type i = 0; i < input.size(); ++i) {
      temp[i] = (input[i]) * (ranges[i].second - ranges[i].first) + ranges[i].first;
    }
    return funct(temp);
  };

  std::unique_ptr<sgpp::base::Grid> grid(sgpp::base::Grid::createNakBsplineBoundaryGrid(dim, 3));
  sgpp::base::GridStorage& gridStorage = grid->getStorage();
  grid->getGenerator().regular(1);

  base::DataVector coeffs(gridStorage.getSize());
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
   * Refine adaptively until number of points is reached.
   */
  while (gridStorage.getSize() < n) {
    /**
     * The SurplusRefinementFunctor chooses the grid point with the highest absolute surplus.
     * Refining the point means, that all children of this point (if not already present) are
     * added to the grid. Also all missing parents are added (recursively).
     */
    base::SurplusRefinementFunctor functor(coeffs, 10);
    grid->getGenerator().refine(functor, &addedPoints);

    coeffs.resize(gridStorage.getSize());
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

    coeffs.copyFrom(funEvals);

    // try hierarchisation
    bool succHierarch = false;

    try {
      std::unique_ptr<base::OperationHierarchisation>(
          sgpp::op_factory::createOperationHierarchisation(*grid))
          ->doHierarchisation(coeffs);
      succHierarch = true;
    } catch (...) {
      succHierarch = false;
    }

    if (!succHierarch) {
      std::cout << "Hierarchizing...\n\n";
      sgpp::base::HierarchisationSLE hierSLE(*grid);
      sgpp::base::sle_solver::Eigen sleSolver;

      // solve linear system
      if (!sleSolver.solve(hierSLE, funEvals, coeffs)) {
        std::cout << "Solving failed, exiting.\n";
        return;
      }
    }
    addedPoints.clear();
  }
  std::ofstream fileout;
  fileout.open(tFilename.c_str());
  for (size_t i = 0; i < grid->getSize(); ++i) {
    auto coords = gridStorage.getCoordinates(gridStorage.getPoint(i));
    for (size_t j = 0; j < grid->getDimension(); ++j) {
      fileout << coords.get(j);
      if (j < grid->getDimension() - 1) {
        fileout << ',';
      }
    }
    fileout << '\n';
  }
}

double PolynomialChaosExpansion::sparseGridQuadrature(
    const std::function<double(const base::DataVector&)>& funct, int dim, int n /*,int level*/) {
  auto numfunc = [&funct](const base::DataVector& input,
                          const std::vector<std::pair<double, double>>& ranges) {
    base::DataVector temp(input.size());
    for (base::DataVector::size_type i = 0; i < input.size(); ++i) {
      temp[i] = input[i] * (ranges[i].second - ranges[i].first) + ranges[i].first;
    }
    return funct(temp);
  };
  std::unique_ptr<sgpp::base::Grid> grid(sgpp::base::Grid::createNakBsplineBoundaryGrid(dim, 3));
  sgpp::base::GridStorage& gridStorage = grid->getStorage();
  int i = 0;
  while (gridStorage.getSize() < n) {
    gridStorage.clear();
    grid->getGenerator().regular(i);
    ++i;
  }
  /**
   * Calculate the surplus vector alpha for the interpolant of \f$
   * f(x)\f$.  Since the function can be evaluated at any
   * point. Hence. we simply evaluate it at the coordinates of the
   * grid points to obtain the nodal values. Then we use
   * hierarchization to obtain the surplus value.
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
  bool succHierarch = false;

  try {
    std::unique_ptr<base::OperationHierarchisation>(
        sgpp::op_factory::createOperationHierarchisation(*grid))
        ->doHierarchisation(evals);
    succHierarch = true;
  } catch (...) {
    succHierarch = false;
  }

  base::DataVector coeffs(evals.getSize());
  if (!succHierarch) {
    std::cout << "Hierarchizing...\n\n";
    sgpp::base::HierarchisationSLE hierSLE(*grid);
    sgpp::base::sle_solver::Eigen sleSolver;

    // solve linear system
    if (!sleSolver.solve(hierSLE, evals, coeffs)) {
      std::cout << "Solving failed, exiting.\n";
      return 1;
    }
    // overwrite evals
    evals = coeffs;
  }

  std::unique_ptr<sgpp::base::OperationQuadrature> opQ(
      sgpp::op_factory::createOperationQuadrature(*grid));
  double res = opQ->doQuadrature(evals);

  double prod = 1;
  for (auto pair : ranges) {
    prod *= (pair.second - pair.first);
  }
  return res * prod;
}

double PolynomialChaosExpansion::adaptiveQuadrature(
    const std::function<double(const base::DataVector&)>& funct, int dim, size_t n) {
  auto numfunc = [&funct](const base::DataVector& input,
                          const std::vector<std::pair<double, double>>& ranges) {
    base::DataVector temp(input.size());
    for (base::DataVector::size_type i = 0; i < input.size(); ++i) {
      temp[i] = (input[i]) * (ranges[i].second - ranges[i].first) + ranges[i].first;
    }
    return funct(temp);
  };

  std::unique_ptr<sgpp::base::Grid> grid(sgpp::base::Grid::createNakBsplineBoundaryGrid(dim, 3));
  sgpp::base::GridStorage& gridStorage = grid->getStorage();
  grid->getGenerator().regular(1);

  base::DataVector coeffs(gridStorage.getSize());
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
  std::vector<size_t> addedPoints;
  /**
   * Refine adaptively until number of points is reached.
   */
  while (gridStorage.getSize() < n) {
    /**
     * The SurplusRefinementFunctor chooses the grid point with the highest absolute surplus.
     * Refining the point means, that all children of this point (if not already present) are
     * added to the grid. Also all missing parents are added (recursively).
     */
    base::SurplusRefinementFunctor functor(coeffs, 10);
    grid->getGenerator().refine(functor, &addedPoints);

    coeffs.resize(gridStorage.getSize());
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

    coeffs.copyFrom(funEvals);

    // try hierarchisation
    bool succHierarch = false;

    try {
      std::unique_ptr<base::OperationHierarchisation>(
          sgpp::op_factory::createOperationHierarchisation(*grid))
          ->doHierarchisation(coeffs);
      succHierarch = true;
    } catch (...) {
      succHierarch = false;
    }

    if (!succHierarch) {
      std::cout << "Hierarchizing...\n\n";
      sgpp::base::HierarchisationSLE hierSLE(*grid);
      sgpp::base::sle_solver::Eigen sleSolver;

      // solve linear system
      if (!sleSolver.solve(hierSLE, funEvals, coeffs)) {
        std::cout << "Solving failed, exiting.\n";
        return 1;
      }
    }
    addedPoints.clear();
  }
  double prod = 1;
  for (auto pair : ranges) {
    prod *= (pair.second - pair.first);
  }
  std::unique_ptr<sgpp::base::OperationQuadrature> opQ(
      sgpp::op_factory::createOperationQuadrature(*grid));
  double res = opQ->doQuadrature(coeffs);
  return res * prod;
}

double PolynomialChaosExpansion::sparseGridQuadratureL2(
    const std::function<double(const base::DataVector&)>& funct, int dim, int n /*,int level*/) {
  auto numfunc = [&funct](const base::DataVector& input,
                          const std::vector<std::pair<double, double>>& ranges) {
    base::DataVector temp(input.size());
    for (base::DataVector::size_type i = 0; i < input.size(); ++i) {
      temp[i] = input[i] * (ranges[i].second - ranges[i].first) + ranges[i].first;
    }
    return funct(temp);
  };
  std::unique_ptr<sgpp::base::Grid> grid(sgpp::base::Grid::createNakBsplineBoundaryGrid(dim, 3));
  sgpp::base::GridStorage& gridStorage = grid->getStorage();

  int i = 0;
  while (gridStorage.getSize() < n) {
    gridStorage.clear();
    grid->getGenerator().regular(i);
    ++i;
  }

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
  bool succHierarch = false;

  try {
    std::unique_ptr<base::OperationHierarchisation>(
        sgpp::op_factory::createOperationHierarchisation(*grid))
        ->doHierarchisation(evals);
    succHierarch = true;
  } catch (...) {
    succHierarch = false;
  }

  base::DataVector coeffs(evals.getSize());
  if (!succHierarch) {
    std::cout << "Hierarchizing...\n\n";
    sgpp::base::HierarchisationSLE hierSLE(*grid);
    sgpp::base::sle_solver::Eigen sleSolver;

    // solve linear system
    if (!sleSolver.solve(hierSLE, evals, coeffs)) {
      std::cout << "Solving failed, exiting.\n";
      return 1;
    }
    evals = coeffs;
  }
  // sample grid at MC points, compare opEval to transformed function(numfunc)
  std::vector<std::uniform_real_distribution<double>> dists(dim);
  std::random_device dev;
  std::mt19937_64 mersenne{dev()};
  for (int i = 0; i < dim; ++i) {
    dists[i] = std::uniform_real_distribution<double>{0.0, 1.0};
  }
  auto gen = [&numfunc, &dists, &mersenne, &grid, &evals, this]() {
    base::DataVector randvec(dists.size());
    for (std::vector<std::uniform_real_distribution<double>>::size_type i = 0; i < dists.size();
         ++i) {
      randvec[i] = dists[i](mersenne);
    }
    std::unique_ptr<sgpp::base::OperationEval> opEval(
        sgpp::op_factory::createOperationEvalNaive(*grid));
    return std::pow(numfunc(base::DataVector(randvec), ranges) - (opEval->eval(evals, randvec)), 2);
  };
  size_t num = 100000;
  base::DataVector results(num);
  std::generate(results.begin(), results.end(), gen);
  return results.sum() / static_cast<double>(num);
}
double PolynomialChaosExpansion::adaptiveQuadratureL2(
    const std::function<double(const base::DataVector&)>& funct, int dim, size_t n) {
  auto numfunc = [&funct](const base::DataVector& input,
                          const std::vector<std::pair<double, double>>& ranges) {
    base::DataVector temp(input.size());
    for (base::DataVector::size_type i = 0; i < input.size(); ++i) {
      temp[i] = (input[i]) * (ranges[i].second - ranges[i].first) + ranges[i].first;
    }
    return funct(temp);
  };

  std::unique_ptr<sgpp::base::Grid> grid(sgpp::base::Grid::createNakBsplineBoundaryGrid(dim, 3));
  sgpp::base::GridStorage& gridStorage = grid->getStorage();
  grid->getGenerator().regular(1);

  base::DataVector coeffs(gridStorage.getSize());

  base::DataVector funEvals(gridStorage.getSize());
  for (size_t i = 0; i < gridStorage.getSize(); i++) {
    base::GridPoint& gp = gridStorage.getPoint(i);
    base::DataVector vec(dim);
    for (int j = 0; j < dim; ++j) {
      vec[j] = gp.getStandardCoordinate(j);
    }
    funEvals[i] = numfunc(vec, ranges);
  }

  std::vector<size_t> addedPoints;

  /**
   * Refine adaptively until number of points is reached.
   */
  while (gridStorage.getSize() < n) {
    /**
     * The SurplusRefinementFunctor chooses the grid point with the highest absolute surplus.
     * Refining the point means, that all children of this point (if not already present) are
     * added to the grid. Also all missing parents are added (recursively).
     */
    base::SurplusRefinementFunctor functor(coeffs, 10);
    grid->getGenerator().refine(functor, &addedPoints);

    coeffs.resize(gridStorage.getSize());
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

    coeffs.copyFrom(funEvals);

    // try hierarchisation
    bool succHierarch = false;

    try {
      std::unique_ptr<base::OperationHierarchisation>(
          sgpp::op_factory::createOperationHierarchisation(*grid))
          ->doHierarchisation(coeffs);
      succHierarch = true;
    } catch (...) {
      succHierarch = false;
    }

    if (!succHierarch) {
      std::cout << "Hierarchizing...\n\n";
      sgpp::base::HierarchisationSLE hierSLE(*grid);
      sgpp::base::sle_solver::Eigen sleSolver;

      // solve linear system
      if (!sleSolver.solve(hierSLE, funEvals, coeffs)) {
        std::cout << "Solving failed, exiting.\n";
        return 1;
      }
    }
    addedPoints.clear();
  }
  // sample grid at MC points, compare opEval to transformed function(numfunc)
  std::vector<std::uniform_real_distribution<double>> dists(dim);
  std::random_device dev;
  std::mt19937_64 mersenne{dev()};
  for (int i = 0; i < dim; ++i) {
    dists[i] = std::uniform_real_distribution<double>{0.0, 1.0};
  }

  auto gen = [&numfunc, &dists, &mersenne, &grid, &coeffs, this]() {
    base::DataVector randvec(dists.size());
    for (std::vector<std::uniform_real_distribution<double>>::size_type i = 0; i < dists.size();
         ++i) {
      randvec[i] = dists[i](mersenne);
    }
    std::unique_ptr<sgpp::base::OperationEval> opEval(
        sgpp::op_factory::createOperationEvalNaive(*grid));
    return std::pow(numfunc(base::DataVector(randvec), ranges) - (opEval->eval(coeffs, randvec)),
                    2);
  };
  size_t num = 100000;
  base::DataVector results(num);
  std::generate(results.begin(), results.end(), gen);
  return results.sum() / static_cast<double>(num);
}

std::vector<std::vector<int>> PolynomialChaosExpansion::multiIndex(int dimension, int order) {
  std::vector<std::vector<int>> index(static_cast<int>(std::pow(order + 1, dimension)));
  std::vector<int> curr(dimension);
  for (size_t j = 0; j < index.size(); ++j) {
    index[j] = curr;
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
base::DataVector PolynomialChaosExpansion::calculateCoefficients(int n, std::string method) {
  auto index = multiIndex(static_cast<int>(types.size()), order);
  base::DataVector result(index.size());
  // calculate aj for each entry in the multiIndex
  for (std::vector<std::vector<int>>::size_type j = 0; j < index.size(); ++j) {
    auto numfunc = [this, &index, &j](const base::DataVector& vec) {
      double prd = 1;
      for (base::DataVector::size_type i = 0; i < vec.getSize(); ++i) {
        prd *= evals[static_cast<int>(types[i])](index[j][i], vec[i]) *
               weights[static_cast<int>(types[i])](vec[i]);
      }
      return prd;
    };
    auto intfunc = [this, &numfunc](const base::DataVector& vec) {
      return numfunc(vec) * func(vec);
    };
    double num;
    if (method == "MC") {
      num = monteCarloQuad(intfunc, n);
    } else if (method == "sparseGrid") {
      num = sparseGridQuadrature(intfunc, static_cast<int>(types.size()), n);
    } else if (method == "adaptiveGrid") {
      num = adaptiveQuadrature(intfunc, static_cast<int>(types.size()), n);
    }
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

void PolynomialChaosExpansion::clearCoefficients() { this->coefficients.clear(); }

double PolynomialChaosExpansion::evalExpansion(const base::DataVector& xi, int n,
                                               std::string method) {
  if (coefficients.empty()) {
    calculateCoefficients(n, method);
  }
  auto index = multiIndex(static_cast<int>(types.size()), order);
  double sum = 0.0;
  for (std::vector<std::vector<int>>::size_type j = 0; j < index.size(); ++j) {
    double prod = 1.0;
    for (std::vector<int>::size_type i = 0; i < index[j].size(); ++i) {
      prod *= evals[static_cast<int>(types[i])](index[j][i], xi[i]);
    }
    sum += prod * coefficients[j];
  }
  return sum;
}

// sample response and compare to pce eval
double PolynomialChaosExpansion::getL2Error(int n, std::string method) {
  int dim = static_cast<int>(types.size());
  std::vector<std::uniform_real_distribution<double>> dists(dim);
  std::random_device dev;
  std::mt19937_64 mersenne{dev()};
  for (int i = 0; i < dim; ++i) {
    dists[i] = std::uniform_real_distribution<double>{ranges[i].first, ranges[i].second};
  }
  auto gen = [this, &dists, &mersenne, &n, &method]() {
    base::DataVector randvec(dists.size());
    for (std::vector<std::uniform_real_distribution<double>>::size_type i = 0; i < dists.size();
         ++i) {
      randvec[i] = dists[i](mersenne);
    }
    return std::pow(func(randvec) - (evalExpansion(randvec, n, method)), 2);
  };
  size_t num = 100000;
  base::DataVector results(num);
  std::generate(results.begin(), results.end(), gen);
  return results.sum() / static_cast<double>(num);
}
double PolynomialChaosExpansion::getMean(int n, std::string method) {
  if (coefficients.empty()) {
    calculateCoefficients(n, method);
  }
  return coefficients[0];
}
double PolynomialChaosExpansion::getVariance(int n, std::string method) {
  if (coefficients.empty()) {
    calculateCoefficients(n, method);
  }
  base::DataVector temp(coefficients.getSize());
  temp.copyFrom(coefficients);
  temp.set(0, 0.0);
  temp.sqr();
  return temp.sum();
}
}  // namespace datadriven
}  // namespace sgpp
