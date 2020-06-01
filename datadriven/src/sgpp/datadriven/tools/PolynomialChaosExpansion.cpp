// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org
#include <algorithm>
#include <cmath>
#include <functional>
#include <map>
#include <numeric>
#include <random>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/datadriven/tools/PolynomialChaosExpansion.hpp>
#include <sgpp/globaldef.hpp>
#include <string>
#include <utility>
#include <vector>
namespace sgpp {
namespace datadriven {
PolynomialChaosExpansion::PolynomialChaosExpansion(std::function<double(const base::DataVector&)> f,
                                                   int order, std::vector<std::string> types,
                                                   std::vector<std::pair<double, double>> ranges,
                                                   double alpha, double beta)
    : func(f), order(order), types(types), ranges(ranges), alpha(alpha), beta(beta) {
  this->weights = std::map<std::string, std::function<double(double)>>{
      {"hermite", [](double x) { return std::exp(-std::pow(x, 2) / 2); }},
      {"jacobi",
       [this](double x) -> double {
         return std::pow((1 - x), this->alpha) * std::pow((1 + x), this->beta);
       }},
      {"legendre", [](double x) { return 1.0; }},
      {"laguerre", [](double x) { return std::exp(-x); }},
      {"genlaguerre",
       [this](double x) -> double { return std::pow(x, this->alpha) * std::exp(-x); }}};
  this->denoms = std::map<std::string, std::function<double(double)>>{
      {"hermite", [](double j) { return std::sqrt(2 * M_PI) * std::tgamma(j + 1); }},
      {"jacobi",
       [this](double j) {
         return (
             (std::pow(2, this->alpha + this->beta + 1) / (2 * j + this->alpha + this->beta + 1)) *
             ((std::tgamma(j + this->alpha + 1) * std::tgamma(j + this->beta + 1)) /
              (std::tgamma(j + this->alpha + this->beta + 1) * std::tgamma(j + 1))));
       }},
      {"legendre", [](double j) { return 2 / ((2 * j) + 1); }},
      {"laguerre", [](double j) { return 1.0; }},
      {"genlaguerre",
       [this](double j) { return std::tgamma(j + this->alpha + 1) / std::tgamma(j + 1); }}};
  this->evals = std::map<std::string, std::function<double(double, double)>>{
      {"hermite", [this](int n, double x) { return evalHermite(n, x); }},
      {"jacobi", [this](int n, double x) { return evalJacobi(n, x); }},
      {"legendre", [this](int n, double x) { return evalLegendre(n, x); }},
      {"laguerre", [this](int n, double x) { return evalLaguerre(n, x); }},
      {"genlaguerre", [this](int n, double x) { return evalGenLaguerre(n, x); }}};
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
// todo?
double PolynomialChaosExpansion::monteCarloQuad(
    std::function<double(const base::DataVector&)> funct, long n) {
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
  auto gen = [funct, &dists, &mersenne]() {
    std::vector<double> randvec(dists.size());
    for (std::vector<std::uniform_real_distribution<double>>::size_type j = 0; j < dists.size();
         ++j) {
      randvec[j] = dists[j](mersenne);
    }
    return funct(base::DataVector(randvec));
  };
  std::vector<double> results(n);
  std::generate(results.begin(), results.end(), gen);
  double sum = std::accumulate(results.begin(), results.end(), 0.0);
  return factor * sum;
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
  base::DataVector result;
  auto index = PolynomialChaosExpansion::multiIndex(static_cast<int>(types.size()), order);
  // calculate aj for each entry in the multiindex
  for (auto entry : index) {
    // lambdas for composite function to be integrated
    auto numfunc = [this, entry](const base::DataVector& vec) {
      double prd = 1;
      for (std::vector<double>::size_type i = 0; i < vec.getSize(); ++i) {
        prd *= evals[types[i]](entry[i], vec[i]) * weights[types[i]](vec[i]);
      }
      return prd;
    };
    auto intfunc = [this, numfunc](const base::DataVector& vec) {
      return numfunc(vec) * func(vec);
    };
    // integrate above function using mc
    double num = monteCarloQuad(intfunc, 10000000);
    // calculate denominator
    double denom = 1.0;
    for (std::vector<std::string>::size_type i = 0; i < types.size(); ++i) {
      denom *= denoms[types[i]](entry[i]);
    }
    double aj = num / denom;
    result.push_back(aj);
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
      prod *= evals[types[i]](index[j][i], xi[i]);
    }
    sum += prod * coefficients[j];
  }
  return sum;
}
}  // namespace datadriven
}  // namespace sgpp
