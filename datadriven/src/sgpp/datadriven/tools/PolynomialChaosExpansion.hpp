// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once
#include <functional>
#include <map>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/globaldef.hpp>
#include <string>
#include <vector>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
namespace sgpp {
namespace datadriven {
class PolynomialChaosExpansion {
  std::function<double(const base::DataVector&)> func;
  int order;
  std::vector<std::string> types;
  std::vector<std::pair<double, double>> ranges;
  double alpha;
  double beta;
  base::DataVector coefficients;

 private:
  double evalLegendre(int n, double x);
  double evalHermite(int n, double x);
  double evalLaguerre(int n, double x);
  double evalJacobi(int n, double x);
  double evalGenLaguerre(int n, double x);

 protected:
  std::map<std::string, std::function<double(double)>> weights;
  std::map<std::string, std::function<double(double)>> denoms;
  std::map<std::string, std::function<double(double, double)>> evals;

 public:
  PolynomialChaosExpansion(std::function<double(const base::DataVector&)> func, int order,
                           std::vector<std::string> types,
                           std::vector<std::pair<double, double>> ranges, double alpha = 0.0,
                           double beta = 0.0);
  ~PolynomialChaosExpansion();
  // move to private

  std::vector<std::vector<int>> multiIndex(int dimension, int order);
  double monteCarloQuad(std::function<double(const base::DataVector&)> func, long n);
  base::DataVector calculateCoefficients();
  base::DataVector getCoefficients();
  double evalExpansion();
};
}  // namespace datadriven
}  // namespace sgpp
