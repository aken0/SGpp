// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once
#include <functional>
#include <map>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/tools/OperationQuadratureMC.hpp>
#include <sgpp/globaldef.hpp>
#include <string>
#include <utility>
#include <vector>

namespace sgpp {
namespace datadriven {
enum class distributionType { Normal, Beta, Uniform, Exponential, Gamma };
class PolynomialChaosExpansion {
  std::function<double(const base::DataVector&)> func;
  int order;
  std::vector<distributionType> types;
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
  std::vector<std::function<double(double)>> weights;
  std::vector<std::function<double(double)>> denoms;
  std::vector<std::function<double(double, double)>> evals;

 public:
  PolynomialChaosExpansion(std::function<double(const base::DataVector&)> func, int order,
                           std::vector<distributionType> types,
                           std::vector<std::pair<double, double>> ranges, double alpha = 0.0,
                           double beta = 0.0);
  ~PolynomialChaosExpansion();

  std::vector<std::vector<int>> multiIndex(int dimension, int order);
  double monteCarloQuad(std::function<double(const base::DataVector&)> funct, size_t n);
  double sparseGridQuadrature(std::function<double(const base::DataVector&)> funct, int dim,
                              int n /*,int level*/);
  double adaptiveQuadrature(std::function<double(const base::DataVector&)> funct, int dim,
                            size_t n);
  double sparseGridQuadratureL2(std::function<double(const base::DataVector&)> funct, int dim,
                                int n /*,int level*/);
  double adaptiveQuadratureL2(std::function<double(const base::DataVector&)> funct, int dim,
                              size_t n);
  base::DataVector calculateCoefficients(int n, std::string method);
  base::DataVector getCoefficients();
  double evalExpansion(const base::DataVector& xi, int n, std::string method);
  double getL2Error();
  void printGrid(int dim, int level, std::string tFilename);
  void printAdaptiveGrid(std::function<double(const base::DataVector&)> funct, int dim, size_t n,
                         std::string tFilename);
};
}  // namespace datadriven
}  // namespace sgpp
