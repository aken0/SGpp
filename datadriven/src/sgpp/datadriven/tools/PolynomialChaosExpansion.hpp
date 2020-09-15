// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once
#include <functional>
#include <map>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/tools/DistributionsVector.hpp>
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
  sgpp::base::DistributionsVector distributions;
  sgpp::base::DistributionsVector standardvec;
  std::vector<std::pair<double, double>> ranges;
  base::DataVector alpha;
  base::DataVector beta;
  base::DataVector coefficients;

 private:
  double evalLegendre(int n, double x);
  double evalHermite(int n, double x);
  double evalLaguerre(int n, double x);
  double evalJacobi(int n, double x, size_t i);
  double evalGenLaguerre(int n, double x, size_t i);
  std::vector<std::function<double(double, size_t)>> weights;
  std::vector<std::function<double(double, size_t)>> denoms;
  std::vector<std::function<double(double, double, size_t)>> evals;
  // multiple functions should be moved to private and are just here for debugging
 public:
  PolynomialChaosExpansion(std::function<double(const base::DataVector&)> func, int order,
                           sgpp::base::DistributionsVector distributions);
  ~PolynomialChaosExpansion();

  std::vector<std::vector<int>> multiIndex(int dimension, int order);
  void printGrid(int dim, int n, std::string tFilename);

  void printAdaptiveGrid(std::function<double(const base::DataVector&)> funct, int dim, size_t n,
                         std::string tFilename);
  double monteCarloQuad(const std::function<double(const base::DataVector&)>& funct,
                        const size_t& n);
  double sparseGridQuadrature(const std::function<double(const base::DataVector&)>& funct, int dim,
                              int n /*,int level*/);
  double adaptiveQuadrature(const std::function<double(const base::DataVector&)>& funct, int dim,
                            size_t n);
  double sparseGridQuadratureL2(const std::function<double(const base::DataVector&)>& funct,
                                int dim, int n /*,int level*/);
  double adaptiveQuadratureL2(const std::function<double(const base::DataVector&)>& funct, int dim,
                              size_t n);
  base::DataVector calculateCoefficients(int n, std::string method);
  base::DataVector getCoefficients();
  void clearCoefficients();
  double evalExpansion(const base::DataVector& xi, int n, std::string method);
  double getMean(int n, std::string method);
  double getVariance(int n, std::string method);
  double getL2Error(int n, std::string method);
};
}  // namespace datadriven
}  // namespace sgpp
