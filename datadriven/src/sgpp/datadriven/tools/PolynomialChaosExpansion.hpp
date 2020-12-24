// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once
#include <functional>
#include <map>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/tools/DistributionsVector.hpp>
//#include <sgpp/globaldef.hpp>
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
  std::vector<std::vector<int>> multiIndex(int dimension, int order);
  // multiple functions should be moved to private and are just here for debugging/plots
 public:
  /*
   *constructs a PCE using total-order expansion for the given function, expansion order and
   *underlying distributions
   */
  PolynomialChaosExpansion(std::function<double(const base::DataVector&)> func, int order,
                           sgpp::base::DistributionsVector distributions);

  ~PolynomialChaosExpansion();

  void printGrid(int dim, int n, std::string tFilename);

  void printAdaptiveGrid(std::function<double(const base::DataVector&)> funct, int dim, size_t n,
                         std::string tFilename);
  double monteCarloQuad(const std::function<double(const base::DataVector&)>& funct,
                        const size_t& n);
  double sparseGridQuadrature(const std::function<double(const base::DataVector&)>& funct, int dim,
                              int n, size_t quadOrder);
  double adaptiveQuadratureWeighted(const std::function<double(const base::DataVector&)>& funct,
                                    int dim, size_t n, size_t quadOrder);
  double sparseGridQuadratureL2(const std::function<double(const base::DataVector&)>& funct,
                                int dim, int n);
  double adaptiveQuadratureL2(const std::function<double(const base::DataVector&)>& funct, int dim,
                              size_t n);

  /*
   * calculates the coefficients using n points and the given method
   */
  base::DataVector calculateCoefficients(int n, std::string method);
  /*
   * returns the calculated coefficients
   */
  base::DataVector getCoefficients();

  void clearCoefficients();
  /*
   * evaluates the PCE at the given point
   */
  double evalExpansion(const base::DataVector& xi, int n, std::string method);
  /*
   * retrns the mean of the expansion
   */
  double getMean(int n, std::string method);
  /*
   * returns the variance of the expansion
   */
  double getVariance(int n, std::string method);
  /*
   * returns the L2 approximation error of the expansion
   */
  double getL2Error(int n, std::string method);
};
}  // namespace datadriven
}  // namespace sgpp
