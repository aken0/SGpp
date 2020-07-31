// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <cmath>
#include <iostream>
#include <numeric>
#include <random>
#include <sgpp/base/tools/Distribution.hpp>
#include <string>

namespace sgpp {
namespace base {

/**
 */
class DistributionBeta : public Distribution {
 public:
  /**
   * Constructor
   */
  explicit DistributionBeta(double alpha = 0, double beta = 0)
      : Distribution(), alpha(alpha), beta(beta), dist(alpha, beta) {}

  /**
   * Destructor
   */
  virtual ~DistributionBeta() {}

  /**
   *
   */
  double sample() { return dist(gen); }

  double eval(double x) {
    return (std::pow(1 - x, alpha) * std::pow(1 + x, beta)) /
           (std::pow(2, alpha + beta + 1) * (std::tgamma(alpha) * std::tgamma(beta)) /
            std::tgamma(alpha + beta));
  }

  sgpp::base::DataVector getBounds() {
    sgpp::base::DataVector bounds{-1, 1};
    return bounds;
  }

  sgpp::base::DistributionType getType() { return sgpp::base::DistributionType::Beta; }

  sgpp::base::DataVector getCharacteristics() { return sgpp::base::DataVector{alpha, beta}; }

 private:
  double alpha;
  double beta;
  std::gamma_distribution<double> dist;
  std::gamma_distribution<double> dist2;
};
}  // namespace base
}  // namespace sgpp
