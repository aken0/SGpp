// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/combigrid/LevelIndexTypes.hpp>
#include <sgpp/combigrid/adaptive/PriorityEstimator.hpp>

#include <map>

namespace sgpp {
namespace combigrid {

/**
 * @brief the AveragingLevelManager from @holzmudd's combigrid module
 */
class AveragingPriorityEstimator : public PriorityEstimator {
 public:
  double estimatePriority(
      const LevelVector& levelVector,
      const std::map<LevelVector, double>& deltasOfDownwardNeighbors) const override;
};

}  // namespace combigrid
}  // namespace sgpp