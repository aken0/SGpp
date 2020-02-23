// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/combigrid/tools/LevelVectorTools.hpp>

#include <algorithm>
#include <cassert>
#include <numeric>
#include <vector>

namespace sgpp {
namespace combigrid {

std::vector<LevelVector> LevelVectorTools::generateHyperCube(
    const LevelVector& minLevel, const LevelVector& maxLevel) {
  LevelVector curLevel = minLevel;
  const size_t curDim = curLevel.size();
  std::vector<LevelVector> result;
  generateHyperCubeRecursive(minLevel, maxLevel, curLevel, curDim, result);
  return result;
}

void LevelVectorTools::generateHyperCubeRecursive(
    const LevelVector& minLevel, const LevelVector& maxLevel, LevelVector& curLevel,
    size_t curDim, std::vector<LevelVector>& result) {
  // if we are at the last dimension, add one one-dimensional row of the hypercube
  if (curDim == 1) {
    for (level_t l = minLevel[0]; l <= maxLevel[0]; ++l) {
      curLevel[0] = l;
      result.push_back(curLevel);
    }
  } else {
    // else, recurse to the next dimension
    for (level_t l = minLevel[curDim - 1]; l <= maxLevel[curDim - 1]; ++l) {
      curLevel[curDim - 1] = l;
      generateHyperCubeRecursive(minLevel, maxLevel, curLevel, curDim - 1, result);
    }
  }
}

std::vector<LevelVector> LevelVectorTools::generateDiagonal(
    const LevelVector& minLevel, level_t levelSum) {
  if (minLevel.empty()) {
    if (levelSum == 0) {
      return std::vector<LevelVector>{LevelVector{}};
    } else {
      return std::vector<LevelVector>{};
    }
  } else {
    const level_t minLevelSum = std::accumulate(minLevel.begin(), minLevel.end(), 0);
    std::vector<LevelVector> result;
    generateDiagonalRecursive(minLevel, minLevelSum, levelSum, LevelVector{}, result);
    return result;
  }
}

void LevelVectorTools::generateDiagonalRecursive(
    const LevelVector& minLevel, level_t minLevelSum, level_t levelSum, const LevelVector& suffix,
    std::vector<LevelVector>& result) {
  const size_t dim = minLevel.size();
  const size_t currentDim = dim - suffix.size();

  if (levelSum < minLevelSum) {
    return;
  } else if (currentDim == 1) {
    LevelVector l{levelSum};
    l.insert(l.end(), suffix.begin(), suffix.end());
    result.push_back(l);
  } else {
    const level_t newMinLevelSum = minLevelSum - minLevel[currentDim - 1];
    LevelVector newSuffix{0};
    newSuffix.insert(newSuffix.end(), suffix.begin(), suffix.end());

    for (level_t l = minLevel[currentDim - 1]; l <= levelSum - newMinLevelSum; l++) {
      newSuffix[0] = l;
      generateDiagonalRecursive(
          minLevel, newMinLevelSum, levelSum - l, newSuffix, result);
    }
  }
}

std::vector<LevelVector> LevelVectorTools::generateDiagonalWithBoundary(
    size_t dim, level_t levelSum) {
  return generateDiagonal(LevelVector(dim, 0), levelSum);
}

std::vector<LevelVector> LevelVectorTools::generateDiagonalWithoutBoundary(
    size_t dim, level_t levelSum) {
  return generateDiagonal(LevelVector(dim, 1), levelSum);
}

std::vector<LevelVector> LevelVectorTools::makeDownwardClosed(
    const std::vector<LevelVector>& subspaceLevels, LevelVector lowestLevelVector) {
  assert(lowestLevelVector.size() == subspaceLevels[0].size());
  std::vector<LevelVector> downwardClosedSet = subspaceLevels;

  // for each subspace level, ...
  for (const LevelVector& subspaceLevel : subspaceLevels) {
    // add the full hypercube of lower levels, if not already present
    for (const LevelVector& level : generateHyperCube(lowestLevelVector, subspaceLevel)) {
      if (std::find(downwardClosedSet.begin(), downwardClosedSet.end(), level) ==
          downwardClosedSet.end()) {
        downwardClosedSet.push_back(level);
      }
    }
  }

  std::sort(downwardClosedSet.begin(), downwardClosedSet.end());
  return downwardClosedSet;
}

}  // namespace combigrid
}  // namespace sgpp
