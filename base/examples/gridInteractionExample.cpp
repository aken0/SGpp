// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

/**
 * \page example_gridInteractionExample_cpp Interaction-Term aware sparse grids.
 * This example shows how grids with more interaction terms differ from simpler grids.
 */

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/grid/generation/hashmap/HashGenerator.hpp>

#include <cstdlib>
#include <tuple>
#include <unordered_set>
#include <vector>

/**
 * @brief decodeCoords
 * @param coords are the coordinates of a grid point
 * @param result is the boolean version of the vector, each entry is true if the
 * corresponding dimension is used, i.e. not equal to 0.5
 */
void decodeCoords(sgpp::base::DataVector& coords, std::vector<bool>& result) {
  for (size_t i = 0; i < coords.getSize(); ++i) {
    result[i] = coords[i] != 0.5;
  }
}

/**
 * @brief main creates a grid with dimension three and level five and
 * adds more and more interaction points.
 * @details This is done by using interaction-term aware sparse grids.
 * First, we create a grid with interaction terms, that model only
 * one feature, then we add the pairwise interactions and finally
 * all interactions. The differences between the grids are printed.
 */
int main(int argc, char** argv) {
  auto dimensions = 3;
  auto level = 5;

  auto terms = std::unordered_set<std::vector<bool>>();
  auto realCoords = std::unordered_set<std::vector<bool>>();
  auto coords = sgpp::base::DataVector(dimensions);
  auto boolCoords = std::vector<bool>(dimensions);
  bool isInserted = false;

  // Add intercept.
  terms.insert(std::vector<bool>(dimensions, false));

  // Add all variables, without interaction
  for (auto dim = 0; dim < dimensions; ++dim) {
    auto vec = std::vector<bool>(dimensions, false);
    vec[dim] = true;
    terms.insert(vec);
  }

  {
    /**
    * Create our first grid.
    */
    auto grid = sgpp::base::Grid::createModLinearGrid(dimensions);
    auto& storage = grid->getStorage();
    auto generator = sgpp::base::HashGenerator();

    generator.regular_inter(storage, level, terms);

    std::cout << "Grid size with one level terms is " << grid->getSize() << std::endl;
    for (size_t i = 0; i < grid->getSize(); ++i) {
      sgpp::base::HashGridPoint gridIndex = storage.getPoint(i);
      gridIndex.getStandardCoordinates(coords);
      decodeCoords(coords, boolCoords);
      std::tie(std::ignore, isInserted) = realCoords.insert(boolCoords);
      if (isInserted) {
        std::cout << "New point" << coords.toString() << std::endl;
      }
    }
  }

  /**
   * Add all two-level-interactions
   */
  for (auto i = 0; i < dimensions; ++i) {
    for (auto j = 0; j < dimensions; ++j) {
      auto vec = std::vector<bool>(dimensions, false);
      vec[i] = true;
      vec[j] = true;
      terms.insert(vec);
    }
  }

  {
    auto grid = sgpp::base::Grid::createModLinearGrid(dimensions);
    auto& storage = grid->getStorage();
    auto generator = sgpp::base::HashGenerator();

    generator.regular_inter(storage, level, terms);

    std::cout << "Grid size with two level terms is " << grid->getSize() << std::endl;
    for (size_t i = 0; i < grid->getSize(); ++i) {
      sgpp::base::HashGridPoint gridIndex = storage.getPoint(i);
      gridIndex.getStandardCoordinates(coords);
      decodeCoords(coords, boolCoords);
      std::tie(std::ignore, isInserted) = realCoords.insert(boolCoords);
      if (isInserted) {
        std::cout << "New point" << coords.toString() << std::endl;
      }
    }
  }

  /**
   * Add all three-level-interactions
   */
  auto vec = std::vector<bool>(dimensions, true);
  terms.insert(vec);

  {
    auto grid = sgpp::base::Grid::createModLinearGrid(dimensions);
    auto& storage = grid->getStorage();
    auto generator = sgpp::base::HashGenerator();

    generator.regular_inter(storage, level, terms);

    std::cout << "Grid size with three level terms is " << grid->getSize() << std::endl;

    for (size_t i = 0; i < grid->getSize(); ++i) {
      sgpp::base::HashGridPoint gridIndex = storage.getPoint(i);
      gridIndex.getStandardCoordinates(coords);
      decodeCoords(coords, boolCoords);
      std::tie(std::ignore, isInserted) = realCoords.insert(boolCoords);
      if (isInserted) {
        std::cout << "New point" << coords.toString() << std::endl;
      }
    }
  }
}
