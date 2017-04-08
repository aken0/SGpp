/*
 * DBMatOfflineLU.hpp
 *
 *  Created on: 02.03.2017
 *      Author: michael
 */

#pragma once

#include <sgpp/datadriven/algorithm/DBMatOfflineGE.hpp>

namespace sgpp {
namespace datadriven {

class DBMatOfflineLU : public DBMatOfflineGE {
 public:
  DBMatOfflineLU(const DBMatDensityConfiguration& oc);

  void decomposeMatrix() override;

  void permuteVector(DataVector& b);

  /**
   * Store the decomposed matrix, the permutation and configuration.
   *
   * @param fname the file name
   */
  void store(const std::string& fname) override;

 private:
  std::unique_ptr<gsl_permutation> permutation;  // Stores the permutation that was
                                                 // applied on the matrix during decomposition
};

} /* namespace datadriven */
} /* namespace sgpp */
