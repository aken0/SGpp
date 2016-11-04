// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/operation/hash/OperationEval.hpp>
#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/operation/hash/common/basis/PolyModifiedClenshawCurtisBasis.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace base {

/**
 * This class implements OperationEval for a grids with poly basis ansatzfunctions with
 *
 */
class OperationEvalModPolyClenshawCurtis : public OperationEval {
 public:
  /**
   * Constructor
   *
   * @param storage the grid's GridStorage object
   * @param degree the polynom's max. degree
   */
  OperationEvalModPolyClenshawCurtis(GridStorage& storage, size_t degree)
      : storage(storage), base(degree) {}

  /**
   * Destructor
   */
  ~OperationEvalModPolyClenshawCurtis() override {}

  double eval(const DataVector& alpha, const DataVector& point) override;

 protected:
  /// Pointer to GridStorage object
  GridStorage& storage;
  /// PolyClenshawCurtis Basis object
  SPolyModifiedClenshawCurtisBase base;
};

}  // namespace base
}  // namespace sgpp
