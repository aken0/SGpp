// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/operation/hash/OperationHierarchisationPolyBoundary.hpp>

#include <sgpp/base/operation/hash/common/algorithm_sweep/HierarchisationPolyBoundary.hpp>
#include <sgpp/base/operation/hash/common/algorithm_sweep/DehierarchisationPolyBoundary.hpp>

#include <sgpp/base/algorithm/sweep.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace base {

void OperationHierarchisationPolyBoundary::doHierarchisation(
  DataVector& node_values) {
  HierarchisationPolyBoundary func(storage, &base);
  sweep<HierarchisationPolyBoundary> s(func, storage);

  // Execute hierarchisation in every dimension of the grid
  for (size_t i = 0; i < this->storage.getDimension(); i++) {
    s.sweep1D_Boundary(node_values, node_values, i);
  }
}

void OperationHierarchisationPolyBoundary::doDehierarchisation(
  DataVector& alpha) {
  DehierarchisationPolyBoundary func(storage, &base);
  sweep<DehierarchisationPolyBoundary> s(func, storage);

  // Execute hierarchisation in every dimension of the grid
  for (size_t i = 0; i < this->storage.getDimension(); i++) {
    DataVector source(alpha);
    s.sweep1D_Boundary(source, alpha, i);
  }
}

}  // namespace base
}  // namespace sgpp
