// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>

#include <sgpp/base/operation/hash/common/algorithm_sweep/HierarchisationModLinear.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
namespace base {



HierarchisationModLinear::HierarchisationModLinear(GridStorage* storage) :
  storage(storage) {
}

HierarchisationModLinear::~HierarchisationModLinear() {
}

void HierarchisationModLinear::operator()(DataVector& source,
    DataVector& result, grid_iterator& index, size_t dim) {
  rec(source, result, index, dim, 0.0, 0.0);
}

void HierarchisationModLinear::rec(DataVector& source, DataVector& result,
                                   grid_iterator& index, size_t dim, float_t fl, float_t fr) {
  // current position on the grid
  size_t seq = index.seq();
  // value in the middle, needed for recursive call and calculation of the hierarchical surplus
  float_t fm = source[seq];

  GridStorage::index_type::level_type l;
  GridStorage::index_type::index_type i;

  index.get(dim, l, i);

  // recursive calls for the right and left side of the current node
  if (index.hint() == false) {
    float_t fltemp = fl;
    float_t frtemp = fr;

    // When we descend the hierarchical basis we have to modify the boundary values
    // in case the index is 1 or (2^l)-1 or we are on the first level
    // level 1, constant function
    if (l == 1) {
      // constant function
      fltemp = fm;
      frtemp = fm;
    }
    // left boundary
    else if (i == 1) {
      float_t ftemp;
      ftemp = fr - fm;
      fltemp = fm - ftemp;
    }
    // right boundary
    else if (static_cast<int>(i) == static_cast<int>((1 << l) - 1)) {
      float_t ftemp;
      ftemp = fl - fm;
      frtemp = fm - ftemp;
    }
    // inner functions
    else {
    }

    // descend left
    index.leftChild(dim);

    if (!storage->end(index.seq())) {
      rec(source, result, index, dim, fltemp, fm);
    }

    // descend right
    index.stepRight(dim);

    if (!storage->end(index.seq())) {
      rec(source, result, index, dim, fm, frtemp);
    }

    // ascend
    index.up(dim);
  }

  // hierarchisation
  result[seq] = fm - ((fl + fr) / 2.0);
}

// namespace detail

} // namespace SGPP
}