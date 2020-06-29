// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once
// #ifdef USE_EIGEN

#include <algorithm>
#include <functional>
#include <map>
#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/grid/type/ModNakBsplineGrid.hpp>
#include <sgpp/base/operation/hash/common/basis/NakBsplineBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/NakBsplineBoundaryBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/NakBsplineExtendedBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/NakBsplineModifiedBasis.hpp>
#include <sgpp/datadriven/activeSubspaces/GaussQuadrature.hpp>
#include <sgpp/datadriven/activeSubspaces/MSplineBasis.hpp>
#include <tuple>
#include <utility>

namespace sgpp {
namespace datadriven {

/**
 * Calculates and stores the scalar products of not a knot B-spline functions with M-spline
 * functions
 */
class MSplineNakBsplineScalarProducts {
 public:
  /**
   * Constructor
   *
   * @param gridType          	type of the not a knot B-spline basis
   * @param degree				degree of the not a knot B-spline basis
   * @param quadOrder			order for the quadrature
   */
  MSplineNakBsplineScalarProducts(sgpp::base::GridType gridType, size_t degree,
                                  sgpp::base::DataVector xi, size_t quadOrder)
      : gridType(gridType), degree(degree), quadOrder(quadOrder) {
    // scale xi to [0,1], assuming it is already sorted
    double xiScale = xi.back() - xi[0];
    double xi0 = xi[0];
    for (size_t i = 0; i < xi.getSize(); i++) {
      xi[i] = (xi[i] - xi0) / xiScale;
    }
    this->xi = xi;
    mSplineBasis.setXi(this->xi);
    basis = initializeBasis(gridType, degree);
    base::GaussLegendreQuadRule1D gauss;
    gauss.getLevelPointsAndWeightsNormalized(quadOrder, coordinates, weights);
  }

  /**
   * initializes a nak B spline basis according to gidType and degree
   * @param gridType	type of the basis
   * @param degree		degree of the basis
   * @return 			the nak Bsplien basis
   */
  std::unique_ptr<sgpp::base::SBasis> initializeBasis(sgpp::base::GridType gridType, size_t degree);

  /**
   * Determines the common support of an M-spline (given by the private knots xi) and a B-spline of
   * given level and index
   *
   * @param level	level of the B-spline
   * @param index	index of the B-spline
   *
   * @return knots in the common support
   */
  sgpp::base::DataVector getCommonSupport(unsigned int level, unsigned int index);

  /**
   * calculates the one dimensional integral \int f*g dx where f and g are B-spline basis
   * functions or first derivatives of B-spline basis functions
   *
   * @param level 	level of the not a knot B-spline basis
   * @param index 	index of the not a knot B-spline basis
   * @param xi		knots defining the M-spline
   *
   * @return  scalarProduct of nak B-spline with given level and index with M-spline with given
   * knots xi
   */
  double basisScalarProduct(unsigned int level, unsigned int index);

  /**
   * calculates the scalar product of a sparse grid nak B-spline interpolant given through its grid
   * and coefficients and an M-spline given through its knots
   *
   * @param grid	grid of the nak B-spline interpolant
   * @param coeff	coefficients of the nak B-spline intepolant
   * @param xi		knots of the M-spline
   *
   * @return the scalar product \int \sum coeff1_i b_i(x) M(x) dx
   *
   */
  double calculateScalarProduct(std::shared_ptr<sgpp::base::Grid> grid,
                                sgpp::base::DataVector coeff);

 private:
  // type of not a knot B-spline basis
  sgpp::base::GridType gridType;
  // degrees of the not a knot B-spline basis
  size_t degree;
  // order for the quadrature
  size_t quadOrder;
  // knot sequence defining the M-spline
  sgpp::base::DataVector xi;
  // M-spline basis
  MSplineBasis mSplineBasis;
  // instance of the not a knot B-spline basis
  std::unique_ptr<sgpp::base::SBasis> basis;
  // quadrature coordinates
  sgpp::base::DataVector coordinates;
  // quadrature weights
  sgpp::base::DataVector weights;

  /**
   * used to get the support segments of a not a knot B-spline basis functions.
   *
   * @param level	level of the B-spline basis function
   * @param index	index of the B-spline basis function
   * @param degree	degree of the B-spline basis function
   *
   * @return the indices of the segments of the not a knot B-spline basis functions support
   */
  sgpp::base::DataVector nakBSplineSupport(size_t level, size_t index, size_t degree);
};

}  // namespace datadriven
}  // namespace sgpp

// #endif /* USE_EIGEN */
