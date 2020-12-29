// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once
// #ifdef USE_EIGEN

#include <sgpp/base/function/scalar/InterpolantScalarFunction.hpp>
#include <sgpp/base/function/scalar/InterpolantScalarFunctionGradient.hpp>
#include <sgpp/base/function/scalar/WrapperScalarFunctionGradient.hpp>
#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/grid/generation/functors/SurplusRefinementFunctor.hpp>
#include <sgpp/base/grid/type/NakBsplineExtendedGrid.hpp>
#include <sgpp/base/grid/type/NakBsplineGrid.hpp>
#include <sgpp/base/operation/hash/common/basis/NakBsplineBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/NakBsplineBoundaryBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/NakBsplineExtendedBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/NakBsplineModifiedBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/NakBsplineExtendedBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/NakPBsplineBasis.hpp>
#include <sgpp/base/tools/RandomNumberGenerator.hpp>
#include <sgpp/base/tools/sle/solver/Armadillo.hpp>
#include <sgpp/base/tools/sle/system/HierarchisationSLE.hpp>
#include <sgpp/datadriven/activeSubspaces/ASMatrix.hpp>
#include <sgpp/datadriven/activeSubspaces/ASResponseSurfaceNakBspline.hpp>
#include <sgpp/datadriven/activeSubspaces/GaussQuadrature.hpp>
#include <sgpp/datadriven/activeSubspaces/NakBsplineScalarProducts.hpp>

#include <string>

namespace sgpp {
namespace datadriven {

/**
 * Used to create, store and use the matrix C for the detection of active subspaces using a B-spline
 * interpolant. So C_{i,j} = \int \nabla f \nabla f^T dx \approx \int \nabla \hat{f} \nabla
 * \hat{f}^T dx where \hat{f} is a nak B-spline interpolant for f
 */
class ASMatrixBspline : public ASMatrix {
 public:
  /**
   * Constructor
   *
   * @param numDim			  number of dimensions
   * @param degree            degree for the B-spline basis functions
   * @param gridType          type of the grid for the interpolant
   * @param evaluationPoints  evaluationPoints only needed for ASMAtrixBsplineData
   * @param functionValues 	  functionValues only needed for ASMAtrixBsplineData
   */
  ASMatrixBspline(size_t numDim, size_t degree, sgpp::base::GridType gridType,
                  sgpp::base::DataMatrix evaluationPoints = sgpp::base::DataMatrix(),
                  sgpp::base::DataVector functionValues = sgpp::base::DataVector())
      : ASMatrix(evaluationPoints, functionValues), numDim(numDim), degree(degree) {
    initialize(gridType);
  }

  /**
   * Destructor
   */
  virtual ~ASMatrixBspline() {}

  /**
   * Sets grid and basis according to grid type
   *
   * @param gridType	type for grid and basis
   */
  void initialize(sgpp::base::GridType gridType);

  /**
   * Create a regular interpolant of the objective function f
   *
   * @param level	level of the underlying grid
   */
  virtual void buildRegularInterpolant(size_t level) = 0;

  /**
   * Create a spatially adaptive interpolant of the objective function f
   *
   * @param maxNumGridPoints	upper threshold for the number of grid points
   * @param initialLevel		the refinement needs an initial regular grid of initialLevel
   * @param refinementsNum		maximum number of points refined in one step
   */
  virtual void buildAdaptiveInterpolant(size_t maxNumGridPoints, size_t initialLevel,
                                        size_t refinementsNum) = 0;

  /**
   * calculates the coefficients for the interpolant based on the objective function and grid.
   * Must be called after every change to the grid!
   */
  virtual void calculateCoefficients() = 0;

  /**
   * General routine to create the Matrix C, currently simply wraps createMatrixMonteCarlo.
   * createMatrixGauss() is better and should be preferred
   */
  void createMatrix(size_t numPoints);

  /**
   * creates the matrix C from the interpolant by using Monte Carlo quadrature
   *
   * @param numPoints	number of Monte Carlo points
   */
  void createMatrixMonteCarlo(size_t numMCPoints);

  /**
   * creates the matrix C using the exact integrals of the underlying B-spline functions from a
   * Gauss quadrature of sufficient degree
   */
  void createMatrixGauss();

  // ----------------- auxiliary routines -----------

  /**
   * refine the current interpolant surplus adaptive by ading the children of <= refinementsNum
   * points
   *
   * @param refinementsNum maximum number of points to be refined
   */
  void refineSurplusAdaptive(size_t refinementsNum);

  /**
   * calculates the entry C_{i,j} of the matrix,
   * C_{i,j} = int df/dxi df/dxj dx
   *
   * @param i 				row
   * @param j 				column
   * @param scalarProducts	a NakBsplineScalarProducts instance to calculate the scalar products
   * @return matrix entry C_{i,j}
   */
  double matrixEntryGauss(size_t i, size_t j,
                          sgpp::datadriven::NakBsplineScalarProducts& scalarProducts);

  /**
   * calculates \int d/dx_i b_k(x) d/dx_j b_l(x) dx
   * This is done by evaluating one dimensional integrals of b_{k_d} b_{l_d},
   * d/dx_i b_{k_i} b_{l_i}, b_{k_j} d/dx_j b_{l_j}, d/dx_i b_{k_i} d/dx_i b_{l_i}
   *
   *@param i index of matrix row
   *@param j index of matrix column
   *@param k index of first basis function
   *@param l index of second basis function
   * @param scalarProducts	a NakBsplineScalarProducts instance to calculate the scalar products
   * @return integral \int d/dx_i b_k(x) d/dx_j b_l(x) dx
   */
  double scalarProductDxbiDxbj(size_t i, size_t j, size_t k, size_t l,
                               sgpp::datadriven::NakBsplineScalarProducts& scalarProducts);

  /**
   * returns the interpolation coefficients
   *
   * @return interpolation coefficients
   */
  sgpp::base::DataVector getCoefficients() { return coefficients; }

  /**
   * save the grid and the inteprolation coefficients to file
   *
   * @param path
   *
   */
  void toFile(std::string path);

  /**
   * creates and returns an ASResponseSurfaceNakBspline based on the active subspace calculated with
   * this class
   *
   * @param asDimension dimension of the active subspace
   * @param gridType    underlying grid type for the response surface
   * @param degree      degree of the underlying B-spline basis
   *
   * @return a response surface based on one of the not a knot B-spine basis types exploiting the
   * active subspace
   */
  sgpp::datadriven::ASResponseSurfaceNakBspline getResponseSurfaceInstance(
      size_t asDimension, sgpp::base::GridType gridType, size_t degree = 3);

  double evalInterpolant(sgpp::base::DataVector);

 protected:
  size_t numDim;
  sgpp::base::GridType gridType;
  size_t degree;
  sgpp::base::DataVector coefficients;
  std::shared_ptr<sgpp::base::Grid> grid;
  std::unique_ptr<sgpp::base::SBasis> basis;
};

}  // namespace datadriven
}  // namespace sgpp

// #endif /* USE_EIGEN */
