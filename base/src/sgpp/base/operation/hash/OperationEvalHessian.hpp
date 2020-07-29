// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef OPERATIONEVALHESSIAN_HPP
#define OPERATIONEVALHESSIAN_HPP

#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/datatypes/DataMatrix.hpp>

#include <vector>

namespace sgpp {
namespace base {

/**
 * Abstract operation for evaluating a linear combination of basis functions, its gradient
 * and its Hessian.
 */
class OperationEvalHessian {
 public:
  /**
   * Constructor.
   */
  OperationEvalHessian() {
  }

  /**
   * Destructor.
   */
  virtual ~OperationEvalHessian() {
  }

  /**
   * @param       alpha     coefficient vector
   * @param       point     evaluation point
   * @param[out]  gradient  gradient vector of the linear combination
   * @param[out]  hessian   Hessian matrix of the linear combination
   * @return                value of the linear combination
   */
  virtual double evalHessian(const DataVector& alpha,
                              const DataVector& point,
                              DataVector& gradient,
                              DataMatrix& hessian) = 0;

  /**
   * @param       alpha     coefficient matrix (each column is a coefficient vector)
   * @param       point     evaluation point
   * @param[out]  value     values of the linear combination
   * @param[out]  gradient  Jacobian of the linear combination (each row is a gradient vector)
   * @param[out]  hessian   vector of Hessians of the linear combination
   */
  virtual void evalHessian(const DataMatrix& alpha,
                           const DataVector& point,
                           DataVector& value,
                           DataMatrix& gradient,
                           std::vector<DataMatrix>& hessian) {
    const size_t d = point.getSize();
    const size_t m = alpha.getNcols();
    DataVector curAlpha(alpha.getNrows());
    DataVector curGradient(d);
    DataMatrix curHessian(d, d);

    value.resize(m);
    gradient.resize(d, d);

    if (hessian.size() != m) {
      hessian.resize(m);
    }

    for (size_t j = 0; j < m; j++) {
      DataMatrix& curHessian = hessian[j];
      curHessian.resize(d, d);
      alpha.getColumn(j, curAlpha);
      value[j] = evalHessian(curAlpha, point, curGradient, curHessian);
      gradient.setRow(j, curGradient);
    }
  }
  /// untransformed evaluation point (temporary vector)
  DataVector pointInUnitCube;
};

}  // namespace base
}  // namespace sgpp

#endif /* OPERATIONEVALHESSIAN_HPP */
