// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/combigrid/pce/PolynomialChaosExpansion.hpp>
#include <sgpp/combigrid/pce/SGppToDakota.hpp>
#include <pecos_data_types.hpp>

#include <vector>

namespace sgpp {
namespace combigrid {

PolynomialChaosExpansion::PolynomialChaosExpansion(
    std::shared_ptr<CombigridOperation> combigridOperation,
    std::shared_ptr<AbstractInfiniteFunctionBasis1D> functionBasis)
    : numDims(combigridOperation->numDims()),
      functionBasis(functionBasis),
      combigridOperation(combigridOperation),
      combigridMultiOperation(nullptr),
      expansionCoefficientsFlag(false),
      sobolIndicesFlag(false) {
  // create tensor operation for pce transformation
  tensorOperation =
      sgpp::combigrid::CombigridTensorOperation::createOperationTensorPolynomialInterpolation(
          combigridOperation->getPointHierarchies(), combigridOperation->getStorage(),
          combigridOperation->getLevelManager(), functionBasis);
}

PolynomialChaosExpansion::PolynomialChaosExpansion(
    std::shared_ptr<CombigridMultiOperation> combigridMultiOperation,
    std::shared_ptr<AbstractInfiniteFunctionBasis1D> functionBasis)
    : numDims(combigridMultiOperation->numDims()),
      functionBasis(functionBasis),
      combigridOperation(nullptr),
      combigridMultiOperation(combigridMultiOperation),
      expansionCoefficientsFlag(false),
      sobolIndicesFlag(false) {
  // create tensor operation for pce transformation
  tensorOperation =
      sgpp::combigrid::CombigridTensorOperation::createOperationTensorPolynomialInterpolation(
          combigridMultiOperation->getPointHierarchies(), combigridMultiOperation->getStorage(),
          combigridMultiOperation->getLevelManager(), functionBasis);
}

PolynomialChaosExpansion::~PolynomialChaosExpansion() {}

double PolynomialChaosExpansion::mean() {
  //  if (orthogPoly != nullptr) {
  // #ifdef USE_DAKOTA
  //    return orthogPoly->mean();
  // #endif
  //  }
  if (!expansionCoefficientsFlag) {
    expansionCoefficients = tensorOperation->getResult();
    expansionCoefficientsFlag = true;
  }
  return expansionCoefficients.get(sgpp::combigrid::MultiIndex(numDims, 0)).getValue();
}

double PolynomialChaosExpansion::variance() {
  //  if (orthogPoly != nullptr) {
  // #ifdef USE_DAKOTA
  //    return orthogPoly->variance();
  // #endif
  //  }
  // PCE norm, i.e. variance (if using appropriate normalized orthogonal polynomials)
  if (!expansionCoefficientsFlag) {
    expansionCoefficients = tensorOperation->getResult();
    expansionCoefficientsFlag = true;
  }

  double var = 0.0;
  auto it = expansionCoefficients.getValues()->getStoredDataIterator();
  if (!it->isValid()) {
    return 0.0;
  }
  it->moveToNext();  // ignore first entry (belonging to mean)

  for (; it->isValid(); it->moveToNext()) {
    double coeff = it->value().value();
    var += coeff * coeff;
  }
  return var;
}

void PolynomialChaosExpansion::getComponentSobolIndices(
    sgpp::base::DataVector& componentSobolIndices, bool normalized) {
  if (!sobolIndicesFlag) {
    computeComponentSobolIndices();
    sobolIndicesFlag = true;
  }
  // copy sobol indices to output vector
  componentSobolIndices.resize(sobolIndices.getSize());
  componentSobolIndices.copyFrom(sobolIndices);

  // divide all the entries by the variance to obtain the Sobol indices
  if (normalized) {
    double var = variance();
    if (var > 1e-14) {
      componentSobolIndices.mult(1. / var);
    }
  }
}

void PolynomialChaosExpansion::computeComponentSobolIndices() {
  size_t numSobolIndices = static_cast<size_t>(std::pow(2, numDims) - 1);
  sobolIndices.resizeZero(numSobolIndices);

  // load index vectors
  auto it_coeffs = expansionCoefficients.getValues()->getStoredDataIterator();
  std::vector<MultiIndex> indexList;
  for (; it_coeffs->isValid(); it_coeffs->moveToNext()) {
    indexList.push_back(it_coeffs->getMultiIndex());
  }

  for (size_t i = 0; i < numSobolIndices; i++) {
    // loop over all the remaining basis functions
    for (std::vector<MultiIndex>::iterator it_ixlist = indexList.begin();
         it_ixlist != indexList.end();) {
      MultiIndex multiIndex = *it_ixlist;

      // check if all the current dimensions are set
      bool dimsAreSet = true;
      size_t idim = 0;
      while (dimsAreSet && idim < numDims) {
        // get kth bit from permuation mask
        size_t isSet = ((i + 1) & (1 << idim)) >> idim;  // in {0, 1}
        // check if the degree of the basis in the current dimension
        // is > 0 if set and equal to 0 if not set
        if (isSet == 0) {
          dimsAreSet &= multiIndex[idim] == 0;
        } else {
          dimsAreSet &= multiIndex[idim] > 0;
        }
        idim++;
      }

      // add the squared coefficient of the current term to the
      // sobol index matrix if all dimensions are present
      // in the current term
      if (dimsAreSet) {
        // load index and delete it from list since every term
        // can just contribute to one sobol index
        indexList.erase(it_ixlist);
        double coefficient = expansionCoefficients.get(multiIndex).value();
        sobolIndices[i] += coefficient * coefficient;
      } else {
        it_ixlist++;
      }
    }
  }
}

void PolynomialChaosExpansion::getTotalSobolIndices(sgpp::base::DataVector& totalSobolIndices,
                                                    bool normalized) {
  // compute the component sobol indices
  if (!sobolIndicesFlag) {
    computeComponentSobolIndices();
    sobolIndicesFlag = true;
  }

  totalSobolIndices.resizeZero(numDims);
  for (size_t idim = 0; idim < numDims; idim++) {
    for (size_t iperm = 0; iperm < sobolIndices.size(); iperm++) {
      // check if the current dimension is set in the current key of the sobol index
      size_t isSet = ((iperm + 1) & (1 << idim)) >> idim;  // in {0, 1}
      if (isSet == 1) {
        totalSobolIndices[idim] += sobolIndices[iperm];
      }
    }
  }

  // divide all the entries by the variance to obtain the Sobol indices
  if (normalized) {
    double var = variance();
    if (var > 1e-14) {
      totalSobolIndices.mult(1. / variance());
    }
  }
}

} /* namespace combigrid */
} /* namespace sgpp */
