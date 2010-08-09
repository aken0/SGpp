/******************************************************************************
* Copyright (C) 2010 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef UPDPHIDPHIDOWNBBITERATIVELINEARBOUNDARY_HPP
#define UPDPHIDPHIDOWNBBITERATIVELINEARBOUNDARY_HPP

#include "grid/GridStorage.hpp"
#include "data/DataVector.hpp"

namespace sg
{

/**
 * This class is helper class to implement the complete Up
 * of following bilinearform \f$\int_{x} \frac{\partial \phi(x)}{x} \frac{\partial \phi(x)}{x} dx\f$
 * for a given dimension by an iterative algorithms on adaptive
 * Sparse Grids with linear ansatzfunctions with boundaries.
 *
 * This is possible due to the fact that the operator's
 * matrix has only entries on the diagonal.
 *
 * -> the Up/Down can be implemented by iterating over
 * all ansatzfunctions
 */
class UpdPhidPhiBBIterativeLinearBoundary
{
private:
	/// Pointer to the grid's storage object
	GridStorage* storage;

public:
	/**
	 * Constructor
	 *
	 * @param storage Pointer to the grid's storage object
	 */
	UpdPhidPhiBBIterativeLinearBoundary(GridStorage* storage);

	/**
	 * Destructor
	 */
	~UpdPhidPhiBBIterativeLinearBoundary();

	/**
	 * This operations performs the calculation of Up in the direction of dimension <i>dim</i>
	 * of following bilinearform: \f$\int_{x} \frac{\partial \phi(x)}{x} \frac{\partial \phi(x)}{x} dx\f$
	 *
	 * @param alpha DataVector that contains the gridpoint's coefficients
	 * @param result DataVector that contains the result of the down operation
	 * @param dim current fixed dimension of the 'execution direction'
	 */
	virtual void operator()(DataVector& alpha, DataVector& result, size_t dim);
};

}

#endif /* UPDPHIDPHIDOWNBBITERATIVELINEARBOUNDARY_HPP */
