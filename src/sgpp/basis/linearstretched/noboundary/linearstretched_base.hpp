/******************************************************************************
 * Copyright (C) 2011 Technische Universitaet Muenchen                         *
 * This file is part of the SG++ project. For conditions of distribution and   *
 * use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
 ******************************************************************************/
// @author Jörg Blank (blankj@in.tum.de), Alexander Heinecke (Alexander.Heinecke@mytum.de), Sarpkan Selcuk (Sarpkan.Selcuk@mytum.de)


#ifndef LINEARSTRETCHED_BASE_HPP
#define LINEARSTRETCHED_BASE_HPP

#include <cmath>
#include <algorithm>
#include "grid/common/Stretching.hpp"

namespace sg
{

/**
 * linearstretched base functions.
 * And here we have another implicit dependence on tensor products
 */
template<class LT, class IT>
class linearstretched_base
{
public:
	/**
	 * Evaluate a base functions.
	 * Has a dependence on the absolute position of grid point and support.

	double eval(LT level, IT index, double p)
	{
		return 1.0 - fabs((1<<level) * p - index);
	}

	 *
	 * Evaluate a base functions.
	 * Has a dependence on the absolute position of grid point and support.
	 *
	 * This version catches errors, that occur if a basis function
	 * is evaluated outside its domain

	double evalSave(LT level, IT index, double p)
	{
		return std::max(1.0 - fabs((1<<level) * p - index), 0.0);
	}*/

	/*
	 * evaluate a base functions
	 * Has a dependence on the position of two grid points with values 1 and 0 and the
	 * support position
	 */
	double eval(double p, double pos0, double pos1){
		return (p-pos0)/(pos1-pos0);
	}
};

}

#endif /* LINEARSTRETCHED_BASE_HPP */
