/*****************************************************************************/
/* This file is part of sgpp, a program package making use of spatially      */
/* adaptive sparse grids to solve numerical problems                         */
/*                                                                           */
/* Copyright (C) 2008 Jörg Blank (blankj@in.tum.de)                          */
/* Copyright (C) 2009 Alexander Heinecke (Alexander.Heinecke@mytum.de)       */
/*                                                                           */
/* sgpp is free software; you can redistribute it and/or modify              */
/* it under the terms of the GNU Lesser General Public License as published  */
/* by the Free Software Foundation; either version 3 of the License, or      */
/* (at your option) any later version.                                       */
/*                                                                           */
/* sgpp is distributed in the hope that it will be useful,                   */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of            */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             */
/* GNU Lesser General Public License for more details.                       */
/*                                                                           */
/* You should have received a copy of the GNU Lesser General Public License  */
/* along with sgpp; if not, write to the Free Software                       */
/* Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA */
/* or see <http://www.gnu.org/licenses/>.                                    */
/*****************************************************************************/

#ifndef UNIDIRGRADIENT_HPP
#define UNIDIRGRADIENT_HPP

#include "grid/GridStorage.hpp"
#include "data/DataVector.h"

#include <vector>

#ifdef USEOMP
#include <omp.h>
#endif

namespace sg
{

/**
 * Unidirectional scheme with gradient used for building the laplace operator
 */
class UnidirGradient
{
public:
	/**
	 * Constructor
	 *
	 * @param storage a GridStorage object that contains the gird points
	 */
	UnidirGradient(GridStorage* storage) : storage(storage) {}

	/**
	 * Destructor
	 */
	virtual ~UnidirGradient() {}


	/**
	 * Starting point of the complete up-down scheme to calculate
	 * the laplace operator. The updown operation is started in
	 * every dimension of the grid.
	 *
	 * @param alpha contains the grid points coefficients
	 * @param result contains the result of the laplace operator
	 */
	virtual void updown(DataVector& alpha, DataVector& result)
	{
#ifdef USEOMP
#ifndef USEOMPTHREE
		result.setAll(0.0);

		#pragma omp parallel shared(result)
		{
			#pragma omp for schedule(static)
			for(size_t i = 0; i < storage->dim(); i++)
			{
				DataVector beta(result.getSize());

				this->updown(alpha, beta, storage->dim() - 1, i);

				#pragma omp critical
				{
					result.add(beta);
				}
			}
		}
#endif
#endif
#ifdef USEOMP
#ifdef USEOMPTHREE
		DataVector beta(result.getSize());
		result.setAll(0.0);

		for(size_t i = 0; i < storage->dim(); i++)
		{
			#pragma omp parallel
			{
				#pragma omp single nowait
				{
					this->updown_parallel(alpha, beta, storage->dim() - 1, i);
				}
			}
			result.add(beta);
		}
#endif
#endif
#ifndef USEOMP
		DataVector beta(result.getSize());
		result.setAll(0.0);

		for(size_t i = 0; i < storage->dim(); i++)
		{
			this->updown(alpha, beta, storage->dim() - 1, i);
			result.add(beta);
		}
#endif
	}

protected:
	typedef GridStorage::grid_iterator grid_iterator;

	/// Pointer to the grid's storage object
	GridStorage* storage;

#ifndef USEOMPTHREE
	/**
	 * Recursive procedure for updown(). In dimension <i>gradient_dim</i> the L2 scalar product of the
	 * gradients is used. In all other dimensions only the L2 scalar product.
	 *
	 * @param dim the current dimension
	 * @param gradient_dim the dimension in which to use the gradient
	 * @param alpha vector of coefficients
	 * @param result vector to store the results in
	 */
	virtual void updown(DataVector& alpha, DataVector& result, size_t dim, size_t gradient_dim)
	{
		if(dim == gradient_dim)
		{
			// this got its own function so we can use partial template specialization
			// if more than just down is needed
			gradient(alpha, result, dim, gradient_dim);
		}
		else
		{
			//Unidirectional scheme
			if(dim > 0)
			{
				// Reordering ups and downs
				// Use previously calculated ups for all future calculations
				// U* -> UU* and UD*

				DataVector temp(alpha.getSize());
				up(alpha, temp, dim);
				updown(temp, result, dim-1, gradient_dim);


				// Same from the other direction:
				// *D -> *UD and *DD

				DataVector result_temp(alpha.getSize());
				updown(alpha, temp, dim-1, gradient_dim);
				down(temp, result_temp, dim);


				//Overall memory use: 2*|alpha|*(d-1)

				result.add(result_temp);
			}
			else
			{
				// Terminates dimension recursion
				up(alpha, result, dim);

				DataVector temp(alpha.getSize());
				down(alpha, temp, dim);

				result.add(temp);
			}

		}
	}

	/**
	 * All calculations for gradient_dim. The gradient is recursivly applied to
	 * all dimension of the grid
	 *
	 * @todo add mathematical description
	 *
	 * @param alpha the coefficients of the grid points
	 * @param result the result of the operations
	 * @param dim the current dimension in the recursion
	 * @param gradient_dim the dimension in that the gradient is applied
	 */
	virtual void gradient(DataVector& alpha, DataVector& result, size_t dim, size_t gradient_dim)
	{
		//Unidirectional scheme
		if(dim > 0)
		{
			// Reordering ups and downs
			// Use previously calculated ups for all future calculations
			// U* -> UU* and UD*

			DataVector temp(alpha.getSize());
			upGradient(alpha, temp, dim);
			updown(temp, result, dim-1, gradient_dim);


			// Same from the other direction:
			// *D -> *UD and *DD

			DataVector result_temp(alpha.getSize());
			updown(alpha, temp, dim-1, gradient_dim);
			downGradient(temp, result_temp, dim);


			//Overall memory use: 2*|alpha|*(d-1)

			result.add(result_temp);
		}
		else
		{
			// Terminates dimension recursion
			upGradient(alpha, result, dim);

			DataVector temp(alpha.getSize());
			downGradient(alpha, temp, dim);

			result.add(temp);
		}

	}
#endif

#ifdef USEOMPTHREE
	/**
	 * Recursive procedure for updown(). In dimension <i>gradient_dim</i> the L2 scalar product of the
	 * gradients is used. In all other dimensions only the L2 scalar product.
	 *
	 * This version is parallelizes using the OpenMP 3 Task concept
	 *
	 * @param dim the current dimension
	 * @param gradient_dim the dimension in which to use the gradient
	 * @param alpha vector of coefficients
	 * @param result vector to store the results in
	 */
	virtual void updown_parallel(DataVector& alpha, DataVector& result, size_t dim, size_t gradient_dim)
	{
		if(dim == gradient_dim)
		{
			// this got its own function so we can use partial template specialization
			// if more than just down is needed
			gradient_parallel(alpha, result, dim, gradient_dim);
		}
		else
		{
			//Unidirectional scheme
			if(dim > 0)
			{
				// Reordering ups and downs
				// Use previously calculated ups for all future calculations
				// U* -> UU* and UD*

				DataVector temp(alpha.getSize());
				DataVector temp_two(alpha.getSize());
				DataVector result_temp(alpha.getSize());

				#pragma omp task
				{
					up(alpha, temp, dim);
					updown_parallel(temp, result, dim-1, gradient_dim);
				}

				// Same from the other direction:
				// *D -> *UD and *DD

				#pragma omp task
				{
					updown_parallel(alpha, temp_two, dim-1, gradient_dim);
					down(temp_two, result_temp, dim);
				}

				#pragma omp taskwait

				//Overall memory use: 3*|alpha|*(d-1)
				result.add(result_temp);
			}
			else
			{
				DataVector temp(alpha.getSize());

				// Terminates dimension recursion

				#pragma omp task
				up(alpha, result, dim);

				#pragma omp task
				down(alpha, temp, dim);

				#pragma omp taskwait

				result.add(temp);
			}

		}
	}

	/**
	 * All calculations for gradient_dim. The gradient is recursivly applied to
	 * all dimension of the grid
	 *
	 * This version is parallelizes using the OpenMP 3 Task concept
	 *
	 * @todo add mathematical description
	 *
	 * @param alpha the coefficients of the grid points
	 * @param result the result of the operations
	 * @param dim the current dimension in the recursion
	 * @param gradient_dim the dimension in that the gradient is applied
	 */
	virtual void gradient_parallel(DataVector& alpha, DataVector& result, size_t dim, size_t gradient_dim)
	{
		//Unidirectional scheme
		if(dim > 0)
		{
			// Reordering ups and downs
			// Use previously calculated ups for all future calculations
			// U* -> UU* and UD*

			DataVector temp(alpha.getSize());
			DataVector temp_two(alpha.getSize());
			DataVector result_temp(alpha.getSize());

			#pragma omp task
			{
				upGradient(alpha, temp, dim);
				updown_parallel(temp, result, dim-1, gradient_dim);
			}

			// Same from the other direction:
			// *D -> *UD and *DD

			#pragma omp task
			{
				updown_parallel(alpha, temp_two, dim-1, gradient_dim);
				downGradient(temp_two, result_temp, dim);
			}

			//Overall memory use: 3*|alpha|*(d-1)

			#pragma omp taskwait

			result.add(result_temp);
		}
		else
		{
			DataVector temp(alpha.getSize());
			// Terminates dimension recursion

			#pragma omp task
			upGradient(alpha, result, dim);

			#pragma omp task
			downGradient(alpha, temp, dim);

			#pragma omp taskwait

			result.add(temp);
		}
	}
#endif

	/**
	 * Up-step in dimension <i>dim</i> for \f$(\phi_i(x),\phi_j(x))_{L_2}\f$.
	 * Applies the up-part of the one-dimensional mass matrix in one dimension.
	 * Computes \f[\int_{x=0}^1  \phi_i(x) \sum_{j, l_i < l_j} \alpha_j \phi_j(x) dx.\f]
	 *
	 * @param dim dimension in which to apply the up-part
	 * @param alpha vector of coefficients
	 * @param result vector to store the results in
	 */
	virtual void up(DataVector& alpha, DataVector& result, size_t dim) = 0;

	/**
	 * Down-step in dimension <i>dim</i> for \f$(\phi_i(x),\phi_j(x))_{L_2}\f$.
	 * Applies the down-part of the one-dimensional mass matrix in one dimension.
	 * Computes \f[\int_{x=0}^1  \phi_i(x) \sum_{j, l_i\geq l_j} \alpha_j \phi_j(x) dx.\f]
	 *
	 * @param dim dimension in which to apply the down-part
	 * @param alpha vector of coefficients
	 * @param result vector to store the results in
	 */
	virtual void down(DataVector& alpha, DataVector& result, size_t dim) = 0;

	/**
	 * down-Gradient step in dimension <i>dim</i> applies the gradient of the mass matrix
	 * in one dimension
	 *
	 * @todo complete mathematical description
	 *
	 * @param alpha the coefficients of the gridpoints
	 * @param result vector with the result of this operation
	 * @param dim the dimension in that down-Gradient is applied
	 */
	virtual void downGradient(DataVector& alpha, DataVector& result, size_t dim) = 0;

	/**
	 * up-Gradient step in dimension <i>dim</i> applies the gradient of the mass matrix
	 * in one dimension
	 *
	 * @todo complete mathematical description
	 *
	 * @param alpha the coefficients of the gridpoints
	 * @param result vector with the result of this operation
	 * @param dim the dimension in that up-Gradient is applied
	 */
	virtual void upGradient(DataVector& alpha, DataVector& result, size_t dim) = 0;
};

}

#endif /* UNIDIRGRADIENT_HPP */
