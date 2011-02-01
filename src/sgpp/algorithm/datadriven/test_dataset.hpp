/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Jörg Blank (blankj@in.tum.de), Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef TEST_DATASET_HPP
#define TEST_DATASET_HPP

#include "grid/GridStorage.hpp"
#include "data/DataVector.hpp"

#include <vector>
#include <utility>
#include <iostream>

namespace sg {

/**
 * Returns the number of correctly classified instances in data without boundaries
 *
 * @param storage GridStorage object that contains the grid points
 * @param basis reference to class that implements to current basis
 * @param alpha the coefficients of the grid points
 * @param data the data the should be tested
 * @param classes the reference classes
 */
template<class BASIS>
double test_dataset( GridStorage* storage, BASIS& basis, DataVector& alpha, DataMatrix& data, DataVector& classes)
{
	typedef std::vector<std::pair<size_t, double> > IndexValVector;

	double correct = 0;

	#pragma omp parallel shared(correct)
	{
		size_t size = data.getNrows();

		std::vector<double> point;

		GetAffectedBasisFunctions<BASIS> ga(storage);

		#pragma omp for schedule(static)
		for(size_t i = 0; i < size; i++)
		{

			IndexValVector vec;
			double result = 0;

			data.getRow(i, point);

			ga(basis, point, vec);

			for(IndexValVector::iterator iter = vec.begin(); iter != vec.end(); iter++)
			{
				result += iter->second * alpha[iter->first];
			}

			if( (result >= 0 && classes[i] >= 0) || (result < 0 && classes[i] < 0) )
			{
				#pragma omp critical
				{
					correct++;
				}
			}
		}
	}

	return correct;
}

/**
 * Returns the MSE for given functions values at the evaluation points
 *
 * @param storage GridStorage object that contains the grid points
 * @param basis reference to class that implements to current basis
 * @param alpha the coefficients of the grid points
 * @param data the data the should be tested
 * @param refValues the function values at the evaluation points
 */
template<class BASIS>
double test_dataset_mse( GridStorage* storage, BASIS& basis, DataVector& alpha, DataMatrix& data, DataVector& refValues)
{
	typedef std::vector<std::pair<size_t, double> > IndexValVector;
	DataVector result(refValues.getSize());
	double mse = 0;

	#pragma omp parallel shared(result)
	{
		size_t size = data.getNrows();
		std::vector<double> point;
		GetAffectedBasisFunctions<BASIS> ga(storage);

		#pragma omp for schedule(static)
		for(size_t i = 0; i < size; i++)
		{

			IndexValVector vec;
			double res = 0;

			data.getRow(i, point);

			ga(basis, point, vec);

			for(IndexValVector::iterator iter = vec.begin(); iter != vec.end(); iter++)
			{
				res += iter->second * alpha[iter->first];
			}

			result[i] = res;
		}
	}

	result.sub(refValues);
	result.sqr();
	mse = result.sum();
	mse /= static_cast<double>(result.getSize());

	return mse;
}

/**
 * Returns the number of correctly classified instances in data without boundaries
 *
 * @param storage GridStorage object that contains the grid points
 * @param basis reference to class that implements to current basis
 * @param alpha the coefficients of the grid points
 * @param data the data the should be tested
 * @param classes the reference classes
 * @param charaNumbers the number of true positives, true negatives, false positives, false negatives (Vector of length 4)
 */
template<class BASIS>
double test_datasetWithCharacteristicNumber( GridStorage* storage, BASIS& basis, DataVector& alpha, DataMatrix& data, DataVector& classes, DataVector& charaNumbers)
{
	typedef std::vector<std::pair<size_t, double> > IndexValVector;

	double correct = 0;
	double tp = 0;
	double tn = 0;
	double fp = 0;
	double fn = 0;

	#pragma omp parallel shared(correct, tp, tn, fp, fn)
	{
		size_t size = data.getNrows();

		std::vector<double> point;

		GetAffectedBasisFunctions<BASIS> ga(storage);

		#pragma omp for schedule(static)
		for(size_t i = 0; i < size; i++)
		{

			IndexValVector vec;
			double result = 0;

			data.getRow(i, point);

			ga(basis, point, vec);

			for(IndexValVector::iterator iter = vec.begin(); iter != vec.end(); iter++)
			{
				result += iter->second * alpha[iter->first];
			}

			if( (result >= 0 && classes[i] >= 0) )
			{
				#pragma omp critical
				{
					tp++;
					correct++;
				}
			}
			else if( (result < 0 && classes[i] < 0) )
			{
				#pragma omp critical
				{
					tn++;
					correct++;
				}
			}
			else if( (result >= 0 && classes[i] < 0) )
			{
				#pragma omp critical
				{
					fp++;
				}
			}
			else // ( (result < 0 && classes[i] >= 0) )
			{
				#pragma omp critical
				{
					fn++;
				}
			}
		}
	}

	if (charaNumbers.getSize() < 4)
	{
		charaNumbers.resize(4);
	}

	charaNumbers.set(0, tp);
	charaNumbers.set(1, tn);
	charaNumbers.set(2, fp);
	charaNumbers.set(3, fn);

	return correct;
}

}

#endif /* TEST_DATASET_HPP */
