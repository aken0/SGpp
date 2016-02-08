// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef WAVELETGRID_HPP
#define WAVELETGRID_HPP

#include <sgpp/base/grid/Grid.hpp>

#include <iostream>

#include <sgpp/globaldef.hpp>


namespace SGPP {
namespace base {

/**
 * grid with wavelet base functions
 */
class WaveletGrid : public Grid {
 protected:
  WaveletGrid(std::istream& istr);

 public:
  /**
   * Constructor of grid with wavelet base functions
   *
   * @param dim the dimension of the grid
   */
  WaveletGrid(size_t dim);

  /**
   * Destructor
   */
  virtual ~WaveletGrid() override;

  virtual SGPP::base::GridType getType() override;

  virtual const SBasis& getBasis() override;

  virtual GridGenerator* createGridGenerator() override;

  static Grid* unserialize(std::istream& istr);

};

}  // namespace base
}  // namespace SGPP

#endif /* WAVELETGRID_HPP */
