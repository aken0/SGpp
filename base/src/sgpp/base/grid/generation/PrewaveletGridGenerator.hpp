// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef PREWAVELETGRIDGENERATOR_HPP
#define PREWAVELETGRIDGENERATOR_HPP

#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/grid/generation/GridGenerator.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
namespace base {

/**
 * GridGenerator for prewavelet grids without boundaries
 */
class PrewaveletGridGenerator: public GridGenerator {
 protected:
  /// pointer to the storage object
  GridStorage* storage;
  GridStorage* shadowstorage;
  typedef GridStorage::index_type index_type;
  typedef index_type::index_type index_t;
  typedef index_type::level_type level_t;

 private:

  /**
   * For the shadow storage, the two left and two right neighbors in each
   * dimension of the refined point are required. This function only adds
   * the point which are not in the actual grid to the shadow storage.
   *
   * @param index point added during refinement
   * @param current_dim current dimension
   * @param target_level target level
   * @param iter iterator for grid storage
   * @param shadowIter iterator for shadow storage
   */
  void addNeighbours(index_type& index, size_t current_dim,
                     level_t target_level, GridStorage::grid_iterator& iter,
                     GridStorage::grid_iterator& shadowIter);

  /**
   * This function ensures that the special adaptive prewavelet grid points have parents.
   *
   * @param iter iterator for grid storage
   * @param shadowIter iterator for shadow storage
   */
  void insertParents(GridStorage::grid_iterator& iter,
                     GridStorage::grid_iterator& shadowIter);

  /**
   * If during the refinement one or more points of the shadow register are added
   * to the actual grid then we have to remove these points from the shadow storage.
   */
  void consolidateShadow();

 public:
  /**
   * Constructor
   * An adaptive grid with prewavelet ansatz functions requires for operations
   * using the up-down algorithm shadow points. These shadow points a needed just
   * for data transport, thus they do not have an influence on the final function.
   * Please refer to SGPP::pde::UpDownOneOpDimWithShadow for more information.
   *
   * @param storage the grid storage object of the the grid, on which the hierarchisation should be executed
   * @param shadowstorage shadow points (see detailed description)
   */
  PrewaveletGridGenerator(GridStorage* storage, GridStorage* shadowstorage);

  /**
   * Destructor
   */
  virtual ~PrewaveletGridGenerator() override;

  virtual void regular(size_t level) override;
  virtual void cliques(size_t level, size_t clique_size) override;
  virtual void full(size_t level) override;
  virtual void refine(RefinementFunctor* func) override;
  virtual size_t getNumberOfRefinablePoints() override;

  virtual void coarsen(CoarseningFunctor* func, DataVector* alpha) override;
  virtual void coarsenNFirstOnly(CoarseningFunctor* func, DataVector* alpha,
                                 size_t numFirstOnly) override;
  virtual size_t getNumberOfRemovablePoints() override;

  virtual void refineMaxLevel(RefinementFunctor* func, size_t maxLevel) override;
  virtual size_t getNumberOfRefinablePointsToMaxLevel(size_t maxLevel) override;

};

}  // namespace base
}  // namespace SGPP

#endif /* PREWAVELETGRIDGENERATOR_HPP */
