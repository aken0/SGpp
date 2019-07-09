/*
 * Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * DensityRatioEstimationMinerFactory.hpp
 *
 * Author: Paul Sarbu
 */

#pragma once

#include <sgpp/datadriven/datamining/builder/MinerFactory.hpp>
#include <sgpp/datadriven/datamining/configuration/DataMiningConfigParser.hpp>
#include <sgpp/datadriven/datamining/modules/dataSource/DataSource.hpp>
#include <sgpp/datadriven/datamining/modules/scoring/Scorer.hpp>
#include <sgpp/datadriven/datamining/modules/dataSource/DataSourceSplitting.hpp>
#include <sgpp/datadriven/datamining/modules/dataSource/DataSourceCrossValidation.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingBase.hpp>

#include <string>

namespace sgpp {
namespace datadriven {

/**
 * Concrete Factory that builds an instance of #sgpp::datadriven::SparseGridMiner for Least Squares
 * Regression
 */
class DensityRatioEstimationMinerFactory : public MinerFactory {
 public:
  /**
   * Default constructor
   */
  DensityRatioEstimationMinerFactory() = default;

  /**
   * Factory method to build a miner object based on a configuration file.
   * @param path Path to a configuration file that defines the structure of the miner object.
   */
  SparseGridMiner* buildMiner(const std::string& path) const override;

 private:
  /**
   * Factory method to build a splitting based data source, i.e. a data source that splits
   * data into validation and training data.
   * @param parser the datamining configuration parser instance to create the data source from
   * @return the data source instances (2 instances)
   */
  std::vector<DataSourceSplitting*> createDataSourceSplitting_TwoDatasets(
      const DataMiningConfigParser& parser) const override;

  /**
   * Factory method to build a cross validation data source, i.e. a data source that can separate
   * one fold from the data as validation set and use the rest for training
   * @param parser the datamining configuration parser instance to create the data source from
   * @return the data source instances (2 instances)
   */
  //  std::vector<DataSourceCrossValidation*> createDataSourceCrossValidation_TwoDatasets(
  //      const DataMiningConfigParser& parser) const;
  /**
   * Build an instance of a #sgpp::datadriven::ModelFittingBase_TwoDatasets object as specified in
   * the
   * configuration file.
   * @param parser parser object that provides methods to query the configuration file.
   * @return Fully configured fitter (instance of a #sgpp::datadriven::ModelFittingBase_TwoDatasets
   * object)
   * as specified in the configuration file.
   */
  ModelFittingBase* createFitter(const DataMiningConfigParser& parser) const override;

  FitterFactory* createFitterFactory(const DataMiningConfigParser& parser) const override {
    throw base::application_exception("HPO is not enabled for this model");
  }
};
} /* namespace datadriven */
} /* namespace sgpp */
