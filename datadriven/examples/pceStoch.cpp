#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/tools/DistributionBeta.hpp>
#include <sgpp/base/tools/DistributionTruncExponential.hpp>
#include <sgpp/base/tools/DistributionTruncGamma.hpp>
#include <sgpp/base/tools/DistributionTruncNormal.hpp>
#include <sgpp/datadriven/tools/PolynomialChaosExpansion.hpp>
#include <sgpp/globaldef.hpp>
#include <sgpp_base.hpp>
#include <sgpp_optimization.hpp>
#include <string>
#include <utility>
#include <vector>
// functions to be integrated
double e(const sgpp::base::DataVector& vec) {
  return vec[0] * vec[0] * vec[0] /*- vec[1] * vec[1]*/;
}
double f(const sgpp::base::DataVector& vec) { return 1.0; }
double g(const sgpp::base::DataVector& vec) {
  return 1 + (std::sin(vec[0]) + std::cos(vec[1])) / std::exp(vec[1]);
}
class Functf : public sgpp::base::ScalarFunction {
 public:
  explicit Functf() : sgpp::base::ScalarFunction(1) {}
  double eval(const sgpp::base::DataVector& x) { return 1.0; }
  virtual void clone(std::unique_ptr<sgpp::base::ScalarFunction>& clone) const {
    clone = std::unique_ptr<sgpp::base::ScalarFunction>(new Functf(*this));
  }
};
class Functe : public sgpp::base::ScalarFunction {
 public:
  explicit Functe() : sgpp::base::ScalarFunction(1) {}
  double eval(const sgpp::base::DataVector& vec) {
    return vec[0] * vec[0] * vec[0] /*- vec[1] * vec[1]*/;
  }
  virtual void clone(std::unique_ptr<sgpp::base::ScalarFunction>& clone) const {
    clone = std::unique_ptr<sgpp::base::ScalarFunction>(new Functe(*this));
  }
};

int main() {
  std::fstream of;
  of.open("plot_pce/stochasticTest.txt", std::ios::out | std::ios::trunc);
  of << std::fixed;
  of << std::setprecision(9);
  sgpp::datadriven::PolynomialChaosExpansion ee = sgpp::datadriven::PolynomialChaosExpansion(
      f, 1,
      std::vector<sgpp::datadriven::distributionType>{sgpp::datadriven::distributionType::Normal},
      std::vector<std::pair<double, double>>{
          std::pair<double, double>{-20, 20},
      },
      2, 2);
  std::cout << "-----------------------------------------------------------------------------------"
            << '\n';
  std::cout << "-----------------------------------------------------------------------------------"
            << '\n';

  std::cout << ee.getMean(200, "adaptiveGrid") << " pce mean" << '\n';
  std::cout << ee.getVariance(200, "adaptiveGrid") << " pce variance" << '\n';

  auto stuff = ee.getCoefficients();
  for (auto entry : stuff) {
    std::cout << entry << ", ";
  }
  std::cout << '\n';

  auto f1 = std::make_shared<Functf>();
  auto e1 = std::make_shared<Functe>();
  sgpp::optimization::SplineResponseSurface surface(f1, sgpp::base::DataVector{-20},
                                                    sgpp::base::DataVector{20},
                                                    sgpp::base::GridType::NakBsplineBoundary);
  surface.surplusAdaptive(200, 1);
  sgpp::base::DistributionsVector dists;
  auto dist1 = std::make_shared<sgpp::base::DistributionNormal>(0, 1);
  auto dist2 = std::make_shared<sgpp::base::DistributionNormal>(0, 1);
  dists.push_back(dist1);
  // dists.push_back(dist2);
  std::cout << surface.eval(sgpp::base::DataVector{0}) << ' ';
  std::cout << "surface eval" << '\n';
  std::cout << surface.getIntegral() << ' ';
  std::cout << "surface integral" << '\n';
  std::cout << surface.getMean(dists, 100) << ' ';
  std::cout << "surface mean" << '\n';
  std::cout << surface.getVariance(dists, 100) << ' ';
  std::cout << "surface variance" << '\n';

  std::cout << "-----------------------------------------------------------------------------------"
            << '\n';

  sgpp::datadriven::PolynomialChaosExpansion ee1 = sgpp::datadriven::PolynomialChaosExpansion(
      e, 3,
      std::vector<sgpp::datadriven::distributionType>{sgpp::datadriven::distributionType::Normal},
      std::vector<std::pair<double, double>>{
          std::pair<double, double>{-20, 20},
      },
      2, 2);
  std::cout << ee1.getMean(200, "adaptiveGrid") << " pce mean" << '\n';
  std::cout << ee1.getVariance(200, "adaptiveGrid") << " pce variance" << '\n';
  auto stuff2 = ee1.getCoefficients();
  for (auto entry : stuff2) {
    std::cout << entry << ", ";
  }
  std::cout << '\n';
  sgpp::optimization::SplineResponseSurface surface2(e1, sgpp::base::DataVector{-20},
                                                     sgpp::base::DataVector{20},
                                                     sgpp::base::GridType::NakBsplineBoundary);
  sgpp::optimization::SplineResponseSurface surface3(e1, sgpp::base::DataVector{-20},
                                                     sgpp::base::DataVector{20},
                                                     sgpp::base::GridType::NakBsplineBoundary);
  surface2.surplusAdaptive(200, 1);
  std::cout << surface2.eval(sgpp::base::DataVector{0}) << ' ';
  std::cout << "surface eval" << '\n';
  std::cout << surface2.getIntegral() << ' ';
  std::cout << "surface integral" << '\n';
  std::cout << surface2.getMean(dists, 100) << ' ';
  std::cout << "surface mean" << '\n';
  std::cout << surface2.getVariance(dists, 100) << ' ';
  std::cout << "surface variance" << '\n';
  std::cout << "-----------------------------------------------------------------------------------"
            << '\n';
  std::cout << "-----------------------------------------------------------------------------------"
            << '\n';
  of.close();
}
