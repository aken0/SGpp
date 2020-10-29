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
// Stress:
double e(const sgpp::base::DataVector& vec) {
  return 1 - ((4 * vec[1]) / (5 * 225 * vec[0])) -
         ((vec[2] * vec[2]) / (25 * 225 * vec[0] * vec[0]));
}
class Functe : public sgpp::base::ScalarFunction {
 public:
  explicit Functe() : sgpp::base::ScalarFunction(3) {}
  double eval(const sgpp::base::DataVector& vec) { return e(vec); }
  virtual void clone(std::unique_ptr<sgpp::base::ScalarFunction>& clone) const {
    clone = std::unique_ptr<sgpp::base::ScalarFunction>(new Functe(*this));
  }
};

int main() {
  std::fstream of;
  of.open("plot_pce/stochasticTest.txt", std::ios::out | std::ios::trunc);
  of << std::fixed;
  of << std::setprecision(14);
  std::cout << std::fixed;
  std::cout << std::setprecision(14);
  sgpp::base::DistributionsVector dists;
  int dim = 3;
  // sgpp::base::DataVector lb{exp(.5), -1600, -400};
  // sgpp::base::DataVector ub{exp(9.5), 5600, 1400};
  sgpp::base::DataVector lb{exp(.5), -1600, -400};
  sgpp::base::DataVector ub{exp(9.5), 5600, 1400};
  auto dist1 = std::make_shared<sgpp::base::DistributionLogNormal>(5, .5);
  auto dist2 = std::make_shared<sgpp::base::DistributionNormal>(2000, 400);
  auto dist3 = std::make_shared<sgpp::base::DistributionNormal>(500, 100);
  dists.push_back(dist1);
  dists.push_back(dist2);
  dists.push_back(dist3);
  auto e1 = std::make_shared<Functe>();
  auto evalVec = dists.sample();
  std::cout << "-----------------------------------------------------------------------------------"
            << '\n';
  sgpp::datadriven::PolynomialChaosExpansion ee1 =
      sgpp::datadriven::PolynomialChaosExpansion(e, 3, dists);
  std::cout << ee1.getL2Error(2000, "adaptiveWeighted") << " pce L2" << '\n';
  std::cout << ee1.getMean(600, "adaptiveGrid") << " pce mean" << '\n';
  std::cout << ee1.getVariance(200, "adaptiveGrid") << " pce variance" << '\n';
  std::cout << ee1.evalExpansion(evalVec, 200, "adaptiveGrid") << " pce eval" << '\n';
  auto stuff2 = ee1.getCoefficients();
  for (auto entry : stuff2) {
    std::cout << entry << ", ";
  }
  std::cout << '\n';
  sgpp::optimization::SplineResponseSurface surface2(e1, lb, ub,
                                                     sgpp::base::GridType::NakBsplineExtended);
  surface2.surplusAdaptive(2000, 2);
  std::cout << surface2.eval(evalVec) << ' ';
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
