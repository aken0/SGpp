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
  auto a = 7;
  auto b = .1;
  return sin(vec[0]) + a * pow(sin(vec[1]), 2) + b * pow(vec[2], 4) * sin(vec[0]);
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
  of << std::setprecision(9);
  std::cout << std::fixed;
  std::cout << std::setprecision(9);
  sgpp::base::DistributionsVector dists;
  int dim = 3;
  sgpp::base::DataVector lb(dim, -M_PI);
  sgpp::base::DataVector ub(dim, M_PI);
  auto dist1 = std::make_shared<sgpp::base::DistributionUniform>(-M_PI, M_PI);
  auto dist2 = std::make_shared<sgpp::base::DistributionUniform>(-M_PI, M_PI);
  auto dist3 = std::make_shared<sgpp::base::DistributionUniform>(-M_PI, M_PI);
  dists.push_back(dist1);
  dists.push_back(dist2);
  dists.push_back(dist3);
  auto e1 = std::make_shared<Functe>();
  auto evalVec = dists.sample();
  std::cout << "-----------------------------------------------------------------------------------"
            << '\n';
  sgpp::datadriven::PolynomialChaosExpansion ee1 =
      sgpp::datadriven::PolynomialChaosExpansion(e, 1, dists);
  std::cout << ee1.getL2Error(1000, "adaptiveWeighted") << " pce L2" << '\n';
  std::cout << ee1.getMean(1600, "adaptiveGrid") << " pce mean" << '\n';
  std::cout << ee1.getVariance(200, "adaptiveGrid") << " pce variance" << '\n';
  std::cout << ee1.evalExpansion(evalVec, 200, "adaptiveGrid") << " pce eval" << '\n';
  auto stuff2 = ee1.getCoefficients();
  for (auto entry : stuff2) {
    std::cout << entry << ", ";
  }
  std::cout << '\n';
  sgpp::optimization::SplineResponseSurface surface2(e1, lb, ub,
                                                     sgpp::base::GridType::NakBsplineBoundary, 3);
  surface2.surplusAdaptive(1600, 1);
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
