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
  return std::pow(vec[0], 3) + vec[1] * vec[1] + vec[2];
}
double f(const sgpp::base::DataVector& vec) { return 1.0; }
double g(const sgpp::base::DataVector& vec) {
  return 1 + (std::sin(vec[0]) + std::cos(vec[1])) / std::exp(vec[1]);
}
class Functf : public sgpp::base::ScalarFunction {
 public:
  explicit Functf(int dim) : sgpp::base::ScalarFunction(dim) {}
  double eval(const sgpp::base::DataVector& x) { return f(x); }
  virtual void clone(std::unique_ptr<sgpp::base::ScalarFunction>& clone) const {
    clone = std::unique_ptr<sgpp::base::ScalarFunction>(new Functf(*this));
  }
};
class Functe : public sgpp::base::ScalarFunction {
 public:
  explicit Functe(int dim) : sgpp::base::ScalarFunction(dim) {}
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
  double mean = 0;
  double sigma = 1;
  double l = (mean - 9 * sigma);
  double r = (mean + 9 * sigma);
  // double l = -1;
  // double r = 1;
  auto dist1 = std::make_shared<sgpp::base::DistributionNormal>(mean, sigma);
  auto dist2 = std::make_shared<sgpp::base::DistributionNormal>(mean, sigma);
  auto dist3 = std::make_shared<sgpp::base::DistributionNormal>(mean, sigma);
  // auto dist2 = std::make_shared<sgpp::base::DistributionUniform>(l, r);
  dists.push_back(dist1);
  dists.push_back(dist1);
  dists.push_back(dist3);
  sgpp::datadriven::PolynomialChaosExpansion ee =
      sgpp::datadriven::PolynomialChaosExpansion(f, 3, dists);
  std::cout << "-----------------------------------------------------------------------------------"
            << '\n';
  std::cout << "-----------------------------------------------------------------------------------"
            << '\n';

  std::cout << ee.getMean(300, "adaptiveGrid") << " pce mean" << '\n';
  std::cout << ee.getVariance(400, "adaptiveGrid") << " pce variance" << '\n';
  std::cout << ee.evalExpansion(sgpp::base::DataVector(dim, 0), 200, "adaptiveGrid") << " pce eval"
            << '\n';
  std::cout << ee.evalExpansion(sgpp::base::DataVector(dim, mean), 200, "adaptiveGrid")
            << " pce eval @mean" << '\n';
  auto stuff = ee.getCoefficients();
  for (auto entry : stuff) {
    std::cout << entry << ", ";
  }
  std::cout << '\n';

  auto f1 = std::make_shared<Functf>(dim);
  auto e1 = std::make_shared<Functe>(dim);
  sgpp::optimization::SplineResponseSurface surface(f1, sgpp::base::DataVector(dim, l),
                                                    sgpp::base::DataVector(dim, r),
                                                    sgpp::base::GridType::NakBsplineBoundary);
  surface.surplusAdaptive(200, 1);
  std::cout << surface.eval(sgpp::base::DataVector(dim, 0)) << ' ';
  std::cout << "surface eval" << '\n';
  std::cout << surface.eval(sgpp::base::DataVector(dim, mean)) << ' ';
  std::cout << "surface eval @ mean" << '\n';
  std::cout << surface.getIntegral() << ' ';
  std::cout << "surface integral" << '\n';
  std::cout << surface.getMean(dists, 100) << ' ';
  std::cout << "surface mean" << '\n';
  std::cout << surface.getVariance(dists, 100) << ' ';
  std::cout << "surface variance" << '\n';

  std::cout << "-----------------------------------------------------------------------------------"
            << '\n';

  sgpp::datadriven::PolynomialChaosExpansion ee1 =
      sgpp::datadriven::PolynomialChaosExpansion(e, 5, dists);
  std::cout << ee1.getMean(400, "adaptiveGrid") << " pce mean" << '\n';
  std::cout << ee1.getVariance(600, "adaptiveGrid") << " pce variance" << '\n';
  std::cout << ee1.evalExpansion(sgpp::base::DataVector(dim, 0), 600, "adaptiveGrid") << " pce eval"
            << '\n';
  std::cout << ee1.evalExpansion(sgpp::base::DataVector(dim, mean), 600, "adaptiveGrid")
            << " pce eval@mean" << '\n';
  auto stuff2 = ee1.getCoefficients();
  for (auto entry : stuff2) {
    std::cout << entry << ", ";
  }
  std::cout << '\n';
  sgpp::optimization::SplineResponseSurface surface2(e1, sgpp::base::DataVector(dim, l),
                                                     sgpp::base::DataVector(dim, r),
                                                     sgpp::base::GridType::NakBsplineBoundary);
  surface2.surplusAdaptive(200, 1);
  std::cout << surface2.eval(sgpp::base::DataVector(dim, 0)) << ' ';
  std::cout << "surface eval" << '\n';
  std::cout << surface2.eval(sgpp::base::DataVector(dim, mean)) << ' ';
  std::cout << "surface eval @ mean" << '\n';
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
