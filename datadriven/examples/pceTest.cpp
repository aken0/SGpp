#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/tools/DistributionTruncExponential.hpp>
#include <sgpp/base/tools/DistributionTruncGamma.hpp>
#include <sgpp/base/tools/DistributionTruncNormal.hpp>
#include <sgpp/datadriven/tools/PolynomialChaosExpansion.hpp>
#include <sgpp/globaldef.hpp>
#include <sgpp_base.hpp>
#include <string>
#include <utility>
#include <vector>
// functions to be integrated
double e(const sgpp::base::DataVector& vec) { return vec[0] * vec[0] * vec[0] - vec[1] * vec[1]; }
double f(const sgpp::base::DataVector& vec) { return 1.0; }
double g(const sgpp::base::DataVector& vec) {
  return 1 + (std::sin(vec[0]) + std::cos(vec[1])) / std::exp(vec[1]);
}
int main() {
  std::fstream of;
  std::string path;
  std::cout << "enter path: " << '\n';
  std::cin >> path;
  of.open("plot_pce/sampleDist.txt", std::ios::out | std::ios::trunc);
  of << std::fixed;
  of << std::setprecision(9);
  sgpp::base::DistributionTruncExponential distN(100);
  for (int i = 0; i < 1000000; ++i) {
    of << distN.sample() << ',';
  }
  /*
  sgpp::datadriven::PolynomialChaosExpansion ee = sgpp::datadriven::PolynomialChaosExpansion(
      g, 10,
      std::vector<sgpp::datadriven::distributionType>{sgpp::datadriven::distributionType::Uniform,
                                                      sgpp::datadriven::distributionType::Uniform},
      std::vector<std::pair<double, double>>{
          std::pair<double, double>{-1, 1},
          std::pair<double, double>{-1, 1},
      },
      0, 0);
  std::cout << ee.sparseGridQuadrature(e, 2, 60000) << '\n';
  */
  /*
  auto re = ee.calculateCoefficients();
  for (auto entry : re) {
    of << entry << ',';
  }
  */
  of << '\n';
  of.close();
}
