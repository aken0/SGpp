#pragma once
#include <sgpp/globaldef.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <map>
#include <functional>
#include <string>
#include <cmath>
#include <utility>
#include <vector>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
namespace sgpp{
  namespace datadriven{
    class PolynomialChaosExpansion{
      double alpha;
      double beta;
      int order;
      std::vector<std::string> types;
      std::vector<std::pair<double, double>> ranges;
      std::function<double(const base::DataVector&)> func;
      base::DataVector coefficients;

      private:
      double evalLegendre(int n, double x);
      double evalHermite(int n, double x);
      double evalLaguerre(int n, double x);
      double evalJacobi(int n, double x);
      double evalGenLaguerre(int n, double x);

      protected:
      std::map<std::string, std::function<double(double)>> weights;
      std::map<std::string, std::function<double(double)>> denoms;
      std::map<std::string, std::function<double(double, double)>> evals;

      public:
      PolynomialChaosExpansion(std::function<double(const base::DataVector&)> func, int order, std::vector<std::string> types, std::vector<std::pair<double, double>> ranges, double alpha, double beta);
      ~PolynomialChaosExpansion();
      // move to private

      std::vector<std::vector<int>> multiIndex(int dimension, int order);
      double monteCarloQuad(std::function<double(const base::DataVector&)> func,long n);
      void calculateCoefficients();
      base::DataVector getCoefficients();
      double evalExpansion();
    };
  }  // namespace datadriven
}  // namespace sgpp
