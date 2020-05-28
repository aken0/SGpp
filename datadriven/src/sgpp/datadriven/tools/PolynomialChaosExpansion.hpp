#pragma once
#include <sgpp/globaldef.hpp>
#include <map>
#include <functional>
#include <string>
#include <cmath>
#include <utility>
#include <vector>
#include "sgpp/base/datatypes/DataVector.hpp"
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
      //double evalHermite(int n, double x);
      //double evalLegendre(double n, double x);
      //double evalLaguerre(int n, double x);
      //double evalJacobi(int n, double x, double alpha, double beta);
      //double evalGenLaguerre(int n, double x, double alpha);

      //double weightHermite(double n, double x);
      //double weightLegendre(double n, double x);
      //double weightLaguerre(double n, double x);
      //double weightJacobi(double n, double x, double alpha, double beta);
      //double weightGenlaguerre(double n, double x, double alpha);

      //      const std::map<std::string,std::function<double(double)>> weights {
      //       {"hermite",[](double x){return std::exp(-std::pow(x,2)/2);}},
      //        {"jacobi",[this](double x)->double{return std::pow((1-x),alpha)*std::pow((1+x),beta);}},
      //       {"legendre",[](double x){return 1;}},
      //      {"laguerre",[](double x){return std::exp(-x);}},
      //     {"genlaguerre",[this](double x) ->double {return std::pow(x,alpha)*std::exp(-x);}}
      //};
      //FIX VALUES!!!
      //      const std::map<std::string, std::function<double(double)>> denoms{
      //        {"hermite",[](double j){return std::sqrt(2*M_PI);}},
      //         {"jacobi",[this](double j){return 0;}},
      //          {"legendre",[](double j){return 2/((2*j)+1);}},
      //          {"laguerre",[](double j){return 1.0;}},
      //          {"genlaguerre",[this](double j){return std::pow(j,alpha)*std::exp(-j);}}
      //      };
      //      const std::map<std::string, std::function<double(double,double)>> evals {
      //        {"hermite",[this](int n, double x){return evalHermite(n, x);}},
      //          {"jacobi",[this](int n, double x){return evalJacobi(n, x,alpha,beta);}},
      //          {"legendre",[this](int n, double x){return evalLegendre(n, x);}},
      //          {"laguerre",[this](int n, double x){return evalLaguerre(n, x);}},
      //          {"genlaguerre",[this](int n, double x){return evalGenLaguerre(n, x,alpha);}}};
      protected:

      std::map<std::string,std::function<double(double)>> returnmap();
      std::map<std::string,std::function<double(double)>> weights;  
      std::map<std::string,std::function<double(double)>> denoms;  
      std::map<std::string,std::function<double(double,double)>> evals;  

      public:
      PolynomialChaosExpansion(std::function<double(const base::DataVector&)> func,int order, std::vector<std::string> types,std::vector<std::pair<double, double>> ranges,double alpha,double beta);
      ~PolynomialChaosExpansion();
      //move to private
      double evalLegendre(int n, double x);
      double evalHermite(int n, double x);
      double evalLaguerre(int n, double x);
      double evalJacobi(int n, double x, double alpha, double beta);
      double evalGenLaguerre(int n, double x, double alpha);
      std::vector<std::vector<int>> multiIndex(int dimension, int order);
      double monteCarloQuad(std::function<double(const base::DataVector&)> func, std::vector<std::pair<double, double>>, long n);
      void calculateCoefficients();
    };
  }//namespace datadriven
}//namespace sgpp
