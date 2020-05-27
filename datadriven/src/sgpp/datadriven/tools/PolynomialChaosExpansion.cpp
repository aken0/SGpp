#include <cmath>
#include <functional>
#include <sgpp/datadriven/tools/PolynomialChaosExpansion.hpp>
#include <sgpp/globaldef.hpp>
#include "sgpp/base/datatypes/DataVector.hpp"
#include <random>
#include <string>
#include <vector>
#include <algorithm>
namespace sgpp{
  namespace datadriven{
    PolynomialChaosExpansion::PolynomialChaosExpansion(std::function<double(const base::DataVector&)> f,int order, std::vector<std::string> types, double alpha,double beta){
      func = f;
      this->alpha=alpha;
      this->beta=beta;
      this->order = order;
      this->types = types;
    }
    PolynomialChaosExpansion::~PolynomialChaosExpansion(){}

    std::map<std::string,std::function<double(double)>> PolynomialChaosExpansion::returnmap(){
      return PolynomialChaosExpansion::weights; 
    }

    double PolynomialChaosExpansion::evalHermite(int n, double x){
      if(n==0){
        return 1.0;
      }
      else if(n==1){
        return x;
      }
      else{
        double next=0.0;
        double last = 1.0;
        double curr = x;
        for (int i = 1.0; i<n; ++i) {
          next = x* curr - i* last; 
          last =curr;
          curr =next; 
        }
        return next;
      }
    }

    double PolynomialChaosExpansion::evalLegendre(int n, double x){
      if(n<0){
        n= -n -1;
      }
      if(n==0){
        return 1.0;
      }
      else if (n==1) {
        return x;
      }
      else {
        double next=0.0;
        //n-1
        double last=1.0;
        //n
        double curr=x;
        for (double i = 1.0; i<n; ++i) {
          next= ((2.0*i+1.0)/(i+1.0))*x*curr - (i/(i+1.0))*last;
          last = curr;
          curr = next;
        }
        return next;
      }

    }
    double PolynomialChaosExpansion::evalLaguerre(int n, double x){
      if(n==0){
        return 1.0;
      }
      else if (n==1) {
        return 1-x;
      }
      else {
        double next;
        //n-1
        double last=1.0;
        //n
        double curr=1-x;
        for (double i = 1.0; i<n; ++i) {
          next= ((2.0*i+1.0-x)*curr - i*last)/(i+1);
          last = curr;
          curr = next;
        }
        return next;
      }
      return 0;
    }
    //TODO FIX THIS
    double PolynomialChaosExpansion::evalJacobi(int n, double x, double alpha, double beta){
      if(n==0){
        return 1.0;
      }
      else if(n==1){
        return (alpha + 1)*(alpha+beta+2)*((x-1)/2);
      }
      else{
        double next=0.0;
        double last = 1.0;
        double curr=(alpha + 1)*(alpha+beta+2)*((x-1)/2);
        for (double i = 2.0; i<=n; ++i) {
          double q1= ((2*i*alpha+beta-1)*((2*i+alpha+beta)*(2*i+alpha+beta-2)*x+std::pow(alpha,2)-std::pow(beta,2)))/
            (2*i*(i+alpha+beta)*(2*i+alpha+beta-2));
          double q2= (2*(i+alpha-1)*(i+beta-1)*(2*i+alpha+beta))/
            (2*i*(i+alpha+beta)*(2*i+alpha+beta-2));
          next= q1*curr - q2*last;
          last = curr;
          curr = next;
        }
        return next;

      }
      return 0;
    }

    double PolynomialChaosExpansion::evalGenLaguerre(int n, double x, double alpha){
      if(n==0){
        return 1.0;
      }
      else if (n==1) {
        return 1+alpha-x;
      }
      else {
        double next=0.0;
        //n-1
        double last=1.0;
        //n
        double curr=1+alpha-x;
        for (double i = 1.0; i<n; ++i) {
          next= ((2.0*i+1.0+alpha-x)*curr - (i+alpha)*last)/(i+1);
          last = curr;
          curr = next;
        }
        return next;
      }
      return 0;
    }
    //TODO? reasonable accuracy
    double PolynomialChaosExpansion::monteCarloQuad(std::function<double(double)> func, double left, double right, long N){
      double factor=(right-left)*(1.0/(double(N)));
      std::random_device dev;
      std::mt19937 mersenne {dev()};
      std::uniform_real_distribution<double> dis {left,right};
      auto gen = [&func,&dis, &mersenne](){
        return func(dis(mersenne));};
      std::vector<double> vec(N);
      std::generate(vec.begin(),vec.end(),gen);
      double sum=std::accumulate(vec.begin(),vec.end(),0.0);
      return factor*sum;
    }
    std::vector<std::vector<int>> PolynomialChaosExpansion::multiIndex(int dimension, int order){
      std::vector<std::vector<int>> index;
      std::vector<int> curr(dimension);
      while (true){
        index.push_back(curr); 
        curr[0]++;
        for (int i = 0; i<dimension; ++i) {
          if (curr[i]>order){
            curr[i]-=order;
            curr[i+1]++;
          }
        }
        if(curr[dimension-1]>order){
          break;
        }
      }
      return index;
    }
    void PolynomialChaosExpansion::calculateCoefficients(){
      base::DataVector result;
      auto index =PolynomialChaosExpansion::multiIndex(this->types.size(),this->order);
      for (auto i= 0; i<index.size(); ++i) {

        //lambda for composite function to be integrated

        //integrate above function using MC
        double num;
        double denom;
        //calculate denominator

        double aj=num/denom; 
        result.push_back(aj);

      }
      


      this->coefficients=result;
    }


  }//namespace datadriven
}//namespace sgpp
