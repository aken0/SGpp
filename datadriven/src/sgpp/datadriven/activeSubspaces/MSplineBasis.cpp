// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/activeSubspaces/MSplineBasis.hpp>

namespace sgpp {
namespace datadriven {

double MSplineBasis::eval(size_t degree, size_t index, double x) {
  //  const mSplineHashType hashKey = std::make_tuple(degree, index, x);
  //  std::map<mSplineHashType, double>::iterator it;
  //  it = precalculatedValues.find(hashKey);

  //  if (it != precalculatedValues.end()) {
  //    return it->second;
  //  } else {
  if (degree == 1) {
    if ((xi[index] <= x) && (x < xi[index + 1])) {
      double result = 1.0 / (xi[index + 1] - xi[index]);
      //        precalculatedValues.insert(std::pair<mSplineHashType, double>(hashKey, result));
      return result;
    } else {
      //        precalculatedValues.insert(std::pair<mSplineHashType, double>(hashKey, 0.0));
      return 0.0;
    }
  } else {
    if (xi[index + degree] != xi[index]) {
      double result = (static_cast<double>(degree) *
                       ((x - xi[index]) * eval(degree - 1, index, x) +
                        (xi[index + degree] - x) * eval(degree - 1, index + 1, x))) /
                      (static_cast<double>(degree - 1) * (xi[index + degree] - xi[index]));
      //      precalculatedValues.insert(std::pair<mSplineHashType, double>(hashKey, result));
      return result;
    } else {
      //        precalculatedValues.insert(std::pair<mSplineHashType, double>(hashKey, 0.0));
      return 0.0;
    }
  }
  //  }
}

double MSplineBasis::xpowplus(double x, size_t n) {
  if (x >= 0)
    return std::pow(x, static_cast<double>(n));
  else
    return 0.0;
}

double MSplineBasis::w(size_t v) {
  double res = 1.0;
  for (unsigned int i = 0; i < xi.size(); i++) {
    if (i != v) {
      // ToDo (rehmemk) diff can only be 0 if there are multiple knots. I should remove
      // duplicate knots in the projected corners
      double diff = xi[v] - xi[i];
      if (fabs(diff) > 1e-14) {
        res *= diff;
      }
    }
  }
  return res;
}

double MSplineBasis::evalTruncated(double x) {
  double res = 0.0;
  size_t n = xi.size() - 1;
  for (size_t v = 0; v < n + 1; v++) {
    res += (static_cast<double>(n) * xpowplus(xi[v] - x, n - 1)) / w(v);
  }
  return res;
}

}  // namespace datadriven
}  // namespace sgpp
