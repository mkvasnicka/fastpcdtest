#include <Rcpp.h>
using namespace Rcpp;

// version for balanced panels

// [[Rcpp::export]]
double fastpcdtest_balanced_(NumericMatrix M) {
  NumericMatrix u0 = M; // de-meaned errors (u_{it} - \bar{u}_i)
  NumericVector sd(M.nrow()); // something as errors' std.
  // prepare u0 and sd
  for(int k = 0; k < M.nrow(); k++){
    u0(k, _) = u0(k, _) - mean(u0(k, _));
    sd(k) = sqrt(sum(pow(u0(k, _), 2)));
  }
  // calculate sum of correlations
  double CD; // CD statistic (the result)
  CD = 0;
  for(int j = 0; j < M.nrow() - 1; j++){
    for(int k = j + 1; k < M.nrow(); k++){
      CD = CD + sum(u0(j, _) * u0(k, _)) / (sd(j) * sd(k));
    }
  }
  // add correction
  CD = CD * sqrt(2.0 * M.ncol() / (M.nrow() * (M.nrow() - 1)));
  // return the result
  return CD;
}



// version for both balanced and unbalanced panels

//
struct corext {
  double correlation; // the correlation
  int common; // the number of common observations
};

corext quasi_correlation(NumericVector x, NumericVector y){
  static bool TMPprint = true;
  // the fields that are non-NA in both vectors
  LogicalVector common = (!is_na(x)) * (!is_na(y));
  // times  that are non-NA in both vectors
  int T = sum(common);
  // valid fields in x and y
  NumericVector validx = x[common];
  NumericVector validy = y[common];
  // averages of valies fields in x and y
  double mx = sum(validx) / T;
  double my = sum(validy) / T;
  // remove the averages from the valid fields
  validx = validx - mx;
  validy = validy - my;
  // compute the rho
  long double rho; // the result
  rho = sqrt(static_cast<long double>(T)) * sum(validx * validy) / sqrt(sum(pow(validx, 2)) * sum(pow(validy, 2)));
  // return the result
  corext result = {.correlation = rho, .common = T};
  return(result);
}

// [[Rcpp::export]]
NumericVector fastpcdtest_unbalanced_(NumericMatrix M, int min_common) {
  bool TMPprint = true; // TMP for debugging
  double CD = 0; // CD statistic (the result)
  corext quasi_cor; // two vector correlation times sqrt(T) + number of common observations
  int invalid_cases = 0; // the number of invalid cases
  // calculate sum of correlations
  for(int j = 0; j < M.nrow() - 1; j++){
    for(int k = j + 1; k < M.nrow(); k++){
      quasi_cor = quasi_correlation(M(j, _), M(k, _));
      if (quasi_cor.common < min_common)
        invalid_cases++;
      else
        CD = CD + quasi_cor.correlation;
      // debugging print
      if (TMPprint && is_true(any(is_na(NumericVector::create(CD))))) {
        Rcpp::Rcout << "First nan when j = " << j << " and k = " << k << std::endl;
        Rcpp::Rcout << "Quasi cor = " << quasi_cor.correlation << ", " << quasi_cor.common << std::endl;
        NumericVector x = M(j, _);
        NumericVector y = M(k, _);
        Rcpp::Rcout << "x = " << x << std::endl;
        Rcpp::Rcout << "y = " << y << std::endl;
        NumericVector xx = x - sum(x) / x.length();
        NumericVector yy = y - sum(y) / y.length();
        Rcpp::Rcout << "upr = " << sum(pow(xx, 2)) * sum(pow(yy, 2)) << std::endl;
        TMPprint = false;
      }
    }
  }
  // add correction
  Rcpp::Rcout << "Debug: CD pre manipulation " << CD << std::endl;
  CD = CD * sqrt(2.0 / (M.nrow() * (M.nrow() - 1) - 2 * invalid_cases));
  // return the result
  NumericVector result = NumericVector::create(CD, static_cast<double>(invalid_cases));
  return result;
}
