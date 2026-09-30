#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
#include <boost/math/special_functions/beta.hpp>
#include <cmath>
#include "hamming.h"

using namespace Rcpp;

// [[Rcpp::export]]
int distH(const arma::ivec& G_1,
          const arma::ivec& G_2) {
  return arma::accu(arma::abs(G_1 - G_2));
}

double cpp_log_incomplete_beta(const double& a, const double& b, const double& d) {
  
  // log of regularized incomplete beta I_d(A,B)
  double log_I = std::log(boost::math::ibeta(a, b, d));
  
  // log of complete beta B(A,B)
  double log_B = std::lgamma(a) +
    std::lgamma(b) -
    std::lgamma(a + b);
  
  return log_I + log_B;
}

// [[Rcpp::export]]
double cpp_int_ker_basemeas_minus_max( const int N, const int d_bar, 
                                       const double a, const double b, const double d){
  const int M = N * (N - 1) / 2;
  const int bin_coeff = M - d_bar;
  
  arma::vec log_p(bin_coeff + 1);
  
  for(int r = 0; r <= bin_coeff; ++r){
    
    double log_binom =
      std::lgamma(bin_coeff + 1.0)
    - std::lgamma(r + 1.0)
    - std::lgamma(bin_coeff - r + 1.0);
    
    double log_beta =
    cpp_log_incomplete_beta(
      a + 2.0 * r + d_bar,
      b + 2.0 * M - 2.0 * r - d_bar,
      d
    );
    
    log_p(r) = log_binom + log_beta;
  }
  
  double max_p = log_p.max();
  
  double sum_p = 0.0;
  
  for(int r = 0; r <= bin_coeff; ++r){
    sum_p += std::exp(log_p(r) - max_p);
  }
  
  return d_bar * std::log(2.0)
    - cpp_log_incomplete_beta(a, b, d)
    + std::log(sum_p)
    + max_p;
}


double log_choose(const int& n, const int& r) {
  return std::lgamma(n + 1.0) -
    std::lgamma(r + 1.0) -
    std::lgamma(n - r + 1.0);
}

double log_sum_exp(const arma::vec& x){
  double m = x.max();
  return m + std::log(arma::sum(arma::exp(x - m)));
}

arma::ivec cpp_rmvbern(const arma::vec & p){
  arma::vec u = arma::randu<arma::vec>(p.n_elem);      // uniform [0,1)
  return arma::conv_to<arma::ivec>::from(u < p);       // elementwise comparison
}

int cpp_sample_1(const arma::ivec& vec, const arma::colvec& prob) {
  int ret = Rcpp::RcppArmadillo::sample(vec, 1, false, prob)(0);
  return ret;
}


// [[Rcpp::export]]
arma::imat sample_from_CER( const arma::ivec& G, const double alpha, const int n){
  arma::vec p =
    alpha + (1.0 - 2.0 * alpha) *
    arma::conv_to<arma::vec>::from(G);
  
  arma::imat G_group(n, G.n_elem);
  
  for (int i = 0; i < n; ++i) {
    G_group.row(i) = cpp_rmvbern(p).t();
  }
  
  return G_group;
}

// [[Rcpp::export]]
double log_lik_CER(const arma::ivec& G, const arma::ivec& G_m, const double alpha){
  int dH = distH(G, G_m);
  int M = G.n_elem;
  
  return dH * std::log(alpha)
    + (M - dH) * std::log(1.0 - alpha);
}
