#include <RcppArmadillo.h>
#include "struct.h"
#include "wrapper.h"


Rcpp::List wrap_theta_map(
    const std::unordered_map<int, theta_cluster>& theta_map)
{
  int K = theta_map.size();
  
  // assume all G have the same dimension
  int M = theta_map.begin()->second.theta.first.n_elem;
  
  arma::mat theta(K, M + 1);
  arma::ivec size(K);
  arma::ivec label(K);
  
  int k = 0;
  
  for (const auto& kv : theta_map)
  {
    int lab = kv.first;
    
    const theta_cluster& cl = kv.second;
    
    const arma::ivec& G_k = cl.theta.first;
    double alpha_k = cl.theta.second;
    
    label(k) = lab;
    size(k)  = cl.size;
    
    // store G_k
    theta.row(k).cols(0, M-1) =
      arma::conv_to<arma::rowvec>::from(G_k);
    
    // store alpha_k
    theta(k, M) = alpha_k;
    
    k++;
  }
  
  return Rcpp::List::create(
    Rcpp::Named("label") = label.t(),
    Rcpp::Named("size")  = size.t(),
    Rcpp::Named("theta") = theta
  );
}
