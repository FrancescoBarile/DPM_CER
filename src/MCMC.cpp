#include <RcppArmadillo.h>
#include <cmath>
#include <unordered_map>
#include <vector>
#include "hamming.h"
#include "struct.h"
#include "tbeta.h"
#include "update_steps.h"
#include "wrapper.h"

using namespace Rcpp;

// [[Rcpp::export]]
Rcpp::List DPM_CER( const arma::imat  & G_mat,
                    const arma::ivec  & G_0,
                    const double a, const double b, const double c,
                    const int sample = 1000, const int burn = 0, const int thinning = 1){
  
  const int n = G_mat.n_rows; 
  const int M = G_mat.n_cols; 
  
  arma::vec G0 = arma::conv_to<arma::vec>::from(G_0);

  
  const int K_init = n;
  std::unordered_map<int, theta_cluster> theta_map;
  
  arma::ivec rho(n, arma::fill::zeros);
  
  for( int k = 0; k < K_init; ++k ){
    rho(k) = k ;
    
    double alpha_init = cpp_rtbeta(a, b);
    
    arma::vec p_k = arma::exp(
      (1 - G0) * std::log(alpha_init)
      + G0 * std::log(1.0 - alpha_init)
    );
    
    arma::ivec G_init = cpp_rmvbern(p_k);
    
    theta_map.emplace( k, theta_cluster{1, {G_init, alpha_init}} );
  }
  
  std::vector<arma::vec> phi(n);
  std::vector<arma::vec> tbeta_a_par(n);
  std::vector<arma::vec> tbeta_b_par(n);
  arma::vec pi_l0(n);
  for( int ell = 0; ell<n; ++ell){
    arma::ivec G_ell = G_mat.row(ell).t();
    int dH_ell = distH(G_ell, G_0);
    
    int M_d = M - dH_ell;
    
    phi[ell].set_size(M_d + 1);
    tbeta_a_par[ell].set_size(M_d + 1);
    tbeta_b_par[ell].set_size(M_d + 1);
    arma::vec sum_vec(M_d + 1);
    for( int r = 0; r <= M_d; ++r  ){
      double a_ell_r =  a + 2*r + dH_ell;
      double b_ell_r =  b + 2*M - 2*r - dH_ell;
      
      tbeta_a_par[ell](r) = a_ell_r;
      tbeta_b_par[ell](r) = b_ell_r;
      
      sum_vec(r) = (dH_ell* std::log(2.0) +  log_choose(M_d, r)) +
      cpp_log_incomplete_beta(a_ell_r, b_ell_r, 0.5) - cpp_log_incomplete_beta(a, b, 0.5);
      
      phi[ell](r) =
        std::exp( dH_ell * std::log(2.0)
        + log_choose(M_d, r)
        + cpp_log_incomplete_beta(a_ell_r, b_ell_r, 0.5)  );
    }
    pi_l0[ell] = std::exp(log_sum_exp(sum_vec) + std::log(c));
    phi[ell] /= arma::sum(phi[ell]);
  }
  
  Rcpp::List return_list(sample);
  unsigned ss = 0;
  const int S = burn + thinning * sample;
  
  std::vector<int> new_labels;
  
  for( int s = 0; s<S; ++s){
    
    for( int ell = 0; ell<n; ++ell){
      arma::ivec G_ell = G_mat.row(ell).t();
      
      update_allocation(
        ell,
        rho,
        theta_map,
        new_labels,
        G_ell,
        G_0,
        pi_l0,
        tbeta_a_par,
        tbeta_b_par,
        phi
      );
      
    }
    
    for (auto& kv : theta_map){
      int k = kv.first;
      
      theta_cluster& theta_k = kv.second;

      arma::ivec& G_k = theta_k.theta.first;
      double& alpha_k = theta_k.theta.second;

      alpha_k = update_alpha( k, G_k, rho, G_mat, G_0, a, b );
      
      G_k     = update_G( k, alpha_k, rho, G_mat, G_0 );
      
    }
    
    
    // Store
    if( ( (s+1) > burn ) & ((s+1-burn) % thinning == 0)){
     
      Rcpp::List theta_out = wrap_theta_map(theta_map);
      
      return_list[ss] = Rcpp::List::create(
        Rcpp::Named("rho")   = rho.t(),
        Rcpp::Named("label") = theta_out["label"],
        Rcpp::Named("size")  = theta_out["size"],
        Rcpp::Named("theta") = theta_out["theta"]
      );
      ss += 1 ;
    }
    
    
  }
  return return_list;
}