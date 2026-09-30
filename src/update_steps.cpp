#include <RcppArmadillo.h>
#include <cmath>
#include "hamming.h"
#include "struct.h"
#include "tbeta.h"
#include "update_steps.h"

void update_allocation(
    const int ell,
    arma::ivec& rho,
    std::unordered_map<int, theta_cluster>& theta_map,
    std::vector<int>& new_labels,
    const arma::ivec& G_ell,
    const arma::ivec& G0,
    const arma::vec& pi_l0,
    const std::vector<arma::vec>& tbeta_a_par,
    const std::vector<arma::vec>& tbeta_b_par,
    const std::vector<arma::vec>& phi
){
  
  const int rho_old = rho(ell);
  
  // Remove observation ell from its old cluster
  auto it = theta_map.find(rho_old);
  
  if (it->second.size == 1) {
    
    new_labels.push_back(rho_old);
    theta_map.erase(it);
    
  } else {
    
    it->second.size -= 1;
    
  }
  
  // const int rho_old = rho(ell);
  // 
  // // Remove observation ell from its old cluster
  // auto it = theta_map.find(rho_old);
  // 
  // theta_cluster& cl = it->second;
  // 
  // if (cl.size == 1) {
  //   
  //   new_labels.push_back(rho_old);
  //   theta_map.erase(it);
  //   
  // } else {
  //   
  //   cl.size -= 1;
  // }
  
  
  // Compute Polya urn probabilities
  arma::vec log_prob(theta_map.size() + 1);
  arma::ivec labels(theta_map.size() + 1);
  
  
  // New cluster
  log_prob(0) = std::log(pi_l0(ell));
  labels(0) = -1;
  
  
  // Existing clusters
  int j = 1;
  
  for(auto& [k, cl] : theta_map){
    
    int nk = cl.size;
    
    arma::ivec& G_k = cl.theta.first;
    double alpha_k = cl.theta.second;
    
    int dH = distH(G_ell, G_k);
    
    double log_like =
      dH * std::log(alpha_k)
      + (G_ell.n_elem - dH) *
        std::log(1.0 - alpha_k);
    
    
    log_prob(j) =
      std::log(nk) + log_like;
    
    labels(j) = k;
    
    j++;
  }
  
  
  // Normalize probabilities
  log_prob -= arma::max(log_prob);
  
  arma::vec prob = arma::exp(log_prob);
  
  prob /= arma::sum(prob);
  
  
  // Sample new allocation
  int sampled = cpp_sample_1(labels, prob);
  
  
  // New cluster
  if(sampled == -1){
    
    // int new_label = new_labels.back();
    // new_labels.pop_back();
    
    int new_label;
    
    if (!new_labels.empty())
    {
      new_label = new_labels.back();
      new_labels.pop_back();
    }
    else
    {
      new_label = rho.max() + 1;
    }  
      
    
    // Sample r
    arma::ivec r_values =
      arma::regspace<arma::ivec>(
        0,
        phi[ell].n_elem-1
      );
    
    int r =
      cpp_sample_1(r_values, phi[ell]);
    
    
    // Sample alpha | r,G
    double alpha_new =
      cpp_rtbeta(
        tbeta_a_par[ell](r),
        tbeta_b_par[ell](r)
      );
    
    
    // Sample consensus graph G_new
    arma::vec z =
      arma::conv_to<arma::vec>::from(
        G0 + G_ell - 1
      );
    
    double log_ratio =
      std::log(alpha_new)
      - std::log(1.0-alpha_new);
    
    
    arma::vec logit =
    2.0 * z * log_ratio;
    
    
    arma::vec p_new =
      1.0 / (1.0 + arma::exp(logit));
    
    
    arma::ivec G_new =
      cpp_rmvbern(p_new);
    
    
    // Add new cluster
    theta_map.emplace(
      new_label,
      theta_cluster{
        1,
        {G_new, alpha_new}
      }
    );
    
    
    rho(ell) = new_label;
    
  }
  // Existing cluster
  else{
    
    auto it_new = theta_map.find(sampled);
    
    it_new->second.size += 1;
    
    rho(ell) = sampled;
  }
}


double update_alpha(
    const int k,
    const arma::ivec& G_k,
    const arma::ivec& rho,
    const arma::imat& G_mat,
    const arma::ivec& G_0,
    const double a,
    const double b)
{
  
  int M = G_mat.n_cols;
  
  // indices of observations belonging to cluster k
  arma::uvec idx_k = arma::find(rho == k);
  
  // number of observations in cluster k
  int n_k = idx_k.n_elem;
  
  // extract all graphs belonging to cluster k
  arma::imat G_cluster = G_mat.rows(idx_k);
  
  // sum of Hamming distances:
  // sum_i d_H(G_i, G_k)
  int dist_sum = arma::accu(
    arma::abs(G_cluster.each_row() - G_k.t())
  );
  
  dist_sum += arma::accu(
    arma::abs(G_0 - G_k)
  );
  
  // Beta posterior parameters
  double a_star = a + dist_sum;
  
  double b_star = b + (n_k + 1) * M - dist_sum;
  
  // sample from truncated Beta(1/2, a_star, b_star)
  return cpp_rtbeta(a_star, b_star);
}


arma::ivec update_G(
    const int k,
    const double alpha_k,
    const arma::ivec& rho,
    const arma::imat& G_mat,
    const arma::ivec& G_0)
{
  // observations assigned to cluster k
  arma::uvec idx_k = arma::find(rho == k);
  
  int n_k = idx_k.n_elem;
  
  // extract cluster graphs
  arma::imat G_cluster = G_mat.rows(idx_k);
  
  // n_{k,m}: number of ones for each coordinate m
  arma::Row<int> n_km_int = arma::sum(G_cluster, 0) + G_0.t();
  arma::rowvec n_km = arma::conv_to<arma::rowvec>::from(n_km_int);
  
  // exponent:
  // 2 * (n_{k,m} - (n_k+1)/2)
  arma::rowvec exponent =
    2.0 * (n_km - (n_k + 1.0) / 2.0);
  
  // log(alpha/(1-alpha))
  double log_odds = std::log(alpha_k / (1.0 - alpha_k));
  
  // probabilities p_{k,m}*
  arma::rowvec p_row =
    1.0 / (1.0 + arma::exp(exponent * log_odds));
  
  // convert rowvec -> vec for cpp_rmvbern
  arma::vec p = p_row.t();
  
  // sample G_k
  return cpp_rmvbern(p);
}
