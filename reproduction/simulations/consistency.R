library(DPmCER)
library(igraph)

project_path <- paste0(sub(as.character(desc::desc_get("Package")), "", getwd()), "storing_folders")

if (!dir.exists( paste0(project_path, "/simulations") ) ) {
  dir.create( paste0(project_path, "/simulations") )
}

sim_path <- paste0(project_path, "/simulations/sim_consistency")

if (!dir.exists( sim_path )) {
  dir.create( sim_path )
}

input_path  <- paste0(sim_path, "/input")
output_path <- paste0(sim_path, "/output")


if (!dir.exists( input_path )) {
  dir.create(input_path)
}
if (!dir.exists( output_path )) {
  dir.create(output_path)
}

# number of nodes
N<-20
# number of edges
E<-choose(N, 2)

# number of clusters
k<-4
# cluster proportions in the first k-1 blocks
w<-rep(1/k, k-1)

# True parameters 
# Gm_j
set.seed(1)
# Scale-free
G_m_1<- sample_fitness_pl(N, E*0.2, exponent.out=2, exponent.in=-1)
set.seed(2)
# Small-world
G_m_2<-sample_smallworld(1, N, nei=5, p=0.2, 
                         loops = FALSE, multiple = FALSE)
# Stochastic block model
n_block <- 2; # number of blocks
p_betw_blocks<-0.1 # probability to assign an edge between two nodes belonging
# to two different blocks
p_block<- rep(p_betw_blocks, n_block)
if(n_block==2){p_block<-p_betw_blocks}
pm_temp <- from_low_adj_to_full_adj(p_block)
pm <- 1-rowSums(pm_temp)
pm <- `diag<-`(pm_temp, 1-rowSums(pm_temp))
block_sizes<-c(rep(floor(N/n_block)+1, N%%n_block), 
               floor(N*rep(1/n_block, n_block-N%%n_block)) )
set.seed(3)
G_m_3 <- sample_sbm(N, pref.matrix = pm, 
                    block.sizes = block_sizes  )

# ER with prob 0.3 (just sample from a CER with null graph mode so that 
# each G_ij^m_4 ~ Ber(0.3) )
set.seed(4)
G_m_4 <- graph_from_adjacency_matrix(from_low_adj_to_full_adj(
  sample_from_CER( rep(0, E), 0.3, 1 ) ), mode="undirected")


G_m_low_mat<-matrix(NA, k, E)
for(j in 1:k){
  G_m_low_mat[j,]<-from_full_adj_to_low_adj(as.matrix(as_adjacency_matrix(get(paste0("G_m_", j) ))))
}


# True alphas (alpha_j's)
al_vec <- c(0.25, 0.3, 0.35, 0.4)
alphas <- c(al_vec[1], al_vec[3],al_vec[2], al_vec[4])

for(j in 1:length(alphas)){
  assign(paste0("alpha_", j), alphas[j])
}

# Number of repetitions per scenario
N_rep <- 100

# Gibbs sampler iterations
sample <- 1000
burn   <- 200
thinning <- 1
# iterations <- burn+thinning*sample

# PRIOR SETTING
# Prior TBeta
a=b=1
# DP concentration
c=1

# sample size scenarios
n_vec<-c(40, 80, 120, 200)

n_cores = max( (1:(parallel::detectCores()))[ N_rep %% (1:(parallel::detectCores()))==0 ] )

save.image(paste0(input_path, "/common_env.RData"))

library(doParallel)
library(foreach)

myCluster <- makeCluster(n_cores, type = "PSOCK") 
registerDoParallel(myCluster)

st_all<-Sys.time()
for(s in 1:length(n_vec)){
  
  out_s_path <- paste0(output_path, "/scenario_", s)
  
  if (!dir.exists( out_s_path )) {
    dir.create( out_s_path )
  }
  
  n<-n_vec[s]
  n_j<-floor(n*c(w, 1-sum(w))) + c(rep(1, n-sum( floor(n*c(w, 1-sum(w))) )), 
                                   rep(0, k-( n-sum( floor(n*c(w, 1-sum(w))) ) ) ) )

  foreach(r = 1:N_rep,  .packages = c("DPMCER", "igraph") ) %dopar% {

    load(paste0(input_path, "/common_env.RData"))
    
    set.seed(r)
    data_arr<-array(NA, dim = c(N, N, n))
    idx_s<-c(1, (cumsum(n_j)+1)[1:(k-1)])
    idx_e<-cumsum(n_j)
    for(j in 1:k){
      for(l in idx_s[j]:idx_e[j]){
        data_arr[,,l]<-from_low_adj_to_full_adj( 
          sample_from_CER( from_full_adj_to_low_adj( 
            as.matrix(as_adjacency_matrix( get(paste0("G_m_", j) ) ) ) ), 
            get(paste0("alpha_", j)), 1 ) )
      }
    }
    
    G_mat <- t(sapply(1:dim(data_arr)[3], function(l){
      mat<-data_arr[,,l]
      out<-t(mat)[lower.tri(t(mat))]
      return(out)
    } ))
    
    # Prior CER
    G_0<-as.numeric(colMeans(G_mat)>=0.5)
    
    set.seed(r)
    dpm_cer_inference = DPM_CER(G_mat, G_0, a, b, c, sample, burn, thinning)
    
    out <- list(input=list(data=G_mat), output=dpm_cer_inference)
    
    saveRDS(out, paste0(out_s_path, "/output_rep_", r, ".rds"))
    rm( out, dpm_cer_inference, G_0, G_mat, data_arr )
    gc()
  }
}
en_all<-Sys.time()
en_all-st_all
stopCluster(myCluster)


library(doParallel)
library(foreach)

load(paste0(input_path, "/common_env.RData"))

d=0.5

# Sample the graph space G_N according to the true density of the obs
n_MC<-1000
set.seed(n_MC)
cl_lab<-sample(x=k, size=n_MC, replace = T, prob =  c(w, 1-sum(w))  )
prop<-table(cl_lab)

idx_s<-c(1, (cumsum(prop)+1)[1:(k-1)])
idx_e<-cumsum(prop)

G_space<-matrix(NA, n_MC, E)
for(j in 1:k){
  for(l in idx_s[j]:idx_e[j]){
    # set.seed(l)
    G_space[l,]<-sample_from_CER( G_m_low_mat[j,], alphas[j], 1 ) 
  }
}
sum(duplicated(G_space))


outnames <- list.files(output_path)

dist_arr<-array(NA, dim=c(N_rep, 2, length(outnames)))

myCluster <- makeCluster(n_cores, type = "PSOCK") 
registerDoParallel(myCluster)

st_all<-Sys.time()
for(s in 1:length(outnames)){
  
  n<-n_vec[s]
  n_j<-floor(n*c(w, 1-sum(w))) + c(rep(1, n-sum( floor(n*c(w, 1-sum(w))) )), 
                                   rep(0, k-( n-sum( floor(n*c(w, 1-sum(w))) ) ) ) )
  
  out_s_path <- paste0(output_path, "/scenario_", s)
  
  distances_s <- foreach(r = 1:N_rep, 
                         .combine = 'c', 
                         .packages = c("DPMCER") ) %dopar% {
    
    distances_distribution <- function(G_space, dpm_cer_inference, G_m_low_mat, alphas, w, n_j, N, G_0, a, b, c, d) {
      my_choose<-function(n, k){ exp( lgamma(n+1) - lgamma(k+1) - lgamma(n-k+1) ) }
      
      k<-length(alphas)
      E<-length(G_0)
      w_full<-c(w, 1-sum(w))
      n<-sum(n_j)
      
      true_dens<-sapply(1:nrow(G_space), function(g){
        sum(w_full*sapply(1:k, function(j){ 
          exp(log_lik_CER(G_space[g,], G_m_low_mat[j,], alphas[j])) } ))
      })
      
      int_vec<-sapply(1:nrow(G_space), function(g){
        d_bar<-distH(G_space[g,], G_0)
        integ<-exp(cpp_int_ker_basemeas_minus_max( N, d_bar, a, b, d))
      })
      
      data_ev_temp = sapply( 1:length(dpm_cer_inference), function(t){
        
        theta_star <- dpm_cer_inference[[t]][[4]]
        
        ev = sapply(1:nrow(theta_star), function(j){
          Gm_j_t_hat    <- as.numeric(theta_star[j,1:E])
          alpha_j_t_hat <- as.numeric(theta_star[j,E+1])
          
          n_j_t_hat <- as.numeric( dpm_cer_inference[[t]][[3]][j] )
          
          out=sapply(1:nrow(G_space), function(g){
            n_j_t_hat*exp(log_lik_CER(G_space[g,], Gm_j_t_hat, alpha_j_t_hat ))
          })
          return(out)
        })
        
        evid_data<-rowSums(ev )
        print(t)
        return(evid_data)
      } )
      
      data_ev<-rowMeans(data_ev_temp)
      
      est_den<-(int_vec + data_ev)/(c+n)
      
      L1_d <- mean( abs(est_den-true_dens)/true_dens )
      KL_d <- mean(  log(true_dens/est_den)  )
      return(c(L1_distance=L1_d, KL_divergence=KL_d))
    }
    
    load(paste0(input_path, "/common_env.RData"))
    
    out = readRDS( paste0(out_s_path, "/output_rep_", r, ".rds"))
    data_s_r = out[[1]][[1]]
    G_0_s_r<- as.numeric(colMeans( data_s_r ) >= 0.5)
    
    dpm_cer_inference = out$output
    
    dist_r<- distances_distribution(G_space, dpm_cer_inference, G_m_low_mat, alphas, w, n_j, N, G_0_s_r, a, b, c, d)
    return(dist_r)
  }
  dist_arr[,1,s]<-distances_s[which(names(distances_s)=="L1_distance")]
  dist_arr[,2,s]<-distances_s[which(names(distances_s)=="KL_divergence")]
  print(paste0("Scenario ", s))
}
en_all<-Sys.time()
en_all-st_all
stopCluster(myCluster)

save.image(paste0(input_path, "/L1_KL_dist_env.RData"))


library(tidyverse)
# devtools::install_github("teunbrand/ggh4x")
library(ggh4x) 
library(reporter)

load(paste0(input_path, "/L1_KL_dist_env.RData"))

distance_measures<-data.frame()
for(s in 1:dim(dist_arr)[3]){
  distance_measures<-rbind(distance_measures,  
                           cbind(rbind(cbind(data.frame(value=dist_arr[,1,s]), Distance="L1") , 
                                       cbind(data.frame(value=dist_arr[,2,s], Distance="KL") ) ),
                                 Scenario = paste0("Scenario ", s) )
  )
}
distance_measures$Model = "DPM"


pp<-ggplot(distance_measures, xlab=",",
           aes(x=Model, y=value, fill=Distance)) +  theme_bw() +
  theme(axis.text.x=element_blank(), axis.title.x=element_blank(), 
        axis.title.y=element_blank(), legend.position="none",
        legend.title=element_blank() ) +
  geom_violin() + 
  facet_grid2(Distance~Scenario, scales="free_y",  
              labeller = labeller(Distance = c("KL" = "KL divergence", 
                                               "L1" = paste("L" %p% supsc("1") , "distance")),
                                  Scenario = c("Scenario 1"= paste0("n = ", n_vec[1]),
                                               "Scenario 2"= paste0("n = ", n_vec[2]),
                                               "Scenario 3"= paste0("n = ", n_vec[3]),
                                               "Scenario 4"= paste0("n = ", n_vec[4])
                                  ) ) )   + 
  theme( strip.text.x = element_text(size = 15 ),
         strip.text.y.right  = element_text( size = 15 ),
         axis.text=element_text(size = 15) )

pp


ggsave(paste0(input_path, "/L1_KL_dist.pdf"), 
       plot = pp, width = 30, height = 15, units = "cm", dpi = 300)
