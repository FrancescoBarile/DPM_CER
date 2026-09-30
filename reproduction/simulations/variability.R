library(DPmCER)
library(igraph)

project_path <- paste0(sub(as.character(desc::desc_get("Package")), "", getwd()), "storing_folders")

if (!dir.exists( paste0(project_path, "/simulations") ) ) {
  dir.create( paste0(project_path, "/simulations") )
}

sim_path <- paste0(project_path, "/simulations/sim_variability")

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

# sample size
n<-40
# number of clusters
k<-4
# cluster proportions in the first k-1 blocks
w<-rep(1/k, k-1)
# clusters size
n_j<-c(ceiling(n*w), n-sum(ceiling(n*w)) )

true_partition<-as.numeric(unlist(sapply(1:length(n_j), function(j){rep(j, n_j[j])})))

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

# ER with prob alpha_4 (just sample from a CER with null graph mode so that 
# each G_ij^m_4 ~ Ber(0.3) )
set.seed(4)
G_m_4 <- graph_from_adjacency_matrix(from_low_adj_to_full_adj(
  sample_from_CER( rep(0, E), 0.3, 1)), mode="undirected")

# Variability scenarios matrix (alpha_j's)
al_vec<-c(0.25, 0.3, 0.35, 0.4)
scenarios<-matrix(NA, nrow = (length(al_vec)+1), ncol= k)
for(i in 1:(nrow(scenarios)-1)){
  scenarios[i,]<-rep(al_vec[i], k)
}
scenarios[nrow(scenarios),] <- c(al_vec[1], al_vec[3],al_vec[2], al_vec[4])

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

n_cores = max( (1:(parallel::detectCores()))[ N_rep %% (1:(parallel::detectCores()))==0 ] )

save.image(paste0(input_path, "/common_env.RData"))

library(doParallel)
library(foreach)

myCluster <- makeCluster(n_cores, type = "PSOCK") 
registerDoParallel(myCluster)

st_all<-Sys.time()
for(s in 1:nrow(scenarios)){
  
  out_s_path <- paste0(output_path, "/scenario_", s)
  
  if (!dir.exists( out_s_path )) {
    dir.create( out_s_path )
  }
  
  foreach(r = 1:N_rep, .packages = c("DPMCER", "igraph")) %dopar% {

    load(paste0(input_path, "/common_env.RData"))
    
    alpha_s = scenarios[s,]
    
    for(j in 1:length(alpha_s)){
      assign(paste0("alpha_", j), alpha_s[j])
    }  
    
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


library(salso)
library(fossil) # rand.index
library(NMF) # purity and entropy of a clustering

load(paste0(input_path, "/common_env.RData"))

outnames <- list.files(output_path)

cls_measures<-array(NA, dim = c(N_rep, length(outnames), 3 )) # 3 clustering measures
for(s in 1:length(outnames)){
  
  out_s <- list.files( paste0(output_path, "/", outnames[s] ) )
  
  for(r in 1:length(out_s)){
    
    out_list <- readRDS( paste0(output_path, "/", outnames[s], "/", out_s[r]  )  )
    
    dpm_cer_inference = out_list$output
    
    
    PARTITION <- t(sapply(
      1:length(dpm_cer_inference),
      function(t)
        dpm_cer_inference[[t]]$rho
    ))
    
    optimal_clustering <- as.numeric(salso(PARTITION, loss=VI(a=1), 
                                           nRuns=10, nCores=0) )
    
    cls_measures[r,s,1]<-adj.rand.index(true_partition, optimal_clustering)
    cls_measures[r,s,2]<-purity(as.factor(optimal_clustering), true_partition)
    cls_measures[r,s,3]<-entropy(as.factor(optimal_clustering), true_partition)
    
  }  
}
  
cls_measures_df<-data.frame()
for(s in 1:dim(cls_measures)[2]){
  cls_measures_df<-rbind(cls_measures_df, 
                         data.frame(cbind(value_AdjRandIndex=cls_measures[,s,1], 
                                          value_Purity=cls_measures[,s,2], 
                                          value_Entropy=cls_measures[,s,3]),
                                    Scenario=paste0(s))
  )
}

library(tidyverse)
cls_measures_df_longer <- cls_measures_df %>% 
  pivot_longer(cols  = c(value_AdjRandIndex, value_Purity, value_Entropy), 
               names_to = c(".value", "Clustering_Measure"),
               names_sep = "\\_")


cls_measures_df_longer$Model <- "DPM_CER"

# devtools::install_github("teunbrand/ggh4x")
library(ggh4x) 
# library(reporter)

pp<-ggplot(cls_measures_df_longer,
           xlab=",",
           aes(x=Model, y=value, fill=Clustering_Measure)) +   theme_bw() +
  theme(axis.text.x=element_blank(), 
        axis.title.x=element_blank(), 
        axis.title.y=element_blank(), legend.position="none",
        legend.title=element_blank() ) +
  geom_violin( )  + 
  facet_grid2(Clustering_Measure~Scenario, 
              labeller = labeller(Clustering_Measure = c("Entropy" = "Entropy",
                                                         "Purity" = "Purity",
                                                         "AdjRandIndex" = "ARI"),
                                  Scenario = c("1"= "low",
                                               "2"= "medium-low",
                                               "3"= "medium",
                                               "4"= "high", 
                                               "5"= "mixed"
                                  ) 
              )
  ) +
  theme( strip.text.x = element_text(size = 15),
         strip.text.y.right  = element_text(size = 15),
         axis.text=element_text(size=15) ) 
pp

ggsave(paste0(input_path, "/clust_dist.pdf"), 
       plot = pp, width = 30, height = 15, units = "cm", dpi = 300)






