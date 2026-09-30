library(DPmCER)

project_path <- paste0(sub(as.character(desc::desc_get("Package")), "", getwd()), "storing_folders")

if (!dir.exists( paste0(project_path, "/simulations") ) ) {
  dir.create( paste0(project_path, "/simulations") )
}

sim_path <- paste0(project_path, "/simulations/sim_core_periphery")

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

N<-20
# number of edges
E<-choose(N, 2)

# sample size
n<-40
# number of clusters
k<-2
# cluster proportions in the first k-1 blocks
w<-rep(1/k, k-1)
# clusters size
n_j<-c(ceiling(n*w), n-sum(ceiling(n*w)) )

true_partition<-as.numeric(unlist(sapply(1:length(n_j), function(j){rep(j, n_j[j])})))

library(igraph)
# True parameters 
alpha_1 <- 0.4
alpha_2 <- 0.3

# Gm_j
set.seed(1)
# Core-periphery
n_total <- N
n_core <- 5
n_periphery <- N - n_core

# Define group memberships: 1 = core, 2 = periphery
block_membership <- c(rep(1, n_core), rep(2, n_periphery))

# Define SBM connection probability matrix (non-assortative pattern)
# Rows and columns correspond to group indices (1 = core, 2 = periphery)
core2core <- 0.8
core2peri <- peri2core <- 0.3
peri2peri <- 0.05
pref_matrix <- matrix(c(
  core2core, core2peri,  # Core to Core, Core to Periphery
  peri2core, peri2peri   # Periphery to Core, Periphery to Periphery
), nrow = 2, byrow = TRUE)

# Simulate SBM graph
G_m_1 <- sample_sbm(n = n_total, pref.matrix = pref_matrix, 
                    block.sizes = c(n_core, n_periphery),
                    directed = FALSE, loops = FALSE )

# Plotting
V(G_m_1)$color <- ifelse(block_membership == 1, "skyblue", "orange")
V(G_m_1)$label.cex = 25
E(G_m_1)$width = 10

lay <- layout_with_fr(G_m_1)

# plot(G_m_1, vertex.label.color="blue", 
#      vertex.color=V(G_m_1)$color,
#      edge.color=E(G_m_1)$color,
#      mode="undirected",
#      layout = lay)


# # EXAMPLE OF OBSERVATION FROM THE FIRST COMPONENT
plot(graph_from_adjacency_matrix(from_low_adj_to_full_adj(
  sample_from_CER(from_full_adj_to_low_adj(as.matrix(as_adjacency_matrix(G_m_1))), alpha_1, 1)),
  mode="undirected") )


# ER with prob core2core (just sample from a CER with null graph mode so that 
# each G_ij^m ~ Ber(core2core) )
set.seed(2)
G_m_2 <- graph_from_adjacency_matrix(from_low_adj_to_full_adj(
  sample_from_CER( rep(0, E), core2core, 1)), mode="undirected")
V(G_m_2)$label.cex = 25
E(G_m_2)$width = 10

# plot(G_m_2, vertex.label.color="blue",
#      vertex.color=V(G_m_2)$color,
#      edge.color=E(G_m_2)$color,
#      mode="undirected",
#      layout = lay)


# # EXAMPLE OF OBSERVATION FROM THE SECOND COMPONENT
plot(graph_from_adjacency_matrix(from_low_adj_to_full_adj(
  sample_from_CER(from_full_adj_to_low_adj(as.matrix(as_adjacency_matrix(G_m_2))), alpha_2, 1)),
  mode="undirected") )

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

st_s<-Sys.time()
foreach( r = 1:N_rep, .packages = c("DPMCER", "igraph") ) %dopar% {
  
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
  st <- Sys.time()
  dpm_cer_inference = DPM_CER(G_mat, G_0, a, b, c, sample, burn, thinning)
  en <- Sys.time()
  en-st
  
  out <- list(input=list(data=G_mat), output=dpm_cer_inference)
  
  saveRDS(out, paste0(output_path, "/output_rep_", r, ".rds"))
  rm( out, dpm_cer_inference, G_0, G_mat, data_arr )
  gc()
}
stopCluster(myCluster)
en_s<-Sys.time()
en_s-st_s

library(salso)
library(fossil) # rand.index
library(NMF) # purity and entropy of a clustering
load(paste0(input_path, "/common_env.RData"))

out_files <- list.files(output_path)

cls_measures <- matrix(NA, nrow = length(out_files), ncol = 3)
for(r in 1:length(out_files)){
  
  out_list <- readRDS( paste0(output_path, "/", out_files[r] ) )
  
  dpm_cer_inference = out_list$output
  
  PARTITION <- t(sapply(
    1:length(dpm_cer_inference),
    function(t)
      dpm_cer_inference[[t]]$rho
  ))
  
  optimal_clustering <- as.numeric(salso(PARTITION, loss=VI(a=1), 
                                         nRuns=10, nCores=0) )
  
  cls_measures[r,1]<-adj.rand.index(true_partition, optimal_clustering)
  cls_measures[r,2]<-purity(as.factor(optimal_clustering), true_partition)
  cls_measures[r,3]<-entropy(as.factor(optimal_clustering), true_partition)
  print(r)
}

cls_measures_df<-data.frame(cbind(value_AdjRandIndex=cls_measures[,1], 
                                  value_Purity=cls_measures[,2], 
                                  value_Entropy=cls_measures[,3]))

saveRDS(cls_measures_df, file= paste0(input_path, "/cls_measures_df.rds" )) 

cls_measures_df = readRDS(paste0(input_path, "/cls_measures_df.rds" ))

library(tidyverse)
cls_measures_df_longer <- cls_measures_df %>% 
  pivot_longer(cols  = c(value_AdjRandIndex, value_Purity, value_Entropy), 
               names_to = c(".value", "Clustering_Measure"),
               names_sep = "\\_")

cls_measures_df_longer$Model="DPM_CER"

pp <- ggplot(cls_measures_df_longer,  xlab=",",
             aes(x=Clustering_Measure, y=value, fill=Clustering_Measure)) +
  geom_violin(alpha = 0.5, bw = 0.015) + 
  scale_fill_manual(name="Clustering_Measure", values=c( rgb(0.9290, 0.6940, 0.1250),
                                                         4,
                                                         rgb(0, 0, 0, 0.5)),
  ) +
  theme_bw() +
  theme(axis.text.x=element_blank(), axis.title.x=element_blank(),
        axis.text=element_text(size=20),
        axis.title.y=element_blank(), legend.position="none",
        legend.title=element_blank() ) + 
  theme(aspect.ratio = 1) 

pp  

ggsave(paste0(input_path, "/core_peripery_cls.pdf"),
       plot = pp, width = 30, height = 15, units = "cm", dpi = 300)
