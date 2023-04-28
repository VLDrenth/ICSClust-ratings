# simulations for ICS with discrete data

# libraries
library(ICS)
library(ICSClust)
library(simstudy)
library(dplyr)
library(ggplot2)
library(GGally)
library(magrittr)
library(clue)

# load helper functions
source("R/utils.R")

# Seed of rng for simulation
seed <- 123                  

# number of simulation runs
B <- 50

# sample size
n <- 1000
p <- 11
delta <- 10

# number of clusters
nb_clusters <- 2

# probability of belonging to particular cluster
cluster_probs <- list(c(0.5, 0.5))

scatter_list <- list(
  `COV-COV4` = list(S1 = ICS_cov, S2 = ICS_cov4),
  `LCOV-COV` = list(S1 = ICS_lcov, S2 = ICS_cov,
                    S1_args = list(mscatter = "cov", proportion = 0.1)),
  `TCOV-COV` = list(S1 = ICS_tcov, S2 = ICS_cov,
                    S1_args = list(beta = 2)))

set.seed(seed)

# Vary over the different probabilities of being in a certain cluster
results_probs <- lapply(cluster_probs, function (cluster_prob) { 
  
  results_scatter <- lapply(seq_len(length(scatter_list)), function (i) {
    scatter_pair <- scatter_list[[i]]
    scatter_name <- names(scatter_list)[i]
    
    print(paste("Trial with", scatter_name, "and probabilities of",
                cluster_probs))
    
    
    # vary simulation runs
    results_r <- lapply(seq_len(B), function(r) {
      print(paste("Trial", r))
      
      # For each cluster generate data according to its covariance
      data <- ICSClust::mixture_sim(pct_clusters = cluster_prob,
                                    n = n,
                                    p = p,
                                    delta = delta
      )
      
      # Remove clusters from data
      clusters <- data[,1]
      data <- data[,-1]
      
      # Remove 'group' from the cluster column
      clusters <- lapply(clusters, function(x) {
        number <- stringr::str_sub(x, start = -1)
        as.integer(number)
      })
      
      clusters <- as.numeric(clusters)
      
      # Perform ICS
      ICS_out <- ICSClust(data,
                          nb_select = nb_clusters - 1,
                          nb_clusters = nb_clusters,
                          method = 'kmeans_clust',
                          criterion = 'med_crit',
                          ICS_args = scatter_pair,
                          clustering_args = list(iter.max = 100, nstart = 20)
      )
      
      # Compute eta squared
      eta2 <- eta2_power(object = ICS::components(ICS_out$ICS_out),
                         clusters = clusters,
                         select = ICS_out$select)
      
      # Compute adjusted Rand Index
      ARI  <- mclust::adjustedRandIndex(clusters, ICS_out$clusters)
      
      # Store found clusters
      data <- cbind(data, ICS_cluster = ICS_out$clusters, eta2 = eta2,
                    ARI = ARI, scatter = scatter_name)
      data
    })
    do.call(rbind, results_r)
  })
  do.call(rbind, results_scatter)
})

results <- do.call(rbind, results_probs)


df_res <- as.data.frame(results)

df_res[,seq_len(ncol(df_res) - 1)] <- 
  sapply(df_res[,seq_len(ncol(df_res) - 1)], as.numeric)

df_res %>% group_by(scatter) %>%
  summarize(ARI_mean = mean(ARI),
            eta2_mean = mean(eta2))
