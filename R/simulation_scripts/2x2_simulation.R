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
seed <- 42                    

# number of simulation runs
B <- 100

# sample sizes
n <- c(1000)
n_max <- max(n)

# sizes of construct 
construct_size <- 6

# number of constructs
num_constructs <- 2

# correlation between items in the same construct
within_corr <- 0.7

# list containing for each cluster a list of probability types
# (1 per construct)
probability_types <- list(cluster1 = list("agree", "agree"),
                          cluster2 = list( "disagree","disagree"),
                          cluster3 = list("agree", "disagree"),
                          cluster4 = list("disagree", "agree"))

# number of clusters
nb_clusters <- length(probability_types)

# probability of belonging to particular cluster
cluster_probs <- list(c(0.1, 0.1, 0.1, 0.7))

# number of categories per item (K))
num_likert <- 7

# Types of second scatter matrix
scatter_list <- list("cov4" = ICS_cov4,
                     "tcov" = ICS_tcov,
                     "lcov" = ICS_lcov)

# define the correlation matrix of the class populations
Rho <- get_blockmatrix(construct_size = construct_size,
                       num_constructs = num_constructs,
                       within_correlation = rep(within_corr, num_constructs),
                       between_correlation = 0.2)

set.seed(seed)

# Vary over the different probabilities of being in a certain cluster
results_probs <- lapply(cluster_probs, function (cluster_prob) { 
  
  results_scatter <- lapply(seq_len(length(scatter_list)), function (i) {
    scatter <- scatter_list[[i]]
    scatter_name <- names(scatter_list)[i]
    
    print(paste("Trial with", scatter_name, "and probabilities of",
                cluster_probs))
    
    
    # vary simulation runs
    results_r <- lapply(seq_len(B), function(r) {
      print(paste("Trial", r))

      # Generate the clusters
      clusters <- sample(seq_len(length(cluster_prob)), replace = TRUE,
                         prob = cluster_prob, size = n)
      
      # number of clusters in actual data
      nb_clusters <- length(cluster_prob)
      
      # order the cluster labels to make it easy to simulate    
      clusters <- sort(clusters)
      
      # get a list of baseprobs matrixes for all clusters to pass to Simstudy
      baseprob_list <- get_baseprobs(construct_size = construct_size,
                                     nb_clusters = nb_clusters,
                                     probability_types = probability_types,
                                     num_likert = num_likert)
      
      # generates the data for all clusters in a single matrix
      full_data <- generate_clusters(baseprob_list = baseprob_list,
                                     clusters = clusters,
                                     Rho = Rho, trial = r)
      
      # Remove trial and cluster information
      data <- full_data[, seq_len(num_constructs * construct_size)]
      
      # Perform ICS
      ICS_out <- ICSClust(data, nb_select = 2,
                          nb_clusters = nb_clusters,
                          method = 'kmeans_clust', criterion = 'med_crit',
                          ICS_args = list(S1 = scatter, S2 = ICS_cov))
      
      # Compute eta squared
      eta2 <- eta2_power(object = ICS::components(ICS_out$ICS_out), clusters = clusters,
                         select = ICS_out$select)
      
      # Compute adjusted Rand Index
      ARI  <- mclust::adjustedRandIndex(clusters, ICS_out$clusters)
      
      # Store found clusters
      full_data <- cbind(full_data, ICS_cluster = ICS_out$clusters, eta2 = eta2,
                         ARI = ARI, scatter = scatter_name)
      full_data
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


# Take only results from trial 1 to make plots
results_1 <- df_res %>%
  filter(trial == 1, scatter == "cov4") %>%
  dplyr::select(-c(trial))

# get true clusters from data
clusters <- factor(results_1$true_cluster)

results_1 <- results_1 %>% 
  select(-c(ICS_cluster, true_cluster, eta2, ARI, scatter))

# Plotting all pairs
GGally::ggpairs(results_1, aes(colour = factor(clusters)),
                title = "True clusters")

# Performing ICS on simulated data
ICS_components <- ICS(results_1, S1 = ICS_cov, S2 = ICS_cov4)

# Plot IC's
component_plot(ICS_components, clusters = clusters)

# Performing ICSClust on simulated data
ICS_out <- ICSClust(results_1, nb_select = nb_clusters - 1,
                    nb_clusters = nb_clusters,
                    method = 'kmeans_clust',
                    criterion = 'med_crit')

# Plot selected components
component_plot(ICS_out$ICS_out, select = ICS_out$select,
               clusters = as.factor(ICS_out$clusters))
