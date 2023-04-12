# simulations for ICS with discrete data

# libraries
library(ICS)
library(ICSClust)
library(simstudy)
library(dplyr)
library(ggplot2)
library(GGally)
library(magrittr)

# load helper functions
source("R/utils.R")

# Seed of rng for simulation
seed <- 42                    

# number of simulation runs
B <- 5

# sample sizes
n <- c(10000)
n_max <- max(n)

# probability of belonging to particular cluster
cluster_probs <- list(c(0.5, 0.5))

# sizes of construct 
construct_size <- 20

# number of constructs
num_constructs <- 2

# list containing for each cluster a list of probability types
# (1 per construct)
probability_types <- list(list("centered", "agree"),
                          list("centered", "disagree"))

# number of categories per item (K)
num_likert <- 7

# define the correlation matrix of the class populations
Rho <- get_blockmatrix(construct_size = construct_size,
                       num_constructs = num_constructs,
                       within_correlation = c(0.8, 0.8),
                       between_correlation = 0.2)

set.seed(seed)

# Vary over the different probabilities of being in a certain cluster
results_probs <- lapply(cluster_probs, function (cluster_prob) { 
  
  # vary simulation runs
  results_r <- lapply(seq_len(B), function(r) {
    
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
    data <- generate_clusters(baseprob_list = baseprob_list,
                              clusters = clusters,
                              Rho = Rho, trial = r)
    data
  })
  do.call(rbind, results_r)
})

results <- do.call(rbind, results_probs)


# Visualizing the clusters
df_res <- as.data.frame(results)

results_1 <- df_res %>% filter(trial == 1) %>% dplyr::select(-c(trial))
clusters <- factor(results_1$cluster)
results_1 <- results_1 %>% select(-c(cluster))

# Plotting all pairs
GGally::ggpairs(results_1, aes(colour = factor(clusters)), title = "True clusters")

# Performing ICS on simulated data
ICS_out <- ICS(results_1)

# Plot IC's
component_plot(ICS_out, clusters = clusters)

# Performing ICSClust on simulated data
ICS_out <- ICSClust(results_1, nb_select = 1, nb_clusters = 2,
                    method = 'kmeans_clust', criterion = 'med_crit')

table(ICS_out$clusters, clusters)

component_plot(ICS_out$ICS_out, select = ICS_out$select,
               clusters = as.factor(ICS_out$clusters))

# Clustering k-means on the raw data
km <- kmeans(results_1, centers = 2)
table(km$cluster, clusters)

