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

# sample size 
n <- c(1000)
n_max <- max(n)

# probability of belonging to particular cluster
cluster_probs <- list(c(0.5, 0.5))

# sizes of construct 
construct_size <- 3

# number of constructs
num_constructs <- 2

# list containing for each cluster a vector of probability
# types (1 for each construct).
probability_types <- list(list("centered", "agree"),
                          list("centered", "disagree"))

# number of categories per item (K)
num_likert <- 7

# current possible distributions
marginals <- c(c("agree", "disagree"))

# define the correlation matrix for the simulated data

Rho <- get_blockmatrix(construct_size = construct_size,
                       num_constructs = num_constructs,
                       within_correlation = c(0.7, 0.7),
                       between_correlation = 0.2)

set.seed(seed)

# Vary over the different probabilities of being in a certain cluster
results_probs <- lapply(cluster_probs, function (cluster_prob) { 
  
  # vary simulation runs
  results_r <- lapply(seq_len(B), function(r) {
    # Generate the clusters
    clusters <- sample(seq_len(length(cluster_prob)), replace = TRUE,
                       prob = cluster_prob, size = n)
    
    # order the cluster labels to make it easy to simulate    
    clusters <- sort(clusters)
    
    baseprobs <- get_baseprobs(construct_size = construct_size,
                               num_clusters = length(cluster_prob),
                               probability_types = probability_types,
                               num_likert = num_likert)
    
    
    data <- simulate_clusters(baseprobs, clusters, Rho, trial = r)
    data
  })
  do.call(rbind, results_r)
})

results <- do.call(rbind, results_probs)


# Visualizing the cluster
df_res <- as.data.frame(results)

results_1 <- df_res %>% filter(trial == 1) %>% dplyr::select(-c(trial))
clusters <- results_1$cluster
results_1 <- results_1 %>% select(-c(cluster))

GGally::ggpairs(results_1, aes(colour = factor(clusters)), title = "True clusters")

