# Adapted from Max Welz, creates the correlation matrix
get_blockmatrix <- function(construct_size,
                            num_constructs,
                            within_correlation,
                            between_correlation) {
  # assigns unique code to each construct in matrix
  # diagonal is all one
  # within_correlation is a vector of correlations of items in each construct
  # between_correlations is the correlation of items in different constructs
  
  num_items <- construct_size * num_constructs
  out <- matrix(between_correlation, nrow = num_items, ncol = num_items)
  name <- rep(NA_character_, num_items)
  J <- construct_size
  
  for(i in 0:(num_constructs - 1)){
    
    interval <- c(1:construct_size) + i * construct_size
    # Set all correlations within constructs equal to given correlation
    out[interval, interval] <- within_correlation[i + 1]
    name[interval] <- paste0(paste0("X", i + 1), "_", 1:construct_size)
    
  } # FOR
  
  diag(out) <- 1.0 # diagonal is one
  colnames(out) <- rownames(out) <- name
  
  out
} # FUN

# Adapted from Max Welz
# K is number of answer categories
get_distribution <- function(type, K = 5)
{
  stopifnot(K %in% c(3, 5, 7, 9))
  pi <- numeric()
  
  if(identical(type, "centered"))
  {
    if(identical(K, 3))
    {
      pi <- c(0.35, 0.3, 0.35)
    } else if(identical(K, 5))
    {
      pi <- c(0.15, 0.2, 0.3, 0.2, 0.15)
    } else if(identical(K, 7))
    {
      pi <- c(0.05, 0.125, 0.2, 0.25, 0.2, 0.125, 0.05)
    } else{
      pi <- c(0.025, 0.05, 0.15, 0.175, 0.2, 0.175, 0.15, 0.05, 0.025)
    } # IF centered
    
  } else if(identical(type, "agree"))
  {
    if(identical(K, 3))
    {
      pi <- c(0.225, 0.325, 0.45)
    } else if(identical(K, 5))
    {
      pi <- c(0.1, 0.15, 0.2, 0.25, 0.3)
    } else if(identical(K, 7))
    {
      pi <- c(0.05, 0.075, 0.1, 0.125, 0.15, 0.225, 0.275)
    } else{
      pi <- c(0.025, 0.05, 0.075, 0.075, 0.1, 0.125, 0.15, 0.175, 0.225)
    } # IF agree
    
  } else if(identical(type, "polarized"))
  {
    if(identical(K, 3))
    {
      pi <- c(0.475, 0.05, 0.475)
    } else if(identical(K, 5))
    {
      pi <- c(0.3, 0.175, 0.05, 0.175, 0.3)
    } else if(identical(K, 7))
    {
      pi <- c(0.25, 0.15, 0.075, 0.05, 0.075, 0.15, 0.25)
    } else{
      pi <- c(0.24, 0.15, 0.05, 0.05, 0.02, 0.05, 0.05, 0.15, 0.24)
    } # IF polarized
    
  } else if(identical(type, "disagree")) 
  {
    pi <- rev(get_distribution("agree", K))
    
  } else stop("unknown type")
  
  return(pi)
  
} # FUN


get_baseprobs <- function(construct_size, num_clusters, probability_types,
                          num_likert) {
  # gets the baseprobs for all clusters in one matrix
  # probability_types: list containing for each cluster a vector of probability
  # types (1 for each construct).
  # Returns a list of baseprobs matrixes, one for each cluster
  
  # list containing for each cluster, a matrix of marginal probabilities
  baseprobs_list <- lapply(seq_len(num_clusters), function(cluster) {
    
    cluster_prob_types <- probability_types[[cluster]]
    
    construct_probs <- lapply(cluster_prob_types, function(prob_type) {
      # Get vector with marginal probabilities based on type
      probabilities <- get_distribution(type = prob_type,
                                        K = num_likert)
      
      # repeat the probabilities in a matrix for Simstudy
      probability_matrix <- matrix(rep(probabilities, each = construct_size),
                                   nrow = construct_size)
      
      probability_matrix
    })
    
    do.call(rbind, construct_probs)
  })
  
  baseprobs_list
}

simulate_clusters <- function(baseprob_list, clusters, Rho, trial=-1) {
  num_clusters <- length(baseprob_list)
  
  data_clusters <- lapply(seq_len(num_clusters), function(cluster) {
    baseprobs <- baseprob_list[[cluster]]
    n <- sum(clusters == cluster)
    temp <- simstudy::genData(n)
    data <- simstudy::genOrdCat(temp,
                                baseprobs = baseprobs,
                                prefix = "q",
                                corMatrix = Rho)
    

    data <- data[,-1] # drop id
    colnames(data) <- colnames(Rho)
    
    # add identifying variables
    data$trial <- trial
    data$cluster <- cluster
    
    names <- colnames(data)
    
    # Convert to integer matrix
    data <- matrix(as.integer(as.matrix(data)), nrow = n)

    # Reset names
    colnames(data) <- names
    

    data
  })
  do.call(rbind, data_clusters)
}
