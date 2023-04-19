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


get_baseprobs <- function(construct_size, nb_clusters, probability_types,
                          num_likert) {
  # gets the baseprobs for all clusters in one matrix
  # probability_types: list containing for each cluster a vector of probability
  # types (1 for each construct).
  # Returns a list of baseprobs matrixes, one for each cluster
  
  # list containing for each cluster, a matrix of marginal probabilities
  baseprobs_list <- lapply(seq_len(nb_clusters), function(cluster) {
    
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

generate_clusters <- function(baseprob_list, clusters, Rho, trial=-1) {
  nb_clusters <- length(unique(clusters))
  
  data_clusters <- lapply(seq_len(nb_clusters), function(cluster) {
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
    data$true_cluster <- cluster
    
    names <- colnames(data)
    
    # Convert to integer matrix
    data <- matrix(as.integer(as.matrix(data)), nrow = n)
    
    # Reset names
    colnames(data) <- names
    
    data
  })
  do.call(rbind, data_clusters)
}

create_confusion_matrix <- function(predicted_clusters, true_clusters) {
  # function that takes in a vector predicted clusters and true clusters
  # returns the confusion matrix that maximizes correct classification
  
  predicted_clusters <- as.numeric(predicted_clusters)
  true_clusters <- as.numeric(true_clusters)
  N <- max(max(predicted_clusters), max(true_clusters))
  
  # create an empty matrix to store the counts
  count_matrix <- matrix(0, nrow = N, ncol = N)
  
  # fill the matrix with the frequency of each pair of true and predicted clusters
  for (i in 1:length(predicted_clusters)) {
    count_matrix[true_clusters[i], predicted_clusters[i]] <- count_matrix[true_clusters[i], predicted_clusters[i]] + 1
  }
  
  # perform optimal assignment using the Hungarian algorithm
  assignment <- solve_LSAP(count_matrix, maximum = TRUE)
  
  # create a confusion matrix using the optimal assignment
  confusion_matrix <- t(count_matrix[assignment,])
  rownames(confusion_matrix) <- 1:N
  colnames(confusion_matrix) <- 1:N
  
  return(confusion_matrix)
}

# From ICSClust package
# https://github.com/AuroreAA/ICSClust
# TODO: Remove when function is exported
# Object is the ICS data, clusters are all the IC's, select are the IC's selected
eta2_power <- function(object, clusters, select){
  if(is.null(clusters)){
    warning("The 'clusters' argument is mandatory to compute the discriminatory 
            power of the reduced data frame.")
  }else{
    df <- data.frame(clusters = as.factor(clusters), object[, select])
    
    # Univariate case: ANOVA
    if(length(select) == 1){
      
      ICS_mod <- lm(as.formula(paste("cbind(", colnames(df)[-1], ") ~ clusters")),
                    data = df)
      heplots::etasq(ICS_mod)[1,1]
      
    }else{
      # Multivariate case: MANOVO with Wilks test
      
      ICS_mod <- manova(as.formula(paste("cbind(", 
                                         paste(colnames(df)[-1], collapse = ","), 
                                         ") ~ clusters")), data = df)
      heplots::etasq(ICS_mod, test = "Wilks")[1,1]
    }
  }
}


