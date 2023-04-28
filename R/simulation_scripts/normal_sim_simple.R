# Script to show difference LCOV-COV and COV-COV4

data <- mixture_sim(pct_clusters = c(0.5, 0.5),
                    n = 1000,
                    p = 10,
                    delta = 10)

clusters <- data[,1]
data <- data[,-1]

ics_lcov <- ICSClust(data, criterion = "med_crit",
                    ICS_args = list(S1 = ICS_lcov, S2 = ICS_cov,
                                    S1_args = list(mscatter = "cov", proportion = 0.1)),
                    nb_clusters = 2,
                    nb_select = 1)

ics_cov4 <- ICSClust(data, criterion = "med_crit",
                     ICS_args = list(S1 = ICS_cov, S2 = ICS_cov4),
                     nb_clusters = 2,
                     nb_select = 1)

selected_lcov <- ics_lcov$select
selectted_cov4 <- ics_cov4$select

# Plot true clusters
GGally::ggpairs(data, aes(colour = factor(clusters)))


# Plot the IC's for LCOV-COV
component_plot(ics_lcov$ICS_out, clusters = factor(clusters))


# Plot the IC's for COV-COV4
component_plot(ics_cov4$ICS_out, clusters = factor(clusters))