
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Aggregation of the results for computations of the MSE
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# Loading of the libraries
library(tidyverse)


# Loading of the simulation data for absolutely continuous copulas =====================

# Models studied
vec_models = c("gamma_tau_N", "gamma_tau_fam",
               "typeContam", "contamFam", "n")

# Folder to store the results of the simulations
# To be modified in order to load results from elsewhere
my_filepath = here::here()

filenames = paste0(my_filepath, "/Sim_", vec_models, ".csv")

totalData = filenames %>%
  map( function(x){return(if(file.exists(x)){x} else {NULL})} ) %>%
  unlist() %>%
  map_dfr( read.csv, sep = ";", dec = ".", header = TRUE)

# Joining together all simulations without outliers
allTypeOutliers = unique(totalData$typeOutliers)
totalData[which(totalData$nOutliers == 0), "typeOutliers"] = "0"


# Aggregating all the results
totalAggregated = totalData %>%
  group_by(n, family, trueTau,
           typeOutliers, nOutliers,
           gamma, nameEstimator, methodTauInit) %>%
  summarise(MSE = mean((estimatedKendallTau - trueTau)^2, na.rm = TRUE),
            Bias = mean(estimatedKendallTau - trueTau, na.rm = TRUE),
            Var = var(estimatedKendallTau, na.rm = TRUE),
            Sd_MSE = sd((estimatedKendallTau - trueTau)^2, na.rm = TRUE),
            nRep = n(),
            MSE_q05 = MSE - 1.96 * Sd_MSE / nRep,
            MSE_q95 = MSE + 1.96 * Sd_MSE / nRep,
            meanRealtime = mean(realTime),
            sdRealtime = sd(realTime),
            .groups = "keep") %>%
  ungroup() %>%
  mutate(nameEstimator = if_else(nameEstimator == "MLE", "CML", nameEstimator))

totalAggregated_0 = totalAggregated %>% dplyr::filter(typeOutliers == "0")
totalAggregated = totalAggregated %>% dplyr::filter(typeOutliers != "0")

totalAggregated =
  bind_rows(
    totalAggregated,
    allTypeOutliers %>% map_dfr(
      function(x){
        totalAggregated_0 %>% mutate(typeOutliers = x)
      }
    )
)

# Loading of the simulation data for Marshall-Olkin copulas ==========================


# Models studied
vec_models = c("MO_gamma_par", "MO_nOutliers")

filenames = paste0(my_filepath, "/Sim_", vec_models, ".csv")

totalData_MO = filenames %>%
  map( function(x){return(if(file.exists(x)){x} else {NULL})} ) %>%
  unlist() %>%
  map_dfr(read.csv, sep = ";", dec = ".", header = TRUE)

# Aggregating all the results

totalAggregated_MO = totalData_MO %>%
  group_by(n, family, truePar,
           typeOutliers, nOutliers,
           alpha_estimation, niter, ndrawings, naveraging,
           gamma, nameEstimator) %>%
  summarise(MSE = mean((estimatedPar - truePar)^2, na.rm = TRUE),
            Bias = mean(estimatedPar - truePar, na.rm = TRUE),
            Var = var(estimatedPar, na.rm = TRUE),
            Sd_MSE = sd((estimatedPar - truePar)^2, na.rm = TRUE),
            nRep = n(),
            MSE_q05 = MSE - 1.96 * Sd_MSE / nRep,
            MSE_q95 = MSE + 1.96 * Sd_MSE / nRep,
            meanRealtime = mean(realTime),
            sdRealtime = sd(realTime),
            .groups = "keep") %>%
  ungroup()
