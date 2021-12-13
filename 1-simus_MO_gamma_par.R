
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Simulation study for Marshall-Olkin copulas
#
# Influence of the parameter "gamma" and "truePar"
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# Loading of the libraries =============================================

library(VineCopula)
library(MMDCopula)
library(pbapply) # to have nice progressbars


# Creation of the file to store the results ==========================

model = "gamma_par"
fileName = paste0(c("Sim_MO", model), collapse = "_")

# Folder to store the results of the simulations
# To be modified in order to store results elsewhere
my_filepath = here::here()
my_filename = paste0(my_filepath, "/", fileName ,".csv")

if (! file.exists(my_filename)){
  write.table(
    x = paste0(c("n", "family", "truePar",
                 "typeOutliers", "nOutliers", "alpha_estimation",
                 "niter", "ndrawings", "naveraging",
                 "par.start",
                 "gamma", "nameEstimator", "estimatedPar",
                 "userTime", "systemTime", "realTime"), collapse=";") ,

    file = my_filename,
    append = F, sep = ";", col.names = FALSE, row.names = FALSE,
    quote = FALSE)
}

vecNameEstimators = c("", "Itau",
                      "MMD_gaussian", "MMD_gaussian.Phi",
                      "MMD_exp-l2", "MMD_exp-l2.Phi",
                      "MMD_exp-l1", "MMD_exp-l1.Phi" )

# Parameters of the experiment =======================================
NReplications = 200
nSample = 1000
family = "MO"
q = 0.001  # contamination parameter

# Parameters of the estimation
alpha_estimation = 1
niter = 100
ndrawings = 10
naveraging = 1

vecPossibleGammas = c(
  0.05, 0.07, 0.09, 0.1, 0.12, 0.15, 0.18, 0.2, 0.25, 0.3, 0.35,
  0.4, 0.6, 0.8, 0.9, 1, 1.1, 1.2, 1.4, 1.6, 1.8, 2)

vecPossiblePars = seq(0.2, 0.9, by = 0.1)


# Main loop ==========================================================
pb = pbapply::startpb(min = 0, max = 2 * NReplications *
                        length(vecPossibleGammas) * length(vecPossiblePars))
i_pb = 1

for (iRep_ in 1:NReplications)
{
  for (quantity_outliers in c(0, 5 * nSample / 100))
  {
    for (truePar in vecPossiblePars)
    {
      # Simulation =====================================================

      U <- BiCopSim.MO(n = nSample, alpha = truePar)

      typeOutliers = paste0("top_left", q, collapse = "_")
      if (quantity_outliers > 0) {
        U[1:quantity_outliers, 1] = runif(n = quantity_outliers, min = 0, max = q)
        U[1:quantity_outliers, 2] = runif(n = quantity_outliers, min = 1-q, max = 1)
      }

      # Estimation =====================================================

      # Computation of the pseudo-observations
      U = pobs(U)

      # Itau
      time1 = proc.time()
      resultItau <- BiCopEst.MO(u1 = U[,1], u2 = U[,2], method = "itau")
      time2 = proc.time()
      timeItau = (time2 - time1)[1:3]

      # Storing results of the Itau ============================

      vec_params = c(nSample, family, truePar,
                     typeOutliers, quantity_outliers,
                     NA, NA, NA, NA, NA
      )

      toWrite2 = c(vec_params, NA, vecNameEstimators[2],
                   resultItau$par, as.numeric(timeItau))

      write.table(
        x = rbind(toWrite2), file = my_filename ,
        append = T, sep = ";", col.names = FALSE, row.names = FALSE
      )

      # Storing results of the MMD =====================================

      for (i_gamma in 1:length(vecPossibleGammas))
      {
        gamma = vecPossibleGammas[i_gamma]

        for (nameKernel in c("gaussian"))
        {
          par.start = runif(n = 1, min = 0.05, max = 0.95)

          # MMD without Phi
          time1 = proc.time()
          nameKernel1 = nameKernel
          resultMMD_gaussian1 <- BiCopEst.MO(u1 = U[,1], u2 = U[,2], method = "MMD",
                                             kernel = nameKernel1, gamma = gamma,
                                             par.start = par.start, alpha = alpha_estimation,
                                             niter = niter, ndrawings = ndrawings,
                                             naveraging = naveraging
          )
          time2 = proc.time()
          timeMMD_gaussian1 = (time2 - time1)[1:3]

          # MMD with Phi
          time1 = proc.time()
          nameKernel2 = paste0(nameKernel, ".Phi")
          resultMMD_gaussian2 <- BiCopEst.MO(u1 = U[,1], u2 = U[,2], method = "MMD",
                                             kernel = nameKernel2, gamma = gamma,
                                             par.start = par.start, alpha = alpha_estimation,
                                             niter = niter, ndrawings = ndrawings,
                                             naveraging = naveraging
          )
          time2 = proc.time()
          timeMMD_gaussian2 = (time2 - time1)[1:3]


          # Saving the results of the experiments =================================
          vec_params = c(nSample, family, truePar,
                         typeOutliers, quantity_outliers, alpha_estimation,
                         niter, ndrawings, naveraging, par.start )

          toWrite3 = c(vec_params, gamma, paste("MMD", nameKernel1 , sep = '_'),
                       resultMMD_gaussian1$par, as.numeric(timeMMD_gaussian1))

          toWrite4 = c(vec_params, gamma, paste("MMD", nameKernel2 , sep = '_'),
                       resultMMD_gaussian2$par, as.numeric(timeMMD_gaussian2))

          write.table(
            x = rbind( toWrite3, toWrite4) , file = my_filename ,
            append = T, sep = ";", col.names = FALSE, row.names = FALSE
          )


          setpb(pb, i_pb)
          i_pb = i_pb + 1
        }
      }
    }
  }
}

pbapply::closepb(pb)

