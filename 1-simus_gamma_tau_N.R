
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Simulation study
#
# Influence of the parameter "gamma" and "tau"
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# Loading of the libraries =============================================

library(VineCopula)
library(MMDCopula)
library(pbapply) # to have nice progressbars


# Creation of the file to store the results ==========================

model = "gamma_tau_N"
fileName = paste0(c("Sim", model), collapse = "_")


# Folder to store the results of the simulations
# To be modified in order to store results elsewhere
my_filepath = here::here()
my_filename = paste0(my_filepath, "/", fileName ,".csv")

if (! file.exists(my_filename)){
  write.table(
    x = paste0(c("n", "family", "trueTau",
                 "typeOutliers", "nOutliers",
                 "gamma", "nameEstimator", "methodTauInit", "tauInit",
                 "estimatedKendallTau", "userTime", "systemTime", "realTime"), collapse=";") ,

    file = my_filename,
    append = F, sep = ";", col.names = FALSE, row.names = FALSE,
    quote = FALSE)
}

vecNameEstimators = c("MLE", "Itau",
                      "MMD_gaussian", "MMD_gaussian.Phi",
                      "MMD_exp-l2", "MMD_exp-l2.Phi",
                      "MMD_exp-l1", "MMD_exp-l1.Phi" )

# Parameters of the experiment =======================================
NReplications = 200
nSample = 1000
family = 1
vecPossibleTaus = seq(-0.9, 0.9, by = 0.1)

q = 0.001  # contamination parameter


vecPossibleGammas = c(
  0.05, 0.07, 0.09, 0.1, 0.12, 0.15, 0.18, 0.2, 0.25, 0.3, 0.35,
  0.4, 0.6, 0.8, 0.9, 1, 1.1, 1.2, 1.4, 1.6, 1.8, 2)
vecNamesKernel = c("gaussian", "exp-l2", "exp-l1")
vecMethodsTauInit = c("unif(-0.95;0.95)", "empKT")



# Main loop ==========================================================
pb = pbapply::startpb(min = 0, max = 2 * length(vecPossibleTaus) *
                        NReplications * length(vecPossibleGammas) *
                        length(vecNamesKernel) * length(vecMethodsTauInit))

i_pb = 1

for (iRep_ in 1:NReplications)
{
  for (trueTau in vecPossibleTaus)
  {
    for (nOutliers in c(0, 5 * nSample / 100))
    {

      # Simulation =====================================================

      U <- BiCopSim(N = nSample, family = family,
                    par = BiCopTau2Par(family = family, tau = trueTau))

      typeOutliers = paste0("top_left", q, collapse = "_")
      if (nOutliers > 0) {
        U[1:nOutliers, 1] = runif(n = nOutliers, min = 0, max = q)
        U[1:nOutliers, 2] = runif(n = nOutliers, min = 1-q, max = 1)
      }

      vec_params = c(nSample, family, trueTau, typeOutliers, nOutliers)


      # Estimation =====================================================

      # Computation of the pseudo-observations
      U = pobs(U)

      # MLE
      time1 = proc.time()
      resultMLE <- BiCopEst(u1 = U[,1], u2 = U[,2],
                            family = family, method = "mle")$tau
      time2 = proc.time()
      timeMLE = (time2 - time1)[1:3]

      # Itau
      time1 = proc.time()
      resultItau <- BiCopEst(u1 = U[,1], u2 = U[,2],
                             family = family, method = "itau")$tau
      time2 = proc.time()
      timeItau = (time2 - time1)[1:3]


      # Storing results of the MLE and Itau ============================

      toWrite1 = c(vec_params, NA,  vecNameEstimators[1], NA, NA,
                   resultMLE, as.numeric(timeMLE))

      toWrite2 = c(vec_params, NA, vecNameEstimators[2], NA, NA,
                   resultItau, as.numeric(timeItau))

      write.table(
        x = rbind(toWrite1, toWrite2), file = my_filename ,
        append = T, sep = ";", col.names = FALSE, row.names = FALSE
      )

      # Storing results of the MMD =====================================

      for (i_gamma in 1:length(vecPossibleGammas))
      {
        gamma = vecPossibleGammas[i_gamma]

        for (nameKernel in vecNamesKernel)
        {
          for (methodTauInit in vecMethodsTauInit)
          {
            tauInit = switch (
              methodTauInit,

              "unif(-0.95;0.95)" = runif(1, min = -0.95, max = 0.95),

              "unif(0.05;0.95)" = runif(1, min = 0.05, max = 0.95),

              "empKT" = pcaPP::cor.fk(U[,1], U[,2])
            )

            # MMD without Phi
            time1 = proc.time()
            nameKernel1 = nameKernel
            resultMMD_gaussian1 <- BiCopEstMMD(
              u1 = U[,1], u2 = U[,2], family = family,
              kernel = nameKernel1, gamma = gamma, tau = tauInit)$tau
            time2 = proc.time()
            timeMMD_gaussian1 = (time2 - time1)[1:3]

            # MMD with Phi
            time1 = proc.time()
            nameKernel2 = paste0(nameKernel, ".Phi")
            resultMMD_gaussian2 <- BiCopEstMMD(
              u1 = U[,1], u2 = U[,2], family = family,
              kernel = nameKernel2, gamma = gamma, tau = tauInit)$tau
            time2 = proc.time()
            timeMMD_gaussian2 = (time2 - time1)[1:3]


            # Saving the results of the experiments =================================

            toWrite3 = c(vec_params, gamma, paste("MMD", nameKernel1 , sep = '_'),
                         methodTauInit, tauInit,
                         resultMMD_gaussian1, as.numeric(timeMMD_gaussian1))

            toWrite4 = c(vec_params, gamma, paste("MMD", nameKernel2 , sep = '_'),
                         methodTauInit, tauInit,
                         resultMMD_gaussian2, as.numeric(timeMMD_gaussian2))

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
}

pbapply::closepb(pb)
