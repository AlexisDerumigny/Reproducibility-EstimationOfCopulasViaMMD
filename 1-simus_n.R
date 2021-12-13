
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Simulation study
#
# Influence of the parameters "n" and "nOutliers"
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# Loading of the libraries =============================================

library(VineCopula)
library(MMDCopula)
library(pbapply) # to have nice progressbars


# Creation of the file to store the results ==========================

model = "n"
fileName = paste0(c("Sim", model), collapse = "_")


# Folder to store the results of the simulations
# To be modified in order to store results elsewhere
my_filepath = here::here()
my_filename = paste0(my_filepath, "/", fileName, ".csv")

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

# Parameters of the experiment =======================================
NReplications = 200
vec_nSample = c(50, 100, 200, 500, 700, 1000, 2000, 3000, 5000, 7000, 10000)
family = 1
trueTau = 0.5

# Contamination
q = 0.001
pctOutliers = c(0, 5, 10)

# Estimation
methodTauInit = "unif(-0.95;0.95)"
vecNameEstimators = c("MLE", "Itau",
                      "MMD_gaussian", "MMD_gaussian.Phi")
vecGamma = c(0.25, 0.8)


pb = pbapply::startpb(min = 0, max = NReplications * length(pctOutliers) * length(vec_nSample))
i_pb = 1

# Main loop ==========================================================

for (iRep in 1:NReplications)
{
  for (i_n in 1:length(vec_nSample))
  {
    for (i_outliers in 1:length(pctOutliers))
    {
      nSample = vec_nSample[i_n]
      quantity_outliers = pctOutliers[i_outliers] * nSample / 100

      # Simulation =========================================================

      U <- BiCopSim(N = nSample, family = family,
                    par = BiCopTau2Par(family = family, tau = trueTau))


      # Contamination ======================================================

      typeOutliers = paste0("top_left", q, collapse = "_")
      if (quantity_outliers > 0) {
        U[1:quantity_outliers, 1] = runif(n = quantity_outliers, min = 0, max = q)
        U[1:quantity_outliers, 2] = runif(n = quantity_outliers, min = 1-q, max = 1)
      }


      # Estimation =========================================================

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

      # Initialisation of the MMD methods
      tauInit = switch (
        methodTauInit,

        "unif(-0.95;0.95)" = runif(1, min = -0.95, max = 0.95),

        "unif(0.05;0.95)" = runif(1, min = 0.05, max = 0.95),

        "empKT" = pcaPP::cor.fk(U[,1], U[,2])
      )

      # MMD_gaussian
      time1 = proc.time()
      resultMMD_gaussian1 <- BiCopEstMMD(u1 = U[,1], u2 = U[,2], family = family,
                                         kernel = "gaussian",
                                         gamma = vecGamma[1], tau = tauInit)$tau
      time2 = proc.time()
      timeMMD_gaussian1 = (time2 - time1)[1:3]

      # MMD_gaussian
      time1 = proc.time()
      resultMMD_gaussian2 <- BiCopEstMMD(u1 = U[,1], u2 = U[,2], family = family,
                                         kernel = "gaussian.Phi",
                                         gamma = vecGamma[2], tau = tauInit)$tau
      time2 = proc.time()
      timeMMD_gaussian2 = (time2 - time1)[1:3]


      vec_params = c(nSample, family, trueTau,
                     typeOutliers, quantity_outliers)

      # Saving the results of the experiments =========================================

      toWrite1 = c(vec_params, NA,  vecNameEstimators[1], NA, NA,
                   resultMLE, as.numeric(timeMLE))
      toWrite2 = c(vec_params, NA, vecNameEstimators[2], NA, NA,
                   resultItau, as.numeric(timeItau))
      toWrite3 = c(vec_params, vecGamma[1], vecNameEstimators[3],
                   methodTauInit, tauInit,
                   resultMMD_gaussian1, as.numeric(timeMMD_gaussian1))
      toWrite4 = c(vec_params, vecGamma[2], vecNameEstimators[4],
                   methodTauInit, tauInit,
                   resultMMD_gaussian2, as.numeric(timeMMD_gaussian2))

      write.table(
        x = rbind(toWrite1, toWrite2, toWrite3, toWrite4) ,
        file = my_filename ,
        append = T, sep = ";", col.names = FALSE, row.names = FALSE
      )

      setpb(pb, i_pb)
      i_pb = i_pb + 1
    }
  }
}

closepb(pb)
