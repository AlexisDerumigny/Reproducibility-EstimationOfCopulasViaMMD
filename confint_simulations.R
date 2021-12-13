
# Preparation of the setting ===========================================

library(MMDCopula)
library(VineCopula)
library(pbapply)


my_settings = data.frame(
  nSample = rep(c(1000, 1000, 2000, 2000), times = 2),
  family = rep(c(1, 3), times = 4),

  methodTauInit = rep(c("unif(-0.95;0.95)",
                        "unif(0.05;0.95)" ), times = 4),
  gamma = rep(c(0.25, 0.1), times = 4),
  subsamplingSize = c(144, 144, 263, 263,
                      1000, 1000, 2000, 2000)
)


# Main loop ==========================================================
NReplications = 100
pb = startpb(min = 0, max = NReplications * nrow(my_settings))
i_pb = 1
setpb(pb, 0)

for (iSetting in 1:nrow(my_settings))
{

  # Creation of the file to store the results ==========================

  # Folder to store the results of the simulations
  # To be modified in order to store results elsewhere
  my_filepath = here::here()

  nameFile = paste(my_filepath, "/simus_CI_tau_", iSetting, ".csv", sep = "")

  if (!file.exists(nameFile)){

    write.table(
      x = paste0(c("n", "family", "trueTau",
                   "typeOutliers", "nOutliers",
                   "gamma", "nameEstimator", "methodTauInit", "tauInit",
                   "nResampling", "subsamplingSize", "corrSubSampling", "level",
                   "lowerCI.Tau", "upperCI.Tau", "lowerCI.Par", "upperCI.Par",
                   "userTime", "systemTime", "realTime"), collapse=";") ,

      file = nameFile ,
      append = F, sep = ";", col.names = FALSE, row.names = FALSE,
      quote = FALSE)
  }

  for (iReplications in 1:NReplications)
  {
    cat(" Replication ", iReplications, " out of ", NReplications, "...\n")

    ## Parameters of the simulation ====================================
    nSample = my_settings$nSample[iSetting]
    family = my_settings$family[iSetting]
    methodTauInit = my_settings$methodTauInit[iSetting]
    trueTau = 0.5


    ## Parameters of the estimation ====================================
    gamma1 = my_settings$gamma[iSetting]
    nameKernel1 = "gaussian"

    nResampling = 100
    subsamplingSize = my_settings$subsamplingSize[iSetting]
    corrSubSampling = TRUE
    level = 0.95


    ## Simulation ======================================================
    U <- BiCopSim(N = nSample, family = family,
                  par = BiCopTau2Par(family = family, tau = trueTau))


    ## Contamination ===================================================
    q = 0.001
    nOutliers = 0
    typeOutliers = paste0("top_left", q, collapse = "_")
    if (nOutliers > 0) {
      U[1:nOutliers, 1] = runif(n = nOutliers, min = 0  , max = q)
      U[1:nOutliers, 2] = runif(n = nOutliers, min = 1-q, max = 1)
    }


    ## Estimation ======================================================
    if (methodTauInit == "unif(-0.95;0.95)"){
      tauInit = runif(1, min = -0.95, max = 0.95)
    } else if (methodTauInit == "empKT") {
      tauInit = pcaPP::cor.fk(U[,1], U[,2])
    } else if (methodTauInit == "unif(0.05;0.95)"){
      tauInit = runif(1, min = 0.05, max = 0.95)
    }

    time1 = proc.time()
    resultMMD <- BiCopConfIntMMD(
      x1 = U[,1], x2 = U[,2], family = family,
      nResampling = nResampling, subsamplingSize = subsamplingSize,
      corrSubSampling = corrSubSampling , level = level,
      kernel = nameKernel1, gamma = gamma1, tau = tauInit)
    time2 = proc.time()
    timeMMD = (time2 - time1)[1:3]


    ## Saving the results ==============================================
    params_simu = c(nSample, family, trueTau)
    params_contam = c(typeOutliers, nOutliers)
    params_esti = c(gamma1, paste("MMD", nameKernel1 , sep = '_'), methodTauInit, tauInit)
    params_CI = c(nResampling, subsamplingSize, corrSubSampling, level)

    toWrite = c(
      params_simu, params_contam, params_esti, params_CI,
      resultMMD$CI.Tau[1], resultMMD$CI.Tau[2],
      resultMMD$CI.Par[1], resultMMD$CI.Par[2], as.numeric(timeMMD) )

    write.table(
      x = t(toWrite), file = nameFile ,
      append = T, sep = ";", col.names = FALSE, row.names = FALSE
    )

    setpb(pb, i_pb)
    i_pb = i_pb + 1
  }

}
