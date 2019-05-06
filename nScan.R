###### TE Simple Model: TE estimation on original & surrogate data ######
### version for cluster: nohub parallelization for scanning all pairs in the system & m values ###
### scanning parameters mx, my, mz individually/parallely with (un-)conditioned TE ### 

# rm(list = ls())

### where to run tests
## running locally: scanning only pairs, not m values
# testRun <- "local"
## running on cluster: scanning pairs and m values
testRun <- "cluster"
## nohup = T for paralleization via nohup, else F
nohup <- T

### m-Scan for kraskov estimator: which m ("mx","my","mz") shoud be scanned via nohup
# mScan <- "mx"
# mScan <- "my"
mScan <- "mz"
## choose "m" for all incrementing all history lengths parallely
# mScan <- "m"
## range of m-Scan
mScan_c <- 1:4

### range of dt-Scan for determining sampling frequency
dt <- 16
# N <- 5000
nScan <- c(50, 100, 200, 500, 1000, 2000, 4000, 8000, 16000, 32000, 64000, 128000)
N <- nScan[1]
dur_c <- N * dt

### parameters for Kraskov estimator
case_XY <- 3
k_kras <- 3

mx_c <- 2
my <- 2
mz <- 2

### model type: oscillation (os), steady-state(ss)
modelType <- "os"
# modelType <- "ss"

### parallel surrogate testing: True or False
surTest <- T



### list for all combinations of pairs
numSpec <- 8     # number of species in model
# numSpec <- 2     # number of species in model

pairLst <- list()
for (x in 2:(numSpec+1)) {
  for (y in 2:(numSpec+1)) {
    if (x != y)
      pairLst[[length(pairLst) + 1]] <- c(x,y)
  }
}

### parameter scan via nohup: variable 'ind' is scanning variabel given via nohup
if (nohup == T) {
  len_pairs <- length(pairLst)
  len_mScan <- length(mScan_c)
  # len_dtScan <- length(dtScan)
  
  condList <- list()
  for (i in 1:length(nScan)) {
    condList[[i]] <- c("mz", nScan[i], "cond", "os")
    # condList[[i*6-5]] <- c("mx", nScan[i], "uncond", "os")
    # condList[[i*6-4]] <- c("my", nScan[i], "uncond", "os")
    # condList[[i*6-3]] <- c("mz", nScan[i], "cond", "os")
    # condList[[i*6-2]] <- c("mx", nScan[i], "uncond", "ss")
    # condList[[i*6-1]] <- c("my", nScan[i], "uncond", "ss")
    # condList[[i*6-0]] <- c("mz", nScan[i], "cond", "ss")
    # condList[[i*3-2]] <- c("mx", nScan[i], "uncond", "os")
    # condList[[i*3-1]] <- c("my", nScan[i], "uncond", "os")
    # condList[[i*3-0]] <- c("mz", nScan[i], "cond", "os")
  }
  
  lenScan <- len_pairs * len_mScan
  id_cond <- floor((ind-1) / lenScan) + 1
  mScan <- condList[[id_cond]][1]
  N <- as.numeric(condList[[id_cond]][2])
  caseTE <- condList[[id_cond]][3]
  modelType <- condList[[id_cond]][4]
  
  if (caseTE == "cond")
    case_XY <- 3
  if (caseTE == "uncond")
    case_XY <- 1
  
  id_pairs <- ind %% len_pairs
  if (id_pairs == 0)
    id_pairs <- len_pairs
  pairLst <- pairLst[id_pairs]
  
  id_m <- floor((ind-1) / len_pairs) + 1
  id_m <- id_m %% len_mScan
  if (id_m == 0)
    id_m <- len_mScan
  
  if (mScan == "mx") {
    mx_c <- mScan_c[id_m]
    mm <- mx_c
  } else if (mScan == "my") {
    my <- mScan_c[id_m]
    mm <- my
  } else if (mScan == "mz") {
    mz <- mScan_c[id_m]
    mm <- mz
  } else if (mScan == "m") {
    mx_c <- mScan_c[id_m]
    my <- mScan_c[id_m]
    mz <- mScan_c[id_m]
    mm <- mz
  }
} else if (nohup == F | testRun == "local")
  mm <- mz

### settings paths
if (testRun == "local") {
  setwd("")
  mdl_pth <- ""
  res_pth <- ""
  library(CoRC)
} else if (testRun == "cluster") {
  .libPaths("/net/home.isilon/ag-pahle/Sirac/Rlibrary")
  setwd("/net/home.isilon/ag-pahle/Sirac/Master/Kraskov/")
  mdl_pth <- "/net/home.isilon/ag-pahle/Sirac/Master/data/kholodenko2000/"
  
  if (case_XY == 1 ) {
    caseTE <- "uncond"
  } else
    caseTE <- "cond"
  res_BaseNme <- ""
  res_PthNme <- paste0(res_BaseNme, "N", N, "/", modelType, "/", caseTE, "/")
  dir.create(paste0(res_PthNme, "results"), recursive = T)
  # dir.create(paste0(res_PthNme, "files"), recursive = T)
  dir.create(paste0(res_PthNme, "plots"), recursive = T)
  res_pth <- paste0(res_PthNme, "results/")
  
  library("backports", lib.loc = "/net/home.isilon/ag-pahle/Sirac/Rlibrary")
  library("CoRC", lib.loc = "/net/home.isilon/ag-pahle/Sirac/Rlibrary")
}

### loading libraries

# The Kraskov estimator of TE  KraskovC_TE(x,y,k) uses the k-th neighbour method for estimating of
# Transfer entropy of  x and y
system("R CMD SHLIB KNNnumber.c")
source("KraskovTE.R")
# source("TE_KraskovData.R")
### corrected script
source("TE_KraskovData-mod.R")

### loading model CoRC

if (modelType == "os")
  mdl_nme <- "modelM1.cps"
if (modelType == "ss")
  mdl_nme <- "modelM2.cps"
loadModel(paste0(mdl_pth, mdl_nme))

### parameters for generating data out of model

# method for data genereation
tc_mthd <- "directMethod"
# tc_mthd <- "tauLeap"
# tc_mthd <- "deterministic"

## choose interval size (dt) and duration (dur) of simulated data to test
## or number of data points N to generate in simulation:
dur_c <- N * dt

### test parameteres
# how often to repeat data generation for replicate TE estimation
# rpt <- 20
rpt <- 10
# how often to repeat permutation generation for replicate TE estimation
# rptSur <- 1000
rptSur <- 10
# version of test Run
run <- 1

### level of noise addition and data scaling before TE estimation
nlevel <- 1e-9

### running TE estimation

cat("mz=", mz)
for (mx in mx_c) {
  
  cat("mx=", mx)
  # drr <- 1
  for (drr in 1:length(dur_c)) {
    
    dr <- dur_c[drr]
    cat("dur:", dr, "dt:", dt, "N:", N, "\n")
    
    for (pair in 1:length(pairLst)) {
      
      print(Sys.time())
      
      targetId <- pairLst[[pair]][1]
      driverId <- pairLst[[pair]][2]
      
      res <- matrix(NA, rpt, 9)
      res[, 1] <- 1:rpt
      res[, 2] <- rep(N, rpt)
      res[, 3] <- rep(mx, rpt)
      res[, 4] <- rep(my, rpt)
      res[, 5] <- rep(mz, rpt)
      res[, 9] <- rep("org", rpt)
      
      resSrgt <- matrix(NA, rpt * rptSur, 9)
      resSrgt[, 1] <- rep(1:rpt, rptSur)
      resSrgt[, 2] <- rep(N, rpt * rptSur)
      resSrgt[, 3] <- rep(mx, rpt * rptSur)
      resSrgt[, 4] <- rep(my, rpt * rptSur)
      resSrgt[, 5] <- rep(mz, rpt * rptSur)
      resSrgt[, 9] <- rep("sur", rpt * rptSur)
      
      res_rpt <- c()
      res_sur <- c()
      
      for (rr in 1:rpt) {
        
        cat("rep=", rr)
        
        ### generating data
        tc <- runTimeCourse(duration = dr, dt = dt, 
                            method = tc_mthd, save_result_in_memory = T)
        dat <- tc$result
        
        # adding small amount of noise to data
        len <- dim(dat)[1] * dim(dat)[2]
        noise <- matrix(runif(len, min = -nlevel, max = nlevel), dim(dat)[1])
        dat_ns <- dat + noise
        
        # z-scaling data
        dat_ns <- scale(dat_ns)
        
        # target <- dat_ns[,"Erk2-PP"]    # col 9
        # driver <- dat_ns[,"Mos-P"]      # col 8
        target <- dat_ns[, targetId]
        driver <- dat_ns[, driverId]
        
        # permutation surrogate 
        # driver <- sample(driver)
        
        # making sure TE_KraskovData gets matrix structure
        cond <- as.matrix(dat_ns[, c(-1, -targetId, -driverId)])
        
        ###  TE estimation ###
        res_rr <- TE_KraskovData(target, driver, mx = mx, my = my, 
                                 case_XY = case_XY, mz = mz, envZ = cond, 
                                 k = k_kras, rep = 1)
        
        res_rpt[rr] <- res_rr
        
        ### loop for permutation surrogates
        if (surTest == T) {
          for (rrSur in 1:rptSur) {
            
            driver_sur <- sample(driver)
            
            res_rrSur <- TE_KraskovData(target, driver_sur, mx = mx, my = my, 
                                        case_XY = case_XY, mz = mz, envZ = cond, 
                                        k = k_kras, rep = 1)
            
            res_sur[rrSur + rptSur * (rr-1)] <- res_rrSur
          }
        }
      }
      
      targetNm <- colnames(dat)[targetId]
      driverNm <- colnames(dat)[driverId]
      
      res[, 6] <- rep(targetNm, rpt)
      res[, 7] <- rep(driverNm, rpt)
      res[, 8] <- res_rpt
      
      if (surTest == T) {
        resSrgt[, 6] <- rep(targetNm, rpt * rptSur)
        resSrgt[, 7] <- rep(driverNm, rpt * rptSur)
        resSrgt[, 8] <- res_sur
      }
      
      
      print(Sys.time())
      cat("dur:", dr, "dt:", dt, "N:", N, "\n")
      cat("target:", targetNm, "driver:", driverNm)
      
      ### writing results
      colnames(res) <- c("#", "N", "mx", "my", "mz", "target", "driver", "TE", "type")
      # print(res)
      print(formatC(res, digits = 15))
      
      base_nme <- paste0("simpleModel-biokin_", modelType, "-tc_")
      res_nme <- paste0(base_nme, tc_mthd, "_N", N, "_dt", dt,
                        "-teEst_mx", mx, "_my", my,  "_mz", mz, 
                        "_cond", case_XY, "_k", k_kras, 
                        "-X:", targetNm, "_Y:", driverNm, 
                        "-rep", rpt, "_0", run)
      # write.csv(res, file = paste0(res_pth, res_nme))
      write.table(res, file = paste0(res_pth, paste0(res_nme, "-org.txt")))
      if (surTest == T) {
        colnames(resSrgt) <- c("#", "N", "mx", "my", "mz", "target", "driver", "TE", "type")
        write.table(resSrgt, file = paste0(res_pth, paste0(res_nme, "-srgt.txt")))
      }
    }
  }
}

#####
