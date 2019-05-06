###### TEKraskovData on simpleModel: TE estimation on original and surrogate data ######
### data analysis: for scanning all pairs in the system & mz ###


#### loading all original data ####

rm(list = ls())
library(plyr)
# library(dplyr)

loadSrgts <- T

### loading paths: without "results/" at the end!
dt <- 16
nScan <- c(50, 100, 200, 500, 1000, 2000, 4000, 8000, 16000, 32000)
# mScan <- "mz"
# mScan <- "mx"
mScan <- "my"
# te_cond <- "cond"
te_cond <- "uncond"
modelType <- "os"
# modelType <- "ss"

basePth <- ""


resAllData <- as.data.frame(matrix(NA, 1, 9))
colnames(resAllData) <- c("#", "N", "mx", "my", "mz", "target", "driver", "TE", "type")

for (i in 1:length(nScan)) {
  # dt <- nScan[i]
  N <- nScan[i]
  endPth <- paste0("N", N, "/", modelType, "/", te_cond, "/")
  resPth <- paste0(basePth, endPth, "results/")
  
  # file_list <- list.files(resPth, pattern = "simpleModel")
  if (mScan == "mx") {
    file_list <- list.files(resPth, pattern = "my2_mz2")
  } else if (mScan == "my") {
    file_list <- list.files(resPth, pattern = 'mx2.+mz2')
    # file_list <- list.files(resPth, pattern = 'mx2.*mz2')
  } else if (mScan == "mz") {
    file_list <- list.files(resPth, pattern = "mx2_my2")
  }
  
  resDtPairs <- ldply(paste0(resPth, file_list), function(x) read.table(x, skip = 1))
  resDtPairs <- resDtPairs[, -1]
  
  colnames(resDtPairs) <- c("#", "N", "mx", "my", "mz", "target", "driver", "TE", "type")
  resAllData <- rbind(resAllData, resDtPairs)
}
resAllData <- resAllData[-1, ]

resSumData <- ddply(resAllData, c("target", "driver", "type", mScan, "N"), summarise,
                    nn    = length(TE),
                    mean = mean(TE),
                    sd   = sd(TE),
                    se   = sd / sqrt(nn)
)

#### plotting all or summed data: for one target with surrogates - serial plot saving ####

library(ggplot2)
library(dplyr)
# library(ggpubr)
library(zoo)
library(ggrepel)


### path for saving plots
pltPth <- paste0(basePth, "dtPlots/", modelType, "/", te_cond, "/")
dir.create(pltPth, recursive = T)
# pltPth <- paste0(dat_pth, "plots/")

### parameteres of data structure
# Nsur <- 1000                # number of surrogate TE estimations per original data
Nsur <- 10                # number of surrogate TE estimations per original data

### parameters of data processing
sumData <- T              # plotting summed up data: T, else F
mFilter <- c(1,2,3,4)   # filtering for specific m values
# mFilter <- c(3,4,5,6)   # filtering for specific m values
# mFilter <- c(1,2,3,4,5,6)   # filtering for specific m values
# mFilter <- 1:6            # filtering for specific m values

### parameters of significance test
sigPlot <- F              # plotting significance levels to dt-Scan
sigMethod <- "wTest"      # test type: wilcox-test:= "wTest"; t-test:= "tTest"; sigma Test:= "sigTest"
alpha <- 0.01               # alpha level for Holmer-Bonferroni correction for wilcox-test and t-test
sigma <- 2                  # sigma level as threshold in Holmer-Bonferroni correction for sigma Test
dt_c <- c(8,16,32,64)     # which dt values to consider for significance Test
ver <- 1                  # version of plots for plot naming

### plotting parameters
dodging <- T
# xlab <- "N"
xlab <- expression(N)
# xlab2 <- mScan
ptSize <- 2
wdth <- 0.25


# te_cond <- "cond"
# te_cond <- "uncond"
# modelType <- "os"
# modelType <- "lin001"
if (te_cond == "cond") {
  ylab <- expression(list(TE[X %<-% Y ~ "|" ~ ZY[-1]] ~ " [bit]"))
} else if (te_cond == "uncond") {
  ylab <- expression(TE[X %<-% Y] ~ " [bit]")
}
if (mScan == "mz")
  xlab2 <- expression(m[z])
if (mScan == "mx")
  xlab2 <- expression(m[x])
if (mScan == "my")
  xlab2 <- expression(m[y])

if (sumData == T) {
  resAllPlt <- resSumData
  
  if (sigPlot == T) {
    resTestData <- resAllData
    ### preparing data.frame for plot: factorizing columns
    resTestData$target <- as.factor(resTestData$target)
    resTestData$driver <- as.factor(resTestData$driver)
    
    ### substituting species names
    resTestData <- data.frame(lapply(resTestData, function(x) {
      if (is.factor(x)) {
        gsub("Erk2", "MK", x)
      } else
        x <- x
    }))
    resTestData <- data.frame(lapply(resTestData, function(x) {
      if (is.factor(x)) {
        gsub("Mek1", "M2K", x)
      } else
        x <- x
    }))
    resTestData <- data.frame(lapply(resTestData, function(x) {
      if (is.factor(x)) {
        gsub("Mos", "M3K", x)
      } else
        x <- x
    }))
    resTestData <- data.frame(lapply(resTestData, function(x) {
      if (is.factor(x)) {
        gsub("-PP", "pp", x)
      } else
        x <- x
    }))
    resTestData <- data.frame(lapply(resTestData, function(x) {
      if (is.factor(x)) {
        gsub("-P", "p", x)
      } else
        x <- x
    }))
    
    ### filtering out specific m values
    if (mScan == "mz") {
      resTestData$mz <- as.factor(resTestData$mz)
      resTestFltr <- filter(resTestData, mz %in% mFilter)
    } else if (mScan == "mx") {
      resTestData$mx <- as.factor(resTestData$mx)
      resTestFltr <- filter(resTestData, mx %in% mFilter)
    } else if (mScan == "my") {
      resTestData$my <- as.factor(resTestData$my)
      resTestFltr <- filter(resTestData, my %in% mFilter)
    }
  }
} else if (sumData == F) {
  resAllPlt <- resAllData
}

### preparing data.frame for plot: factorizing columns
resAllPlt$target <- as.factor(resAllPlt$target)
resAllPlt$driver <- as.factor(resAllPlt$driver)

### substituting species names
resAllPlt <- data.frame(lapply(resAllPlt, function(x) {
  if (is.factor(x)) {
    gsub("Erk2", "MK", x)
  } else
    x <- x
}))
resAllPlt <- data.frame(lapply(resAllPlt, function(x) {
  if (is.factor(x)) {
    gsub("Mek1", "M2K", x)
  } else
    x <- x
}))
resAllPlt <- data.frame(lapply(resAllPlt, function(x) {
  if (is.factor(x)) {
    gsub("Mos", "M3K", x)
  } else
    x <- x
}))
resAllPlt <- data.frame(lapply(resAllPlt, function(x) {
  if (is.factor(x)) {
    gsub("-PP", "pp", x)
  } else
    x <- x
}))
resAllPlt <- data.frame(lapply(resAllPlt, function(x) {
  if (is.factor(x)) {
    gsub("-P", "p", x)
  } else
    x <- x
}))

### filtering out specific m values
if (mScan == "mz") {
  resAllPlt$mz <- as.factor(resAllPlt$mz)
  resFltrPlt <- filter(resAllPlt, mz %in% mFilter)
} else if (mScan == "mx") {
  resAllPlt$mx <- as.factor(resAllPlt$mx)
  resFltrPlt <- filter(resAllPlt, mx %in% mFilter)
} else if (mScan == "my") {
  resAllPlt$my <- as.factor(resAllPlt$my)
  resFltrPlt <- filter(resAllPlt, my %in% mFilter)
}



### looping through all targets and drivers for each plot
targetList <- as.vector(unique(resFltrPlt$target))
for (trgt in targetList) {
  # for (i in 1:1) {
  #   trgt <- targetList[i]
  #   trgt <- "A"
  
  cat("\n", "target:", trgt, "\n")
  # print(dt)
  resTrgtPlt <- filter(resFltrPlt, target == trgt)
  
  driverList <- as.vector(unique(resTrgtPlt$driver))
  for (drvr in driverList) {
    # for (j in 1:1) {
    #   drvr <- driverList[j]
    #   drvr <- "Ai"
    
    cat("\n", "driver:", drvr, "\n")
    resDrvrPlt <- filter(resTrgtPlt, driver == drvr)

    resDrvrPlt$driver <- factor(resDrvrPlt$driver,
                                levels = c("MK", "MKp", "MKpp",
                                           "M2K", "M2Kp", "M2Kpp",
                                           "M3K", "M3Kpp"))
    
    if (te_cond == "cond") {
      ylab <- bquote(TE[.(drvr) %->% .(trgt) ~ "|" ~ Z ~ .(drvr)[-1]] ~ " [bit]")
    } else if (te_cond == "uncond") {
      ylab <- bquote(TE[.(drvr) %->% .(trgt)] ~ " [bit]")
    }
    
    
    ### calculating for significance test p-value matrix
    if (sigPlot == T) {
      m_c <- unique(resDrvrPlt[, 4])
      dt_full <- unique(resDrvrPlt$d_t)
      len_m <- length(m_c)
      len_dt <- length(dt_c)
      
      sigTest <- matrix(list(), len_dt, len_m)
      
      for (dtt in 1:length(dt_c)) {
        # for (dtt in 1:1) {
        for (mm in 1:length(m_c)) {
          # for (mzz in 1:1) {
          cat(dt_c[dtt], m_c[mm], "\t")
          
          ### preparing data.frame for significance test
          if (sumData == T) {
            resTest <- filter(resTestFltr, target == trgt & driver == drvr)
            resTest <- filter(resTest, d_t == dt_c[dtt] & !!ensym(mScan) == m_c[mm])
          } else if (sumData == F) {
            resTest <- filter(resDrvrPlt, d_t == dt_c[dtt] & !!ensym(mScan) == m_c[mm])
          }
          if (nrow(resTest) == 0)
            break
          
          resTest_org <- filter(resTest, type == "org")
          resTest_sur <- filter(resTest, type == "sur")
          
          ### significance test between original and surrogate estimates
          if (sigMethod == "tTest") {

            pList <- c()
            for (i in 1:nrow(resTest_org)) {
              surSeq <- 1:Nsur + Nsur*(i-1)
              shapTest <- shapiro.test(resTest_sur$TE[surSeq])
              if (shapTest$p.value > 0.05) {
                cat("test for normality passed", shapTest$p.value, "\n")
              } else if (shapTest$p.value <= 0.05) {
                cat("test for normality failed", shapTest$p.value, "\n")
              }
              # tTest <- t.test(resTest_org$TE, resTest_sur$TE, paired = F, alternative = "less")
              tTest <- t.test(resTest_sur$TE[surSeq], mu = resTest_org$TE[i], alternative = "less", conf.level = 0.99)
              pList[i] <- tTest$p.value
            }
            sigTest[dtt, mm][[1]] <- pList
            # browser()
            
          } else if (sigMethod == "wTest") {
            
            pList <- c()
            for (i in 1:nrow(resTest_org)) {
              surSeq <- 1:Nsur + Nsur*(i-1)
              wTest <- wilcox.test(resTest_sur$TE[surSeq], mu = resTest_org$TE[i], alternative = "less", conf.level = 0.99)
              # tTest <- t.test(resTest_org$TE[i], resTest_sur$TE[surSeq], paired = F)
              pList[i] <- wTest$p.value
            }
            # sigTest[dtt, mzz][1] <- as.list(pList)
            sigTest[dtt, mm][[1]] <- pList
            
          } else if (sigMethod == "sigTest") {
            
            pList <- c()
            for (i in 1:nrow(resTest_org)) {
              surSeq <- 1:Nsur + Nsur*(i-1)
              # sTest <- (resTest_org$TE[i] - mean(resTest_sur$TE[surSeq])) / sd(resTest_sur$TE[surSeq])
              sTest <- (resTest_org$TE[i] - mean(resTest_sur$TE[surSeq])) / 
                sd(resTest_sur$TE[surSeq]) / sqrt(length(resTest_sur$TE[surSeq]))
              # tTest <- t.test(resTest_org$TE[i], resTest_sur$TE[surSeq], paired = F)
              pList[i] <- sTest
            }
            sigTest[dtt, mm][[1]] <- pList
          }
          
        }
        if (nrow(resTest) == 0)
          break
      }
      if (nrow(resTest) == 0)
        break
      
      ### Holm-Bonferroni correction for tTest and wTest
      if (sigMethod == "tTest" | sigMethod == "wTest") {
        
        corrSigTest <- matrix(NA, len_dt, len_m)
        rownames(corrSigTest) <- dt_c
        # alpha <- 0.05
        n <- nrow(resTest_org)
        
        for (dtt in 1:length(dt_c)) {
          for (mm in 1:length(m_c)) {
            
            pList <- sigTest[dtt, mm][[1]]
            pList <- sort(pList)
            
            for (k in 1:length(pList)) {
              
              corrAlpha <- alpha / (n-k+1)
              # corrAlpha <- alpha
              
              if (pList[k] > corrAlpha) {
                corrSigTest[dtt, mm] <- ""
                # corrSigTest[dtt, mm] <- paste(dt_c[dtt],m_c[mm])
                break
              } else if (pList[k] <= corrAlpha & k == n) {
                corrSigTest[dtt, mm] <- "*"
                # corrSigTest[dtt, mm] <- paste(dt_c[dtt],m_c[mm])
              }
            }
          }
        }
        
        ### Holmer-Bonferroni correction for sigma Test
      } else if (sigMethod == "sigTest") {
        
        corrSigTest <- matrix(NA, len_dt, len_m)
        rownames(corrSigTest) <- dt_c
        # sigma <- 6
        n <- nrow(resTest_org)
        
        for (dtt in 1:length(dt_c)) {
          for (mm in 1:length(m_c)) {
            
            pList <- sigTest[dtt, mm][[1]]
            pList <- sort(pList)
            
            for (k in 1:length(pList)) {
              
              # corrPval <- pList[k] / (n-k+1)
              corrPval <- pList[k]
              
              if (corrPval <= sigma) {
                corrSigTest[dtt, mm] <- ""
                # corrSigTest[dtt, mm] <- paste(dt_c[dtt],m_c[mm])
                break
              } else if (corrPval > sigma & k == n) {
                corrSigTest[dtt, mm] <- "*"
                # corrSigTest[dtt, mm] <- paste(dt_c[dtt],m_c[mm])
              }
            }
          }
        }
      }
      
      ### completing corrSigTest matrix for subsequent cbind
      corrSigTest_full <- matrix("", length(dt_full), len_m)
      for (dtff in 1:length(dt_full)) {
        for (dtcc in 1:length(dt_c)) {
          if (as.character(dt_full[dtff]) == rownames(corrSigTest)[dtcc]) {
            corrSigTest_full[dtff, ] <- corrSigTest[dtcc, ]
          }
        }
      }
      
      ### adding significance level to data.frame
      sigLevel <- as.vector(corrSigTest_full)
      if (sumData == T) {
        sigLevel <- c(sigLevel, rep(NA, length(sigLevel)))
        resDrvrPlt <- cbind(resDrvrPlt, sigLevel)
        # colnames(resDrvrPlt) <- c(colnames(resDrvrPlt), "sig")
      } else if (sumData == F) {
        # TODO: construct sigLevel for unsummed Data
      }
      
    }
    
    
    ### plotting script
    
    ## determining ylim
    if (sumData == T) {
      minData <- resDrvrPlt$mean - resDrvrPlt$sd
      maxData <- resDrvrPlt$mean + resDrvrPlt$sd
      if (sigPlot == T) {
        maxData <- maxData + 0.010
      }
      if (min(minData >= 0)) {
        ylim <- c(0, max(maxData))
      } else
        ylim <- c(min(minData), max(maxData))
    } else if (sumData == F) {
      if (min(resDrvrPlt$TE >= 0)) {
        ylim <- c(0, max(resDrvrPlt$TE))
      } else
        ylim <- c(min(resDrvrPlt$TE), max(resDrvrPlt$TE))
    }
    
    ## determine type of plot
    if (sumData == T) {
      
      # pl <- ggplot(resDrvrPlt, aes(d_t, mean, group = type, label = type))
      pl <- ggplot(resDrvrPlt, aes(N, mean))
      
      if (dodging == T) {
        # dodge <- position_jitterdodge(jitter.width = wdth, dodge.width = 0.001, seed = 1)
        dodge <- position_dodge(0.025)
        pl <- pl + 
          geom_line(aes(y = mean, col = !!ensym(mScan), linetype = type)) +
          geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd, col = !!ensym(mScan),
                            linetype = type), position = dodge, width = 0.4)
      } else {
        pl <- pl + 
          # geom_point(aes(col = m, shape = type), width = wdth, size = ptSize)
          geom_line(aes(y = mean, col = !!ensym(mScan), linetype = type)) +
          geom_linerange(aes(ymin = mean - sd, ymax = mean + sd, col = !!ensym(mScan), 
                             linetype = type))
      }
      if (sigPlot == T) {
        pl <- pl +
          geom_text_repel(aes(x = (N + N)/1.7, y = mean + 0.0, 
                              label = sigLevel, col = !!ensym(mScan)),
                          nudge_y = 0.01,
                          show.legend = F,
                          direction = "y",
                          segment.alpha = 0)
      }
      pl <- pl +
        scale_shape_manual(values = c(1,4))
      
      
    } else if (sumData == F) {
      
      pl <- ggplot(resDrvrPlt, aes(N, TE))
      
      if (dodging == T) {
        pl <- pl + 
          geom_point(aes(shape = type, col = !!ensym(mScan)), size = ptSize,
                     position = position_jitterdodge(jitter.width = wdth, dodge.width = 0.25))
      } else {
        pl <- pl + 
          geom_point(aes(col = !!ensym(mScan), shape = type), width = wdth, size = ptSize)
      }
      pl <- pl +
        scale_shape_manual(values = c(1,4))
    }
    
    ## color gradient
    pl <- pl +
      scale_color_brewer(palette = "RdYlBu", direction = -1, name = xlab2)
    
    ## theme, legends, margins and text sizes
    pl <- pl + theme_bw() +
      coord_cartesian(ylim = ylim) +
      guides(color = guide_legend(xlab2), shape = "none", linetype = "none") +
      scale_x_log10(breaks = nScan[seq(1, length(nScan), 2)]) +
      labs(x = xlab, y = ylab) +
      theme(text = element_text(size = 20),
            axis.text = element_text(size = 15),
            axis.text.x = element_text(angle = 0, hjust = 0.5),
            axis.title.x = element_text(margin = margin(t = 20)),
            axis.title.y = element_text(margin = margin(r = 10)))
    pl
    
    ### writing pdf
    plt_nme <- paste0("simpleModel_", modelType, 
                      "-teKrasData-dt", dt, "NAll", 
                      "-", trgt, "<-", drvr, "-", 
                      te_cond, "_", mScan, "-", ver)
    ggsave(paste0(plt_nme, ".pdf"), width = 8.5, height = 6.5, path = pltPth, plot = pl)
  }
  if (exists("resTest")) {
    if (nrow(resTest) == 0)
      break
  }
}
