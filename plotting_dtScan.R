###### TEKraskovData on simpleModel: TE estimation on original and surrogate data ######
### data analysis: for scanning all pairs in the system & mz ###


#### loading all original data ####

rm(list = ls())
library(plyr)

loadSrgts <- T

### loading paths: without "results/" at the end!
N <- 20000
dtScan <- c(64, 32, 16, 8, 4, 2, 1, 0.5, 0.1, 0.05, 0.01)
# modelType <- "os"
modelType <- "ss"
# te_cond <- "cond"
te_cond <- "uncond"
# mScan <- "mz"
# mScan <- "mx"
mScan <- "my"

basePth <- ""

resAllData <- as.data.frame(matrix(NA, 1, 10))
colnames(resAllData) <- c("#", "N", "mx", "my", "mz", "target", "driver", "TE", "type", "d_t")

for (i in 1:length(dtScan)) {
  dt <- dtScan[i]
  endPth <- paste0("dt", dt, "/", modelType, "/", te_cond, "/")
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
  resDtPairs <- cbind(resDtPairs, rep(dt, dim(resDtPairs)[1]))
  
  # correction for different column lengths
  if (dim(resDtPairs)[2] != 10) {
    resDtPairs <- cbind(resDtPairs, rep(2, dim(resDtPairs)[1]), rep(2, dim(resDtPairs)[1]))
    if (mScan == "mz") {
      resDtPairs <- resDtPairs[, c(1, 2, 9, 10, 3, 4, 5, 6, 7, 8)]
    } else if (mScan == "mx") {
      resDtPairs <- resDtPairs[, c(1, 2, 3, 9, 10, 4, 5, 6, 7, 8)]
    } else if (mScan == "my") {
      resDtPairs <- resDtPairs[, c(1, 2, 9, 3, 10, 4, 5, 6, 7, 8)]
    }
  }
  
  colnames(resDtPairs) <- c("#", "N", "mx", "my", "mz", "target", "driver", "TE", "type", "d_t")
  resAllData <- rbind(resAllData, resDtPairs)
}
resAllData <- resAllData[-1, ]

resSumData <- ddply(resAllData, c("target", "driver", "type", mScan, "d_t"), summarise,
                    N    = length(TE),
                    mean = mean(TE),
                    sd   = sd(TE),
                    se   = sd / sqrt(N)
)

#### plotting all or summed data: for one target with surrogates - serial plot saving ####

library(ggplot2)
library(dplyr)


### path for saving plots
pltPth <- paste0(basePth, "dtPlots/", modelType, "/", te_cond, "/")
dir.create(pltPth, recursive = T)

### plot version
sumData <- T              # plotting summed up data: T, else F
mFilter <- c(1,2,3)   # filtering for specific m values
ver <- 1                  # version of plots for plot naming

### plotting parameters
dodging <- T
# xlab <- "N"
xlab <- expression(delta ~ t ~ " [a.u.]")
ptSize <- 2
wdth <- 0.25


if (te_cond == "cond") {
  # ylab <- expression(list(TE[X %<-% Y ~ "|" ~ Z]))
  ylab <- expression(list(TE[X %<-% Y ~ "|" ~ ZY[-1]] ~ " [bit]"))
  # ylab <- expression(list(TE[X %<-% Y ~ "|" ~ ZY ~ "⁻"])),
  # ylab <- "TE_X<-Y|XY⁻"
} else if (te_cond == "uncond") {
  ylab <- expression(TE[X %<-% Y] ~ " [bit]")
}
if (mScan == "mz")
  xlab2 <- expression(m[z])
if (mScan == "mx")
  xlab2 <- expression(m[x])
if (mScan == "my")
  xlab2 <- expression(m[y])
# mScan <- expression(m[x/y/z])
# mScan <- expression(m[x] / m[y] / m[z])

if (sumData == T) {
  resAllPlt <- resSumData
} else
  resAllPlt <- resAllData

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
  
  cat("target:", trgt, "\n")
  # print(dt)
  resTrgtPlt <- filter(resFltrPlt, target == trgt)
  
  driverList <- as.vector(unique(resTrgtPlt$driver))
  for (drvr in driverList) {
    
    drvrNme <- drvr
    trgtNme <- trgt
    
    cat("driver:", drvr, "\n")
    resDrvrPlt <- filter(resTrgtPlt, driver == drvr)
    
    ### factorizing columns and ordering factors
    resDrvrPlt$driver <- factor(resDrvrPlt$driver,
                                levels = c("MK", "MKp", "MKpp",
                                           "M2K", "M2Kp", "M2Kpp",
                                           "M3K", "M3Kpp"))
    
    if (te_cond == "cond") {
      ylab <- bquote(TE[.(drvr) %->% .(trgt) ~ "|" ~ Z ~ .(drvr)[-1]] ~ " [bit]")
    } else if (te_cond == "uncond") {
      ylab <- bquote(TE[.(drvr) %->% .(trgt)] ~ " [bit]")
    }
    
    ### plotting script
    
    ## determining ylim
    if (sumData == T) {
      minData <- resDrvrPlt$mean - resDrvrPlt$sd
      maxData <- resDrvrPlt$mean + resDrvrPlt$sd
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
      pl <- ggplot(resDrvrPlt, aes(d_t, mean))
      if (dodging == T) {
        dodge <- position_dodge(0.025)
        pl <- pl + 
          geom_line(aes(y = mean, col = !!ensym(mScan), linetype = type)) +
          geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd, col = !!ensym(mScan),
                            linetype = type), position = dodge, width = 0.4)
        
      } else {
        pl <- pl + 
          geom_line(aes(y = mean, col = !!ensym(mScan), linetype = type)) +
          geom_linerange(aes(ymin = mean - sd, ymax = mean + sd, col = !!ensym(mScan), 
                             linetype = type))
      }
      pl <- pl +
        scale_shape_manual(values = c(1,4))
      
    } else if (sumData == F) {
      pl <- ggplot(resDrvrPlt, aes(d_t, TE))
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
      scale_color_manual(values = c("dodgerblue3", "green4", "brown3"))
      # scale_color_brewer(palette = "RdYlBu", direction = -1, name = xlab2)

    
    ## theme, legends, margins and text sizes
    pl <- pl + theme_bw() +
      coord_cartesian(ylim = ylim) +
      guides(color = guide_legend(xlab2), shape = "none", linetype = "none") +
      scale_x_log10(breaks = dtScan[seq(1, length(dtScan), 2)]) +
      labs(x = xlab, y = ylab) +
      theme(text = element_text(size = 20),
            axis.text = element_text(size = 15),
            axis.text.x = element_text(angle = 0, hjust = 0.5),
            axis.title.x = element_text(margin = margin(t = 20)),
            axis.title.y = element_text(margin = margin(r = 10)))
    pl
    
    ### writing pdf
    plt_nme <- paste0("kholo002_", modelType, 
                      "-teKrasData-N", N, "dtAll", 
                      "-", trgt, "<-", drvr, "-", 
                      te_cond, "_", mScan, "-", ver)
    ggsave(paste0(plt_nme, ".pdf"), width = 8.5, height = 6.5, path = pltPth, plot = pl)
  }
}
# pl

#####
