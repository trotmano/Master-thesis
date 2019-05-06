#########################################
TE_KraskovData <- function(targetX, driverY, mx, my, case_XY = 1,
                           mz = 0, envZ = NULL, 
                           k = 3, rep = 1){
  
  
  #####################################################################
  #  Kraskov estimation param k
  k_Kraskov <- k
  #rep = 10
  
  # ???
  TE_Km <- c()
  
  if (!(is.vector(targetX)) | !(is.vector(driverY))) {
    stop("xData or yData is not a vector!")
  }
  if (length(targetX) != length(driverY)) {
    stop("data length of X and Y differ!")
  }
  N <- length(targetX)
  mk <- max(mx,my,mz)
  
  for(m in 1:rep){
    
    #####################################################################
    
    #  For calculating TE_{X<-Y}
    if(case_XY == 1){
      
      if (mk > (N-1)) {
        stop("mx/my too big for data length")
      }
      
      # Target X
      X <- matrix(ncol = mx, nrow = N-mk+1)
      for(l in 1:mx)
        X[,l] <- targetX[(mk+1-l):(N+1-l)]
    }
    
    
    #  For calculating TE_{X<-Y|Z}
    if(case_XY == 2){
      
      if (mk > (N-1)) {
        stop("mx/my too big for data length")
      }
      
      if (is.null(envZ)) {
        stop("no environment variable Z given!") 
      } else if (is.vector(envZ)) {
        nCon <- 1
        if (length(envZ) != N) {
          stop("data length of targetX and envZ differ!")
        }
      } else if (is.numeric(ncol(envZ))) {
        nCon <- ncol(envZ)
        if (nrow(envZ) != N) {
          stop("data length of targetX and envZ differ!")
        }
      } else {
        stop("check data structure of variabel Z!")
      }
      nncol <- mx + mz * nCon
      
      # browser()
      
      X <- matrix(ncol = nncol, nrow = N-mk+1)
      # X <- matrix(ncol = mx+mz, nrow=N-mk+1)
      # browser()
      
      # Target X
      for(l in 1:mx)
        X[,l] <- targetX[(mk+1-l):(N+1-l)]
      
      
      # Conditionining environment Z
      if (is.vector(envZ)) {
        for(l in 1:mz){ 
          lz <- mx + l
          # browser()
          X[,lz] <- envZ[(mk+1-l):(N+1-l)]
        }
      } else {
        for (cc in 1:nCon) {
          for(l in 1:mz){ 
            lz <- mx + l + mz*(cc-1)
            # browser()
            X[,lz] <- envZ[(mk+1-l):(N+1-l), cc]
          }
        }
      }
      # browser()
    }
    
    
    #  For calculating TE_{X<-Y|ZY-}
    if(case_XY == 3){
      
      if (mk >= (N-1)) {
        stop("mx/my/mz too big for data length")
      }
      
      mk = max(mx, mz, my+1)
      mky = max(mx-1, mz-1, my)
      
      if (mk > N) {
        stop("mx/my/mz too big for data length")
      }
      
      if (is.null(envZ)) {
        stop("no environment variable Z given!") 
      } else if (is.vector(envZ)) {
        nCon <- 1
        if (length(envZ) != N) {
          stop("data length of targetX and envZ differ!")
        }
      } else if (is.numeric(ncol(envZ))) {
        nCon <- ncol(envZ)
        if (nrow(envZ) != N) {
          stop("data length of targetX and envZ differ!")
        }
      } else {
        stop("check data structure of variabel Z!")
      }
      nncol <- mx + mz * nCon + my
      
      X <- matrix(ncol = nncol, nrow = N-mk+1)
      # X <- matrix(ncol = mx+mz, nrow=N-mk+1)
      
      # Target X
      for(l in 1:mx)
        X[,l] <- targetX[(mk+1-l):(N+1-l)]
      
      
      # Conditionining environment Z
      if (is.vector(envZ)) {
        for(l in 1:mz){ 
          lz <- mx + l
          X[,lz] <- envZ[(mk+1-l):(N+1-l)]
        }
      } else {
        for (cc in 1:nCon) {
          for(l in 1:mz){ 
            lz <- mx + l + mz*(cc-1)
            # browser()
            X[,lz] <- envZ[(mk+1-l):(N+1-l), cc]
          }
        }
      }
      
      # Conditionining variable driver Y
      # TODO: Check Here
      for(l in 1:my){
        # TODO: How to choose index lz
        ly <- mx + mz * nCon + l
        # browser()
        X[,ly] <- driverY[(mky+1-l):(N-l)]
      }
    }
    
    # driver Y
    Y <- matrix(ncol = my, nrow=N-mk+1)
    for(l in 1:my)
      Y[,l] <- driverY[(mk+1-l):(N+1-l)]
    
    
    TE_Km[m] <- KraskovTE(X,Y,k_Kraskov)
    # browser()
  }
  
  TEXY_Kr <- mean(TE_Km)
  res <- TEXY_Kr
  # browser()
  return(res)
  # return(X)
}
