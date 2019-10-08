# dat is a matrix of sample X metabolite



#### Mean Centering ####

MeanCenter <- function(dat, batch) { # dat is a matrix of sample X metabolite
  u <- aggregate(dat, by = list(batch), mean) # get mean value for each metabolite, grouped by batch
  u <- u[batch, -1]
  uw <- apply(dat, 2, mean)
  uw <- t(replicate(length(batch), uw))
  dif <- u - uw
  dat <- dat - dif
  dat <- as.data.frame(t(dat))
  return(dat)
}


#### Median Centerng ####

MedianCenter <- function(dat, batch) { # dat is a matrix of sample X metabolite
  u <- aggregate(dat, by = list(batch), median) # get median value for each metabolite, grouped by batch
  u <- u[batch, -1]
  uw <- apply(dat, 2, mean)
  uw <- t(replicate(length(batch), uw))
  dif <- u - uw
  dat <- dat - dif
  dat <- as.data.frame(t(dat))
  return(dat)
}


#### Quantile ####

Quantile <- function(dat) { # dat is a matrix of metabolite X sample
  mzid <- colnames(dat)
  cn <- rownames(dat)
  dat <- as.data.frame(t(dat))
  oo <- apply(dat, 2, function(x) rank(x, ties.method = "first"))
  dat <- apply(dat, 2, sort)
  u <- apply(dat, 1, mean)
  dat <- apply(dat, 2, function(x) x=u)
  dat <- sapply(1:ncol(dat), function(i) dat[oo[, i], i])
  rownames(dat) <- mzid
  colnames(dat) <- cn
  dat <- as.data.frame(dat)
  return(dat)
}


#### Combat ####
# This is a modification of the ComBat function code from the sva package that can be found at
# https://bioconductor.org/packages/release/bioc/html/sva.html 

#' Adjust for batch effects using an empirical Bayes framework
#'
#' ComBat allows users to adjust for batch effects in datasets where the batch covariate is known, using methodology
#' described in Johnson et al. 2007. It uses either parametric or non-parametric empirical Bayes frameworks for adjusting data for
#' batch effects.  Users are returned an expression matrix that has been corrected for batch effects. The input
#' data are assumed to be cleaned and normalized before batch effect removal.
#'
#' @param dat Genomic measure matrix (dimensions probe x sample) - for example, expression matrix
#' @param batch {Batch covariate (only one batch allowed)}
#' @param mod Model matrix for outcome of interest and other covariates besides batch
#' @param par.prior (Optional) TRUE indicates parametric adjustments will be used, FALSE indicates non-parametric adjustments will be used
#' @param prior.plots (Optional) TRUE give prior plots with black as a kernel estimate of the empirical batch effect density and red as the parametric
#' @param mean.only (Optional) FALSE If TRUE ComBat only corrects the mean of the batch effect (no scale adjustment)
#' @param ref.batch (Optional) NULL If given, will use the selected batch as a reference for batch adjustment.
#' @param BPPARAM (Optional) BiocParallelParam for parallel operation
#'
#' @return data A probe x sample genomic measure matrix, adjusted for batch effects.
#'
#' @importFrom graphics lines par
#' @importFrom stats cor density dnorm model.matrix pf ppoints prcomp predict
#' qgamma qnorm qqline qqnorm qqplot smooth.spline var
#' @importFrom utils read.delim
#'
#'

combat <- function(dat, batch, mod=NULL, eb=TRUE, verbose=TRUE){
  aprior <- function(gamma.hat){m=mean(gamma.hat); s2=var(gamma.hat); (2*s2+m^2)/s2}
  bprior <- function(gamma.hat){m=mean(gamma.hat); s2=var(gamma.hat); (m*s2+m^3)/s2}
  postmean <- function(g.hat,g.bar,n,d.star,t2){(t2*n*g.hat+d.star*g.bar)/(t2*n+d.star)}
  postvar <- function(sum2,n,a,b){(.5*sum2+b)/(n/2+a-1)}
  
  it.sol  <- function(sdat,g.hat,d.hat,g.bar,t2,a,b,conv=.0001){
    n <- apply(!is.na(sdat),1,sum)
    g.old <- g.hat
    d.old <- d.hat
    change <- 1
    count <- 0
    while(change>conv){
      g.new <- postmean(g.hat,g.bar,n,d.old,t2)
      sum2 <- apply((sdat-g.new%*%t(rep(1,ncol(sdat))))^2, 1, sum,na.rm=TRUE)
      d.new <- postvar(sum2,n,a,b)
      change <- max(abs(g.new-g.old)/g.old,abs(d.new-d.old)/d.old)
      g.old <- g.new
      d.old <- d.new
      count <- count+1
    }
    #cat("This batch took", count, "iterations until convergence\n")
    adjust <- rbind(g.new, d.new)
    rownames(adjust) <- c("g.star","d.star")
    adjust
  }
  
  dat <- as.data.frame(t(dat))
  
  .checkConstantRows <- function(dat){
    sds <- unlist(apply(dat, 1, sd))
    ns <- sum(sds==0)
    if (ns>0){
      message <- paste0(ns, " rows (features) were found to be constant across samples. Please remove these rows before running ComBat.")
      stop(message)
    }
  }
  .checkConstantRows(dat)
  if (eb){
    if (verbose) cat("[combat] Performing ComBat with empirical Bayes\n")
  } else {
    if (verbose) cat("[combat] Performing ComBat without empirical Bayes (L/S model)\n")
  }
  # make batch a factor and make a set of indicators for batch
  batch <- as.factor(batch)
  batchmod <- model.matrix(~-1+batch)  
  if (verbose) cat("[combat] Found",nlevels(batch),'batches\n')
  
  # A few other characteristics on the batches
  n.batch <- nlevels(batch)
  batches <- list()
  for (i in 1:n.batch){
    batches[[i]] <- which(batch == levels(batch)[i])
  } # list of samples in each batch  
  n.batches <- sapply(batches, length)
  n.array <- sum(n.batches)
  
  #combine batch variable and covariates
  design <- cbind(batchmod,mod)
  # check for intercept in covariates, and drop if present
  check <- apply(design, 2, function(x) all(x == 1))
  design <- as.matrix(design[,!check])
  
  # Number of covariates or covariate levels
  if (verbose) cat("[combat] Adjusting for",ncol(design)-ncol(batchmod),'covariate(s) or covariate level(s)\n')
  
  # Check if the design is confounded
  if(qr(design)$rank<ncol(design)){
    if(ncol(design)==(n.batch+1)){
      stop("[combat] The covariate is confounded with batch. Remove the covariate and rerun ComBat.")
    }
    if(ncol(design)>(n.batch+1)){
      if((qr(design[,-c(1:n.batch)])$rank<ncol(design[,-c(1:n.batch)]))){
        stop('The covariates are confounded. Please remove one or more of the covariates so the design is not confounded.')
      } else {
        stop("At least one covariate is confounded with batch. Please remove confounded covariates and rerun ComBat.")
      }
    }
  }
  
  
  ##Standardize Data across features
  if (verbose) cat('[combat] Standardizing Data across features\n')
  B.hat <- solve(t(design)%*%design)%*%t(design)%*%t(as.matrix(dat))
  
  
  #Standarization Model
  grand.mean <- t(n.batches/n.array)%*%B.hat[1:n.batch,]
  var.pooled <- as.matrix(((dat-t(design%*%B.hat))^2))%*%rep(1/n.array,n.array)
  stand.mean <- t(grand.mean)%*%t(rep(1,n.array))
  if(!is.null(design)){
    tmp <- design;tmp[,c(1:n.batch)] <- 0
    stand.mean <- stand.mean+t(tmp%*%B.hat)
  }	
  s.data <- (dat-stand.mean)/(sqrt(var.pooled)%*%t(rep(1,n.array)))
  
  ##Get regression batch effect parameters
  if (eb){
    if (verbose) cat("[combat] Fitting L/S model and finding priors\n")
  } else {
    if (verbose) cat("[combat] Fitting L/S model\n")
  }
  
  batch.design <- design[,1:n.batch]
  gamma.hat <- solve(t(batch.design)%*%batch.design)%*%t(batch.design)%*%t(as.matrix(s.data))
  
  delta.hat <- NULL
  for (i in batches){
    delta.hat <- rbind(delta.hat,apply(s.data[,i], 1, var,na.rm=T))
  }
  
  # Empirical Bayes correction:
  gamma.star <- delta.star <- NULL
  gamma.bar <- t2 <- a.prior <- b.prior <- NULL
  if (eb){
    ##Find Priors
    gamma.bar <- apply(gamma.hat, 1, mean)
    t2 <- apply(gamma.hat, 1, var)
    a.prior <- apply(delta.hat, 1, aprior)
    b.prior <- apply(delta.hat, 1, bprior)
    
    
    ##Find EB batch adjustments
    if (verbose) cat("[combat] Finding parametric adjustments\n")
    for (i in 1:n.batch){
      temp <- it.sol(s.data[,batches[[i]]],gamma.hat[i,],delta.hat[i,],gamma.bar[i],t2[i],a.prior[i],b.prior[i])
      gamma.star <- rbind(gamma.star,temp[1,])
      delta.star <- rbind(delta.star,temp[2,])
    }
  } 
  
  ### Normalize the Data ###
  if (verbose) cat("[combat] Adjusting the Data\n")
  bayesdata <- s.data
  j <- 1
  for (i in batches){
    if (eb){
      bayesdata[,i] <- (bayesdata[,i]-t(batch.design[i,]%*%gamma.star))/(sqrt(delta.star[j,])%*%t(rep(1,n.batches[j])))
    } else {
      bayesdata[,i] <- (bayesdata[,i]-t(batch.design[i,]%*%gamma.hat))/(sqrt(delta.hat[j,])%*%t(rep(1,n.batches[j])))
    }
    j <- j+1
  }
  
  bayesdata <- (bayesdata*(sqrt(var.pooled)%*%t(rep(1,n.array))))+stand.mean
  #return(list(dat.combat=bayesdata, 
  #           gamma.hat=gamma.hat, delta.hat=delta.hat, 
  #            gamma.star=gamma.star, delta.star=delta.star, 
  #           gamma.bar=gamma.bar, t2=t2, a.prior=a.prior, b.prior=b.prior, batch=batch, mod=mod, 
  #            stand.mean=stand.mean, stand.sd=sqrt(var.pooled)[,1])
  #)
  return(bayesdata)
}



#### RUV ####
# Y is a matrix of sample X metabolite

design.matrix <- function (a, name = "X", remove.collinear = TRUE, include.intercept = TRUE) {
    if (is.vector(a)) 
      a = matrix(a)
    if (is.matrix(a)) {
      if (is.null(colnames(a))) 
        colnames(a) = paste(name, 1:ncol(a), sep = "")
      for (i in 1:ncol(a)) {
        if (is.na(colnames(a)[i])) 
          colnames(a)[i] = paste(name, i, sep = "")
        if (colnames(a)[i] == "") 
          colnames(a)[i] = paste(name, i, sep = "")
      }
    }
    if (is.numeric(a)) 
      A = a
    else {
      if (is.factor(a)) {
        a = data.frame(a)
        names(a)[1] = paste(name, 1, sep = "")
      }
      a = data.frame(a)
      varnames = colnames(a)
      for (i in 1:length(varnames)) if (varnames[i] == paste("X", 
                                                             i, sep = "")) 
        varnames[i] = paste(name, i, sep = "")
      A = design.column(a[, 1], name = varnames[1])
      if (ncol(a) > 1) 
        for (i in 2:ncol(a)) A = cbind(A, design.column(a[, 
                                                          i], name = varnames[i]))
    }
    if (remove.collinear) 
      if (ncol(A) > 1) {
        if (ncol(A) > nrow(A)) 
          A = A[, 1:nrow(A)]
        A0 = scale(A, center = FALSE, scale = TRUE)
        d = svd(A)$d
        if (d[1]/d[length(d)] > 10^9) {
          warning("Collinearity detected.  Removing some columns.")
          toremove = NULL
          for (i in 2:ncol(A0)) {
            A1 = A0[, 1:(i - 1)]
            if (!is.null(toremove)) 
              A1 = A1[, -toremove]
            if (mean(residop(A0[, i, drop = FALSE], A1)^2) < 
                10^(-9)) 
              toremove = c(toremove, i)
          }
          A = A[, -toremove, drop = FALSE]
        }
      }
    if (include.intercept) {
      if (sum(residop(matrix(1, nrow(A), 1), A)^2) > 10^(-8)) {
        A = cbind(1, A)
        colnames(A)[1] = paste(name, "0", sep = "")
      }
    }
    return(A)
  }

getK <- function (Y, X, ctl, Z = 1, eta = NULL, include.intercept = TRUE, 
            fullW0 = NULL, cutoff = NULL, method = "select", l = 1, inputcheck = TRUE) 
  {
    if (is.data.frame(Y)) 
      Y = data.matrix(Y)
    m = nrow(Y)
    n = ncol(Y)
    X = design.matrix(X, include.intercept = FALSE)
    p = ncol(X)
    if (is.numeric(Z)) 
      if (length(Z) == 1) 
        Z = matrix(1, m, 1)
    if (!is.null(Z)) {
      Z = design.matrix(Z, name = "Z", include.intercept = include.intercept)
      q = ncol(Z)
    }
    else q = 0
    ctl = tological(ctl, n)
    nc = sum(ctl)
    if (inputcheck) 
      inputcheck1(Y, X, Z, ctl)
    if (p > 1) 
      return(getK(Y, X[, l, drop = FALSE], ctl, Z = cbind(X[, 
                                                            -l, drop = FALSE], Z), eta = eta, include.intercept = include.intercept, 
                  fullW0 = fullW0, cutoff = cutoff, method = method, 
                  inputcheck = FALSE))
    Y = RUV1(Y, eta, ctl, include.intercept = include.intercept)
    if (q > 0) {
      Y = residop(Y, Z)
      X = residop(X, Z)
    }
    if (method == "select") {
      if (nc <= 1000 & m <= 300) 
        method = "leave1out"
      else method = "fast"
      if (method == "fast" & (nc < 3 * m + 100 | m > 500)) 
        warning("Either m or n_c (or both) are too large.  This algorithm (\"getK\") will be slow, provide poor results, or both.  It is recommended to choose K in some other manner.  To override this warning, set \"method\" to either \"fast\" or \"leave1out.\"  Currently defaulting to method \"fast.\"", 
                immediate. = TRUE)
    }
    K1 = min(m - p - q - 1, nc - 2)
    X = X/sqrt(sum(X^2))
    Yc = Y[, ctl]
    Y0 = residop(Y, X)
    if (is.null(fullW0)) {
      fullW0 = svd(Y0 %*% t(Y0))$u[, 1:(m - p - q), drop = FALSE]
    }
    W0 = fullW0[, 1:K1, drop = FALSE]
    alpha = solve(t(W0) %*% W0) %*% t(W0) %*% Y0
    bycx = solve(t(X) %*% X) %*% t(X) %*% Yc
    ac = alpha[, ctl, drop = FALSE]
    sizeratios = getsizeratios(bycx, ac, method = method)
    if (is.null(cutoff)) 
      cutoff = getKcutoff(m, n)
    keep = sizeratios > cutoff
    K = sum(keep)
    return(list(k = K, cutoff = cutoff, sizeratios = sizeratios, 
                fullW0 = fullW0))
  }

getsizeratios <- function (bycx, ac, method = "fast") 
  {
    K1 = nrow(ac)
    nc = length(bycx)
    if (method == "fast") {
      bwx = bycx %*% t(ac) %*% solve(ac %*% t(ac))
      sizeratios = rep(0, K1)
      for (i in 1:K1) {
        divisor = sqrt(1 + sum(bwx[1:i]^2))
        betachat = bycx - bwx[, 1:i, drop = FALSE] %*% ac[1:i, 
                                                          , drop = FALSE]
        scaledbetac = sqrt(nc/(nc - i)) * betachat/divisor
        sizeratios[i] = median(abs(as.vector(ac[i, ])/as.vector(scaledbetac)))
      }
      return(sizeratios)
    }
    else if (method == "leave1out") {
      bycx = as.vector(bycx)
      A = Aj = ajtB = bwx = t(bycx) %*% t(ac)
      B = Bj = ac %*% t(ac)
      sizeratios = rep(0, K1)
      scaledbetac = rep(0, nc)
      temp = .Call("getsizeratios", as.double(A), as.double(B), 
                as.double(Aj), as.double(Bj), as.double(bycx), as.double(ac), 
                as.double(ajtB), as.double(sizeratios), as.double(scaledbetac), 
                as.double(bwx), as.integer(K1), as.integer(nc))
      return(temp[[8]])
    }
    else return(0)
  }


lsvd <- function (Y) 
  {
    temp = svd(Y %*% t(Y))
    return(list(u = temp$u, d = sqrt(temp$d)))
  }


residop <- function (A, B) 
  {
    return(A - B %*% solve(t(B) %*% B) %*% t(B) %*% A)
  }


tological <- function (ctl, n) {
    ctl2 = rep(FALSE, n)
    ctl2[ctl] = TRUE
    return(ctl2)
  }

inputcheck1 <- function (Y, X, Z, ctl, check.na = TRUE) {
    m = nrow(Y)
    n = ncol(Y)
    if (m > n) 
      warning("m is greater than n!  This is not a problem itself, but may indicate that you need to transpose your data matrix.  Please ensure that rows correspond to observations (e.g. microarrays) and columns correspond to features (e.g. genes).")
    if (check.na & sum(is.na(Y)) > 0) 
      warning("Y contains missing values.  This is not supported.")
    if (sum(Y == Inf, na.rm = TRUE) + sum(Y == -Inf, na.rm = TRUE) > 
        0) 
      warning("Y contains infinities.  This is not supported.")
    XZ = cbind(X, Z)
    d = svd(t(XZ) %*% XZ, nu = 0, nv = 0)$d
    if (d[1]/d[length(d)] > 10^10) 
      warning("There appears to be linear dependence between the columns of X and Z.")
    if (sum(ctl) == 0) 
      warning("No genes are defined as control genes.  This is not supported.")
    return(NULL)
}

projectionplotvariables <- function (Y, X, W) {
    W0 = residop(W, X)
    bwx = solve(t(X) %*% X) %*% t(X) %*% W
    svdw0 = svd(W0)
    vd = t((1/svdw0$d) * t(svdw0$v))
    W0 = W0 %*% vd
    bwx = bwx %*% vd
    W0Y = t(W0) %*% Y
    u = svd(W0Y %*% t(W0Y))$u
    W0 = W0 %*% u
    bwx = bwx %*% u
    projectionplotalpha = t(W0) %*% Y
    byx = solve(t(X) %*% X) %*% t(X) %*% Y
    projectionplotW = W0 + X %*% bwx
    return(list(byx = byx, bwx = bwx, projectionplotalpha = projectionplotalpha, 
                projectionplotW = projectionplotW))
  }

RUV1 <- function (Y, eta, ctl, include.intercept = TRUE) {
  if (is.null(eta)) 
    return(Y)
  m = nrow(Y)
  n = ncol(Y)
  ctl = tological(ctl, n)
  if (is.numeric(eta)) 
    if (length(eta) == 1) 
      eta = matrix(1, n, 1)
  if (is.matrix(eta)) 
    if (nrow(eta) != n) 
      eta = t(eta)
  eta = design.matrix(eta, name = "eta", include.intercept = include.intercept)
  eta = t(eta)
  Yc = Y[, ctl, drop = FALSE]
  etac = eta[, ctl, drop = FALSE]
  if (sum(is.na(Yc)) == 0) 
    return(Y - Yc %*% t(etac) %*% solve(etac %*% t(etac)) %*% 
             eta)
  else {
    for (i in 1:m) {
      keep = !is.na(Yc[i, ])
      Yci = Yc[i, keep, drop = FALSE]
      etaci = etac[, keep, drop = FALSE]
      Y[i, ] = Y[i, ] - Yci %*% t(etaci) %*% solve(etaci %*% 
                                                     t(etaci)) %*% eta
    }
    return(Y)
  }
}

RUV2 <- function (Y, X, ctl, k, Z = 1, eta = NULL, include.intercept = TRUE, 
            fullW = NULL, svdyc = NULL, do_projectionplot = TRUE, inputcheck = TRUE)  { # Y is a matrix of sample X metabolite
    if (is.data.frame(Y)) 
      Y = data.matrix(Y)
    m = nrow(Y)
    n = ncol(Y)
    X = rX = design.matrix(X, include.intercept = FALSE)
    p = ncol(X)
    if (is.numeric(Z)) 
      if (length(Z) == 1) 
        Z = matrix(1, m, 1)
    if (!is.null(Z)) {
      Z = design.matrix(Z, name = "Z", include.intercept = include.intercept)
      q = ncol(Z)
    }
    else q = 0
    ctl = tological(ctl, n)
    if (inputcheck) 
      inputcheck1(Y, X, Z, ctl)
    if (k > sum(ctl)) 
      stop("k must not be larger than the number of controls")
    Y = RUV1(Y, eta, ctl, include.intercept = include.intercept)
    if (q > 0) {
      Y = residop(Y, Z)
      X = residop(X, Z)
    }
    if (is.null(fullW)) {
      if (is.null(svdyc)) 
        svdyc = lsvd(Y[, ctl, drop = FALSE])
      fullW = svdyc$u[, 1:min((m - p - q), sum(ctl)), drop = FALSE]
    }
    W = alpha = byx = bwx = projectionplotW = projectionplotalpha = NULL
    if (k > 0) {
      W = fullW[, 1:k, drop = FALSE]
      XZW = cbind(X, Z, W)
      if (do_projectionplot) {
        ppvars = projectionplotvariables(Y, X, W)
        byx = ppvars$byx
        bwx = ppvars$bwx
        projectionplotalpha = ppvars$projectionplotalpha
        projectionplotW = ppvars$projectionplotW
      }
    }
    else {
      XZW = cbind(X, Z)
    }
    A = solve(t(XZW) %*% XZW)
    AXZW = A %*% t(XZW)
    betagammaalphahat = AXZW %*% Y
    resids = Y - XZW %*% betagammaalphahat
    betahat = betagammaalphahat[1:p, , drop = FALSE]
    if (k > 0) 
      alpha = betagammaalphahat[(p + q + 1):(p + q + k), , 
                                drop = FALSE]
    multiplier = as.matrix(diag(A)[1:p])
    df = m - p - q - k
    sigma2 = apply(resids^2, 2, sum)/df
    sigma2 = as.vector(sigma2)
    se = sqrt(multiplier %*% t(sigma2))
    tvals = betahat/se
    pvals = tvals
    for (i in 1:nrow(pvals)) pvals[i, ] = 2 * pt(-abs(tvals[i, 
                                                            ]), df)
    Fstats = apply(betahat * (solve(AXZW[1:p, , drop = FALSE] %*% 
                                      t(AXZW[1:p, , drop = FALSE])) %*% betahat), 2, sum)/p/sigma2
    Fpvals = pf(Fstats, p, df, lower.tail = FALSE)
    return(list(betahat = betahat, sigma2 = sigma2, t = tvals, 
                p = pvals, Fstats = Fstats, Fpvals = Fpvals, multiplier = multiplier, 
                df = df, W = W, alpha = alpha, byx = byx, bwx = bwx, 
                X = rX, k = k, ctl = ctl, Z = Z, eta = eta, fullW = fullW, 
                projectionplotW = projectionplotW, projectionplotalpha = projectionplotalpha, 
                include.intercept = include.intercept, method = "RUV2"))
  }

RUV4 <- function (Y, X, ctl, k, Z = 1, eta = NULL, include.intercept = TRUE, 
            fullW0 = NULL, inputcheck = TRUE) 
  { # Y is a matrix of sample X metabolite
    if (is.data.frame(Y)) 
      Y = data.matrix(Y)
    m = nrow(Y)
    n = ncol(Y)
    X = rX = design.matrix(X, include.intercept = FALSE)
    p = ncol(X)
    if (is.numeric(Z)) 
      if (length(Z) == 1) 
        Z = matrix(1, m, 1)
    if (!is.null(Z)) {
      Z = design.matrix(Z, name = "Z", include.intercept = include.intercept)
      q = ncol(Z)
    }
    else q = 0
    ctl = tological(ctl, n)
    if (inputcheck) 
      inputcheck1(Y, X, Z, ctl)
    Y = RUV1(Y, eta, ctl, include.intercept = include.intercept)
    if (q > 0) {
      Y = residop(Y, Z)
      X = residop(X, Z)
    }
    Y0 = residop(Y, X)
    if (is.null(fullW0)) {
      fullW0 = svd(Y0 %*% t(Y0))$u[, 1:(m - p - q), drop = FALSE]
    }
    if (k > 0) {
      W0 = fullW0[, 1:k, drop = FALSE]
      alpha = solve(t(W0) %*% W0) %*% t(W0) %*% Y0
      byx = solve(t(X) %*% X) %*% t(X) %*% Y
      bycx = byx[, ctl, drop = FALSE]
      alphac = alpha[, ctl, drop = FALSE]
      bwx = bycx %*% t(alphac) %*% solve(alphac %*% t(alphac))
      W = W0 + X %*% bwx
      XZW = cbind(X, Z, W)
    }
    else {
      XZW = cbind(X, Z)
      W = alpha = byx = bwx = NULL
    }
    A = solve(t(XZW) %*% XZW)
    AXZW = A %*% t(XZW)
    betagammaalphahat = AXZW %*% Y
    resids = Y - XZW %*% betagammaalphahat
    betahat = betagammaalphahat[1:p, , drop = FALSE]
    multiplier = as.matrix(diag(A)[1:p])
    df = m - p - q - k
    sigma2 = apply(resids^2, 2, sum)/df
    sigma2 = as.vector(sigma2)
    se = sqrt(multiplier %*% t(sigma2))
    tvals = betahat/se
    pvals = tvals
    for (i in 1:nrow(pvals)) pvals[i, ] = 2 * pt(-abs(tvals[i, 
                                                            ]), df)
    Fstats = apply(betahat * (solve(AXZW[1:p, , drop = FALSE] %*% 
                                      t(AXZW[1:p, , drop = FALSE])) %*% betahat), 2, sum)/p/sigma2
    Fpvals = pf(Fstats, p, df, lower.tail = FALSE)
    return(list(betahat = betahat, sigma2 = sigma2, t = tvals, 
                p = pvals, Fstats = Fstats, Fpvals = Fpvals, multiplier = multiplier, 
                df = df, W = W, alpha = alpha, byx = byx, bwx = bwx, 
                X = rX, k = k, ctl = ctl, Z = Z, eta = eta, fullW0 = fullW0, 
                include.intercept = include.intercept, method = "RUV4"))
  }


#### EigenMS ####

# EigenMS normalization
# Ref: "Normalization of peak intensities in bottom-up MS-based proteomics using
#       singular value decomposition" Karpievitch YV, Taverner T, Adkins JN,
#       Callister SJ, Anderson GA, Smith RD, Dabney AR. Bioinformatics 2009
# and
# Ref: "Metabolomics data normalization with EigenMS"
#       Karpievitch YK, Nikolic SB, Wilson R, Sharman JE, Edwards LM
#       Submitted to PLoS ONE
# 
# Here we allow multiple facotrs to be preserved in ANOVA model before identifying bais. 
# This requires a bit of thought on the side of the researchers, only few factors should be 'preserved'.
# For example 'treatment group' is important to preserve but 'age' may or may not be important to preserve.
# Here we do not utilize peptides with 1+ grp missing completely, this is a separate problem 
# addressed by Wang et al, 2011
#
# Written by Yuliya Karpievitch, Tom Taverner, and Shelley Herbrich 
# email: yuliya.k@gmail.com
#
# Version has been split into 2 main functions: eig_norm1 finds significant bias trends,
# then the user can decide if he/she wants to use that number.
# For Matabolomics data we suggest useing 20% of the number of samples examined 
# as the number of bias trends to eliminated. By setting 
# ints_eig1$h.c = ceil(.2 * num_samples)
# if the number of automatically identified bias trends are close to that number we suggest 
# using the estimated number of bias trends.
# 
# eig_norm2 function normalizes the data
# Note that rescaling has been abandoned in the latest version due to discussions abot the fact 
# that systematic bias is what has been inadvertently added, thus we need to remove the bias 
# but nto add any additional variation. 
# Metabolomics in particular has a lot of variation that still remains that we concluded 
# there is not need for rescaling. 
#
# EigenMS estimates and preserves fixed effects. 
# Contributions to using Random effects are welcome. 

# HOWTO RUN: 
# source('~/EigenMS/EigenMS.R')
# dd = # read in the data
# grps = # read in group information file 
# logInts = # subset the dd matrix to only hte portion that contains intensities
#           # replace 0's with NAs, do a log-transformation to approximate Normality.
# prot.info = cbind(data.frame(rownames(dd)), data.frame(rownames(dd))) # peptideIDs, no PR IDs here, just duplicate column
#             # in case od protein IDs, those are not importnat for normalization and will be ignored.  
# ints_eig1 = eig_norm1(m=logInts, treatment=grps, prot.info=prot.info)
# ints_norm = eig_norm2(rv=ints_eig1) 
#
# eig_norm1 = function(m, treatment, prot.info, write_to_file='')
#       m - m x n matrix of log-transformed intensities, num peptides x num samples
#       treatment - either a single factor or a data.frame of factors (actual R factors, else code will fail) 
#                   eg:  bl = gl(3,14,168) # 3 factors, repeat every 14 times, total of 168 samples
#                        it = gl(2,7,168)  # 2 factors, repeat every 7 samples, 168 total
#                        Pi = gl(2,42,168) # 2 factors, repeat twice: 42 of PI+, PI-, PI+, PI-
#                       grpFactors = data.frame(bl,it,Pi)# factors we would like to preserve in EigenMS
#       prot.info - 2 column data frame with peptide and protein IDs in that order.
#                   for metabolites both columns should contain metabolite IDs.
#       write_to_file -  if a string is passed in, 'complete' peptides will be written to that file name
#                       Some peptides could be eliminated due to too many missing values (if not able to do ANOVA)                      

# First portion of EigenMS: identify bias trends 
eig_norm1 = function(m, treatment, prot.info, write_to_file=''){
  # Identify significant eigentrends, allow the user to adjust the number (with causion! if desired)
  # before normalizing with eig_norm2
  # 
  # Input:
  #   m: An m x n (peptides x samples) matrix of expression data, log-transformed!
  #      peptide and protein identifiers come from the get.ProtInfo()
  #   treatment:  either a single factor indicating the treatment group of each sample i.e. [1 1 1 1 2 2 2 2...]
  #               or a frame of factors:  treatment= data.frame(cbind(data.frame(Group), data.frame(Time)) 
  #   prot.info: 2+ colum data frame, pepID, prID columns IN THAT ORDER. 
  #              IMPORTANT: pepIDs must be unique identifiers and will be used as Row Names 
  #              If normalizing non-proteomics data, create a column such as: paste('ID_',seq(1:num_rows), sep='')
  #              Same can be dome for ProtIDs, these are not used for normalization but are kept for future analyses 
  #   write_to_file='' - if a string is passed in, 'complete' peptides (peptides with NO missing observations)
  #              will be written to that file name
  #                    
  # Output: list of:
  #   m, treatment, prot.info, grp - initial parameters returned for futre reference 
  #   my.svd - matrices produced by SVD 
  #   pres - matrix of peptides that can be normalized, i.e. have enough observations for ANOVA, 
  #   n.treatment - number of factors passed in
  #   n.u.treatment - number of unique treatment facotr combinations, eg: 
  #                   Factor A: a a a a c c c c
  #                   Factor B: 1 1 2 2 1 1 2 2
  #                   then:  n.treatment = 2; n.u.treatment = 4
  #   h.c - bias trends 
  #   present - names/IDs of peptides on pres
  #   complete - complete peptides, no missing values, these were used to compute SVD
  #   toplot1 - trends automatically produced, if one wanted to plot at later time. 
  #   Tk - scores for each bias trend 
  #   ncompl - number of complete peptides with no missing observations
  print("Data dimentions: ")  
  print(dim(m))
  # check if treatment is a 'factor' vs data.frame', i.e. single vs multiple factors
  if(class(treatment) == "factor") { # TRUE if one factor
    n.treatment = 1 # length(treatment)
    n.u.treatment = length(unique(treatment))[1]
  } else { # data.frame
    n.treatment = dim(treatment)[2]
    n.u.treatment = dim(unique(treatment))[1] # all possible tretment combinations
  }
  # convert m to a matrix from data.frame
  m = as.matrix(m) # no loss of information
  
  # filter out min.missing, here just counting missing values
  # if 1+ treatment completely missing, cannot do ANOVA, thus cannot preserve grp diff.
  # IMPORTANT: we create a composite grp = number of unique combinations of all groups, only for 
  # 'nested' groups for single layer group is left as it is 
  grpFactors = treatment # temporary var, leftover from old times...
  
  nGrpFactors = n.treatment # length(colnames(treatment)) # not good: dim(grpFactors)
  if(nGrpFactors > 1) { # got nested factors
    ugrps = unique(grpFactors)
    udims = dim(ugrps)
    grp = NULL
    for(ii in 1:udims[1]) {
      pos = grpFactors[,1] == ugrps[ii,1] # set to initial value
      for(jj in 2:udims[2]) { 
        pos = pos & grpFactors[,jj] == ugrps[ii,jj]
      }
      grp[pos] = rep(ii, sum(pos))
    }
    grp = as.factor(grp)
  } else {
    grp = treatment
  }
  nobs = array(NA, c(nrow(m), length(unique(grp)))) # noobs = number of observations 
  
  print('Treatmenet groups:')
  print(grp)
  
  for(ii in 1:nrow(m)) {
    for(jj in 1:length(unique(grp))) {
      nobs[ii,jj] = sum(!is.na(m[ii, grp==unique(grp)[jj]])) # total number of groups num(g1) * num(g2) * ...
    } 
  } 
  # now 'remove' peptides with missing groups
  present.min = apply(nobs, 1, min) # number present in each group
  ii = present.min == 0   # 1+ obs present in ALL of the groups
  nmiss = sum(present.min == 0) # not used, one value of how many peptides have 1+ grp missing completely
  pmiss = rbind(m[ii,]) # these have 1+ grp missing !!!!
  # rownames must be UNIQUE, if have possible duplicates: use 'ii' ?
  rownames(pmiss) = prot.info[ii,1]  # set rownames, 
  
  # create matrix for peptides with enough observations for ANOVA
  # 'present' are names of the peptides (pepID) and 'pres' are abundances
  # NOTE: ! negates the proteins, so we get ones that have 1+ obs in each group 
  present = prot.info[which(!prot.info[,1] %in% rownames(pmiss)), ] # rownames OK
  # pres = m[which(!rownames(m) %in% rownames(pmiss)), ]
  pres = m[which(!prot.info[,1] %in% rownames(pmiss)), ] # is this OK?
  rownames(pres) = prot.info[which(!prot.info[,1] %in% rownames(pmiss)),1]
  
  print('Selecting complete peptides')
  # Should issue an error message if we have NO complete peptides.
  # select only 'complete' peptides, no missing values
  nobs = array(NA, nrow(pres)) # reassign noobs to dims of 'present' 
  numiter = nrow(pres)
  for (ii in 1:numiter) {
    # if(ii %% 100 == 0) { print(ii) }
    nobs[ii] = sum(!is.na(pres[ii,]))
  }
  
  iii = nobs == ncol(pres)
  complete = rbind(pres[iii,])
  
  #  write out a file of complete peptides if file name is passed in
  if(write_to_file != '') {
    write.table(complete, file = write_to_file, append = FALSE,
                quote = FALSE, sep = "\t",
                eol = "\n", na = "NaN", dec = ".", row.names = TRUE,
                col.names = TRUE, qmethod = c("escape", "double"))
  }
  
  # compute bias with 'complete' matrix and residuals from 'present' 
  # calculate eigenpeptides for 'complete' data only
  # if have only 1 group, we do not need to preserve group differernces, everything is the same group, ex: QC samples
  # contrasts will fail if have only 1 group, thus have else
  if(n.u.treatment > 1) { 
    print('Got 2+ treatment grps')
    # check to see if we have multiple factors
    grpdim = dim(treatment)
    
    lm.fm = makeLMFormula(treatment, 'TREAT') # using general function that can accomodate for 1+ number of factors
    TREAT = treatment
    TREAT = data.frame(treatment) # temp var to work if we got only 1 treatment vector.
    if(class(treatment) == "factor") {
      colnames(TREAT) = "TREAT"
    } else {
      colnames(TREAT) = colnames(treatment)
    }     
    attach(TREAT)
    
    mod.c = model.matrix(lm.fm$lm.formula, data=TREAT, eval(parse(text=lm.fm$lm.params))) 
    Y.c = as.matrix(complete)
    options(warn = -1)
    
    # use lm() to get residuals
    formula1 = paste('t(Y.c)~', as.character(lm.fm$lm.formula)[2], sep = '')
    TREAT = treatment
    fit_lmAll = lm(eval(parse(text=formula1)))
    R.c = residuals(fit_lmAll)  # Oct 2 messing with residuals...
  } else {  # 1 group only, set residuals to original matrix
    print('Got 1 treatment grp')
    mod.c = as.numeric(t(treatment))
    R.c = t(as.matrix(complete))  # needs to be transposed to match the matrix returned from lm
    TREAT = treatment
  }
  
  print('Computing SVD, estimating Eigentrends...') # let user know what is going on
  # residuals are centered around 0, here center samples not peptides/metabolites
  # centering is basic normalization
  
  R.c_center = scale(R.c, center = TRUE, scale = FALSE)  # t(scale(t(R.c), center = TRUE, scale = FALSE))
  my.svd = svd(R.c_center)  # can use wrapper below to chek if SVD has a problem...
  temp = my.svd$u
  my.svd$u = my.svd$v
  my.svd$v = temp
  
  #identify number of eigenvalues that account for a significant amount of residual variation
  numcompletepep = dim(complete)[1] # save to return to the user as part of the return list  
  # this is important info for publications
  # tell users how many peptides/metabolites the trends are based on
  # can also be determined by doing dim(return_value_fromEIg_norm1$pres)
  
  print(paste('Number of treatments: ', n.u.treatment))
  h.c = sva.id(complete, treatment, n.u.treatment, lm.fm=lm.fm,seed=1234)$n.sv
  print(paste("Number of significant eigenpeptides/trends", h.c) )
  
  # show RAW trends
  # center each peptide around zero (subtract its mean across samples)
  complete_center = scale(t(complete), center = TRUE, scale = FALSE)
  print('Preparing to plot...')
  
  n.u.treatment
  toplot1 = svd(complete_center) # scales above
  temp = toplot1$u
  toplot1$u = toplot1$v
  toplot1$v = temp
  
  par(mfcol=c(3,2))
  par(mar = c(2,2,2,2))
  #plot.eigentrends(toplot1, "Raw Data")
  #plot.eigentrends(my.svd, "Residual Data")
  
  d = my.svd$d;  ss = d^2;
  Tk = signif(ss/sum(ss)* 100, 2)
  
  retval = list(m=m, treatment=treatment, my.svd=my.svd,
                pres=pres, n.treatment=n.treatment, n.u.treatment=n.u.treatment,
                h.c=h.c, present=present, prot.info=prot.info,
                complete=complete, toplot1=toplot1, Tk=Tk, ncompl=numcompletepep,
                grp=grp) 
  return(retval)
}


# Second portion of EigenMS: remove effects of bias 
# split into 2 functions to allo the user to examine the bias trends and be able to change the number
# by resetting h.c returned by eig_norm1(return_value_from_eig_norm1)
eig_norm2 = function(rv) { 
  # UNPUT:
  #   rv - return value from the eig_norm1
  #   if user wants to change the number of bias trends that will be eliminated h.c in rv should 
  #   be updates to the desired number
  # 
  # OUTPUT: 
  #   normalized - matrix of normalized abundances with 2 columns of protein and peptdie names
  #   norm_m - matrix of normalized abundances, no extra columns  
  #   eigentrends - found in raw data, bias trendsup to h.c
  #   rescrange - rescaling range for the addition of the while noise to avoid overfitting 
  #   norm.svd - trends in normalized data, if one wanted to plot at later time. 
  #   exPeps - excluded peptides - excluded due to exception in fitting a linear model
  
  m = rv$pres # yuliya: use pres matrix, as we cannot deal with m anyways, need to narrow it down to 'complete' peptides
  treatment = rv$treatment
  my.svd = rv$my.svd
  pres = rv$pres
  n.treatment = rv$n.treatment
  n.u.treatment = rv$n.u.treatment 
  numFact = dim(rv$treatment)[2]
  print(paste('Unique number of treatment combinations:', n.u.treatment) )
  h.c = rv$h.c
  present = rv$present
  toplot1 = rv$toplot1
  # vector of indicators of peptides that threw exeptions 
  exPeps = vector(mode = "numeric", length = nrow(pres))
  
  print("Normalizing...")
  treatment = data.frame(treatment) # does this need to be done?
  if(n.u.treatment > 1) {
    lm.fm = makeLMFormula(treatment, 'ftemp')
    mtmp = model.matrix(lm.fm$lm.formula, data=treatment, eval(parse(text=lm.fm$lm.params)))  #contrasts=list(bl="contr.sum", it="contr.sum",Pi="contr.sum", tp="contr.sum"))
  } else {  # have 1 treatment group
    mtmp = treatment # as.numeric(t(treatment)) 
  }
  # above needed to know how many values will get back for some matrices
  # create some variables:
  betahat = matrix(NA,nrow=dim(mtmp)[2],ncol=nrow(pres)) 
  newR = array(NA, c(nrow(pres), ncol(pres))) #, n.treatment))
  norm_m = array(NA, c(nrow(pres), ncol(pres))) # , n.treatment))
  numsamp = dim(pres)[2]
  numpep = dim(pres)[1]
  betahat_n = matrix(NA,nrow=dim(mtmp)[2],ncol=nrow(pres))
  rm(mtmp) 
  
  V0 = my.svd$v[,1:h.c,drop=F]   # residual eigenpeptides
  
  if(n.u.treatment == 1) { # got 1 treatment group
    for (ii in 1:nrow(pres)) {
      if(ii%%250 == 0) { print(paste('Processing peptide ',ii))  }
      pep = pres[ii, ] 
      pos = !is.na(pep)
      peptemp = as.matrix(pep[pos]) # take only the observed values
      resm = rep(NA, numsamp) 
      resm[pos] = as.numeric(pep[pos])
      bias = array(NA, numsamp)
      bias[pos] = resm[pos] %*% V0[pos,] %*% t(V0[pos,])
      norm_m[ii, ] = as.numeric(pep - bias)
    }
    
  } else { # got 2+ treatment groups
    for (ii in 1:nrow(pres)) {
      if(ii %% 100 == 0) { print(paste('Processing peptide ',ii))  }
      pep = pres[ii, ] 
      pos = !is.na(pep)
      peptemp = as.matrix(pep[pos]) # take only the observed values, may not be needed in R? but this works
      ftemp = treatment[pos,]
      ftemp = data.frame(ftemp)
      #### use try, not entirely sure if need for modt, need it for solve lm?!
      options(warn = -1)
      lm.fm = makeLMFormula(ftemp, 'ftemp') # using general function that can accomodate for 1+ number of factors
      modt = try(model.matrix(lm.fm$lm.formula, data=ftemp, eval(parse(text=lm.fm$lm.params))), silent=TRUE)
      options(warn = 0)
      
      if(!inherits(modt, "try-error")) { # do nothing if could not make model matrix
        options(warn = -1)
        # if we are able to solve this, we are able to estimate bias  
        bhat =  try(solve(t(modt) %*% modt) %*% t(modt) %*% peptemp)
        options(warn = 0)
        if(!inherits(bhat, "try-error")) {
          betahat[,ii] = bhat
          ceffects = modt %*% bhat  # these are the group effects, from estimated coefficients betahat
          
          resm = rep(NA, numsamp) # really a vector only, not m 
          resm[pos] = as.numeric(pep[pos] - ceffects)
          bias = array(NA, numsamp)
          bias[pos] = resm[pos] %*% V0[pos,] %*% t(V0[pos,])
          norm_m[ii, ] = as.numeric(pep - bias)
          
          # yuliya:  but newR should be computed on Normalized data
          resm_n = rep(NA, numsamp)
          bhat_n =  solve(t(modt) %*% modt) %*% t(modt) %*% norm_m[ii, pos]
          betahat_n[,ii] = bhat_n
          ceffects_n = modt %*% bhat_n
          resm_n[pos] = norm_m[ii,pos] - ceffects
          newR[ii, ] = resm_n
        } else {
          print(paste('got exception 2 at peptide:', ii, 'should not get here...')) 
          exPeps[ii] = 2 # should not get 2 here ever...
        }
      } else {
        print(paste('got exception at peptide:', ii)) 
        exPeps[ii] = 1 # keep track of peptides that threw exeptions, check why...
      }
    }
  } # end else - got 2+ treatment groups
  
  #####################################################################################
  # rescaling has been eliminated form the code after discussion that bias 
  # adds variation and we remove it, so no need to rescale after as we removed what was introduced
  y_rescaled = norm_m # for 1 group normalization only, we do not rescale
  # add column names to y-rescaled, now X1, X2,...
  colnames(y_rescaled) = colnames(pres) # these have same number of cols
  rownames(y_rescaled) = rownames(pres) 
  y_resc = data.frame(present, y_rescaled)  
  rownames(y_resc) = rownames(pres)  # rownames(rv$normalized)
  final = y_resc # row names are assumed to be UNIQUE, peptide IDs are unique
  
  # rows with all observations present
  complete_all = y_rescaled[rowSums(is.na(y_rescaled))==0,,drop=F]
  
  #  x11() # make R open new figure window
  par(mfcol=c(3,2))
  par(mar = c(2,2,2,2))
  # center each peptide around zero (subtract its mean across samples)
  # note: we are not changing matrix itself, only centerig what we pass to svd
  complete_all_center = t(scale(t(complete_all), center = TRUE, scale = FALSE))
  toplot3 = svd(complete_all_center)
  #plot.eigentrends(toplot1, "Raw Data")
  #plot.eigentrends(toplot3, "Normalized Data")
  
  print("Done with normalization!!!")
  colnames(V0) =  paste("Trend", 1:ncol(V0), sep="_")
  
  maxrange = NULL # no rescaling # data.matrix(maxrange)
  return(list(normalized=final, norm_m=y_rescaled, eigentrends=V0, rescrange=maxrange, 
              norm.svd=toplot3, exPeps=exPeps)) 
} # end function eig_norm2

# make a string formula to use in 'lm' call when computing grp differences to preserve
makeLMFormula = function(eff, var_name='') {
  # eff - effects used in contrasts
  # var_name - for singe factor use var-name that is passed in as variable names, otherwise it has no colnmae
  #           only used for a single factor
  if(is.factor(eff))
  {
    ndims = 1
    cols1 = var_name # ftemp in EigenMS
  }
  else
  {
    ndims = dim(eff)[2] 
    cols1 = colnames(eff)
  }
  lhs = cols1[1]
  lm.fm = NULL
  # check if can have a list if only have 1 factor...
  
  params = paste('contrasts=list(', cols1[1], '=contr.sum', sep=)
  
  if (ndims > 1) { # removed ndims[2] here, now ndims holds only 1 dimention...
    for (ii in 2:length(cols1))
    {
      lhs = paste(lhs, "+", cols1[ii])  # bl="contr.sum",
      params = paste(params, ',', cols1[ii], '=contr.sum', sep='')
    }
  }
  params = paste(params,")") 
  lm.formula = as.formula(paste('~', lhs))
  lm.fm$lm.formula = lm.formula
  lm.fm$lm.params = params
  return(lm.fm)
}	

plot.eigentrends.start = function(svdr, title1, pos1=1){
  # No check for valid range of pos1 is performed!!! 
  v = svdr$v
  d = svdr$d
  ss = d^2
  Tk = signif(ss/sum(ss)* 100, 2)
  #  pe = signif(d/sum(d, na.rm=T)*100, 2)
  titles = paste("Trend ", pos1:(pos1+3), " (", Tk[pos1:(pos1+3)], "%)", sep = "")
  do.text = function(j) mtext(titles[j], cex=0.7, padj=-0.7, adj=1)
  range.y = range(as.numeric(v[,pos1:(pos1+3)]), na.rm=T)
  
  toplot1_1 = as.numeric(v[,pos1])
  toplot1_2 = as.numeric(v[,(pos1+1)])
  toplot1_3 = as.numeric(v[,(pos1+2)])
  
  plot(c(1:length(toplot1_1)), toplot1_1, type='b', ann=F, ylim=range.y)
  do.text(1)
  abline(h=0, lty=3)
  title(title1, cex.main = 1.2, font.main= 1, col.main= "purple", ylab=NULL)
  plot(c(1:length(toplot1_2)), toplot1_2, type='b', ann=F, ylim=range.y)
  do.text(2)
  abline(h=0, lty=3)
  plot(c(1:length(toplot1_3)), toplot1_3, type='b', ann=F, ylim=range.y)
  do.text(3)
  abline(h=0, lty=3)
  return(Tk)
}

# not used in EigenMS 
make.formula.string = function(factors, do.interactions=FALSE){
  fs = "1"
  if(length(factors)){
    fs = paste(factors, collapse=" + ")
    if(do.interactions && length(factors) > 1)
      fs = paste(unlist(lapply(as.data.frame(t(combinations(length(factors), 2, factors)), stringsAsFactors=F), paste, collapse="*")), collapse = " + ")
  }
  return(fs)
}

sva.id = function(dat, treatment, n.u.treatment, lm.fm, B=500, sv.sig=0.05, seed=NULL) {
  # Executes Surrogate Variable Analysis
  # Input:
  #   dat: A m peptides/genes by n samples matrix of expression data
  #   mod: A model matrix for the terms included in the analysis 
  #   n.u.treatment - 0 or 1, if we are normalizing data with NO groups or some groups, QC vs samples
  #   B: The number of null iterations to perform
  #   sv.sig: The significance cutoff for the surrogate variables
  #   seed: A seed value for reproducible results
  # Output
  #    n.sv: Number of significant surrogate variables. 
  #    id: An indicator of the significant surrogate variables
  #    B: number of permutation to do
  #    sv.sig: significance level for surrogate variables
  print("Number of complete peptides (and samples) used in SVD")
  print(dim(dat))
  
  
  if(!is.null(seed))  { set.seed(seed) }
  warn = NULL
  n = ncol(dat)
  m = nrow(dat)
  
  # ncomp = length(as.numeric(n.u.treatment))
  ncomp = n.u.treatment # JULY 2013: as.numeric(n.u.treatment)
  print(paste("Number of treatment groups (in svd.id): ", ncomp))
  # should be true for either case and can be used later
  
  if(ncomp > 1) { #   
    formula1 = paste('t(dat)~', as.character(lm.fm$lm.formula)[2], sep = '')
    fit_lmAll = lm(eval(parse(text=formula1)))
    res = t(residuals(fit_lmAll))
  } else {
    res = dat
  }
  # centering was not done before...
  # center each peptide around zero (subtract its mean across samples)
  # note: we are not changing matrix itself, only centerig what we pass to svd
  res_center = t(scale(t(res), center = TRUE, scale = FALSE))
  
  uu = svd(t(res_center)) # NEED a WRAPPER for t(). the diag is min(n, m)
  temp = uu$u
  uu$u = uu$v
  uu$v = temp
  
  
  # yuliya: Sept 2014: can I get around without using H?? 
  #  ndf = min(n, m) - ceiling(sum(diag(H)))  
  #  dstat = uu$d[1:ndf]^2/sum(uu$d[1:ndf]^2)
  #  dstat0 = matrix(0,nrow=B,ncol=ndf)
  #  s0 = diag(uu$d) # no need for diag here, should investigate why this is a vector already...
  s0 = uu$d
  s0 = s0^2
  dstat = s0/sum(s0)  # this did not have 'diag' in it in Tom's code...
  ndf = length(dstat) # sticking to Tom's variable name
  # print(paste('length(dstat) = ndf = ', ndf))
  dstat0 = matrix(0,nrow=B,ncol=ndf) # num samples (?)
  
  print("Starting Bootstrap.....")
  # this is the Bootstrap procedure that determines the number of significant eigertrends... 
  for(ii in 1:B){
    if(ii %% 50 == 0) { print(paste('Iteration ', ii)) }
    res0 = t(apply(res, 1, sample, replace=FALSE)) # regression
    # yuliya: not sure if this is needed at all
    # not needed for 1 group normalizaiton
    ##### res0 = res0 - t(H %*% t(res0))
    # yuliya: Sept 3, 2014: REMOVED above line. Do not think this needs to be done.. 
    # center each peptide around zero (subtract its mean across samples)
    # note: we are not changing matrix itself, only centerig what we pass to svd
    res0_center = t(scale(t(res0), center = TRUE, scale = FALSE))
    uu0 = svd(res0_center)
    temp = uu0$u  # why did tom do this??
    uu0$u = uu0$v
    uu0$v = temp
    
    ss0 = uu0$d  # no need for diag.... 
    ss0 = ss0^2
    dstat0[ii,] = ss0 / sum(ss0) # Tk0 in Matlab
  }
  
  # yuliya: check p-values here, Tom had mean value...
  psv = rep(1,n)
  for(ii in 1:ndf){
    # psv[ii] = mean(dstat0[,ii] >= dstat[ii])
    # should this be compared to a MEAN?  Should this be dstat0[ii,] ?
    posGreater = dstat0[,ii] > dstat[ii]
    psv[ii] = sum(posGreater) / B
  }
  
  # p-values for peptides have to be in monotonically increasing order, 
  # set equal to previous one if not the case
  for(ii in 2:ndf){
    if(psv[(ii-1)] > psv[ii]) {
      # psv[ii] = max(psv[(ii-1)],psv[ii]) 
      psv[ii] = psv[(ii-1)] 
    }
  }
  nsv = sum(psv <= sv.sig)
  # tom - this should be at least 1
  # nsv = min(sum(psv <= sv.sig), 1, na.rm=T)
  return(list(n.sv = nsv,p.sv=psv))
}
# end sva.id


#### Batch Normalizer ####

### This contains an R function for normalizing GC/MS metabolomics data using Batch Normalizer.
### The Batch Normalizer function was originally published by Wang et al. (2012) in Analytical Chemistry.
### The function below is based on the algorithm as described in the published manuscript.

### Metabs missing >80% of data within phenotype will be excluded.


### batchNorm 

### CONTROL.DATA.LONG is a data frame of observed metabolite values for QC and analytical samples
### CONTROL.DATA.LONG must also include columns for id, batch, and seq1 (run order)
### id will be used as row names for CONTROL.DATA.LONG

### METAB.ID is a vector of metabolite names to be normalized and
### must be a subset of column names of CONTROL.DATA.LONG

### A.PHENOS is a character vector indicating analytical sample types

### QC.PHENOS is a character vector indicating QC sample types
### QC.PHENOS needs to be the same length as A.PHENOS, and order should correspond to A.PHENOS
### QC.PHENOS indicates which QC types correspond to analytical sample type if multiple QCs are involved


batchNormFUN<-function(CONTROL.DATA.LONG, METAB.ID, A.PHENOS, QC.PHENOS){
  
  names(A.PHENOS)<-QC.PHENOS
  
  names(CONTROL.DATA.LONG)<-tolower(names(CONTROL.DATA.LONG))  
  rownames(CONTROL.DATA.LONG)<-CONTROL.DATA.LONG[, "id"]
  metabs<-colnames(CONTROL.DATA.LONG)[grep(METAB.ID, colnames(CONTROL.DATA.LONG))]
  metabdata<-CONTROL.DATA.LONG[, metabs] 
  
  metabdata.t<-as.data.frame(t(metabdata))
  
  normDataList<-lapply(unique(QC.PHENOS), function(z){
    
    qc<-grep(z, rownames(CONTROL.DATA.LONG), value=TRUE)
    a<-grep(paste(A.PHENOS[which(names(A.PHENOS)==z)], collapse="|"), rownames(CONTROL.DATA.LONG), value=TRUE)
    
    ### eq1
    # get total abundance, TAq, for each control sample by adding all of the log2 peak intensities for each qc sample
    TA.q<-apply(metabdata.t[, qc], 2, function(x) sum(x, na.rm=TRUE))
    
    
    # get median total abundance across all qc and real data samples, median(TA*)
    sums<-apply(metabdata.t[, c(qc, a)], 2, function(x) sum(x, na.rm=TRUE))
    
    med.TA.<-median(sums, na.rm=TRUE)
    
    
    # get scaling factor, TA.scale.q, for each control sample
    TA.scale.q<-(med.TA./TA.q)
    
    
    # multiply control-sample-specific scaling factors by observed log2 peak areas for each metab within control sample for Yp,q
    scaled<-data.frame(mapply(`*`, metabdata.t[, qc], TA.scale.q, SIMPLIFY=FALSE))
    colnames(scaled)<-colnames(metabdata.t[, qc])
    rownames(scaled)<-rownames(metabdata.t[, qc])
    
    
    # using Yp,q as outcome variables, create separate linear regression models for each metabolite and each batch.
    # treat the injection order of the control sample as a linear predictor
    scaled.t<-as.data.frame(t(scaled))
    scaled.t$id<-rownames(scaled.t)
    
    batch.pos<-CONTROL.DATA.LONG[qc, c("id", "batch", "seq1")]
    
    merged<-merge(batch.pos, scaled.t, by="id")
    rownames(merged)<-merged$id
    merged$seq1<-as.numeric(merged$seq1)
    
    ### select metab cols
    pcs<-colnames(merged)[grep(METAB.ID, colnames(merged))]
    
    nbatches<-length(unique(CONTROL.DATA.LONG[, "batch"]))
    
    ### create dfs for the resultant values
    results.base<-as.data.frame(matrix(NA, nrow=length(pcs), ncol=nbatches))
    rownames(results.base)<-pcs
    colnames(results.base)<-paste(1:nbatches, sep="")
    results.R<-results.base
    
    for (j in 1:length(pcs)){
      
      myMetab<-pcs[j]
      
      for (i in 1:nbatches){
        
        subset<-merged[which(merged$batch==i), ]
        
        subset$seq1<-as.numeric(subset$seq1)
        
        formula<-as.formula(paste(myMetab, "~ seq1", sep=""))
        
        if (sum(is.na(subset[, myMetab]))>=2){
          
          results.base[myMetab, i]<-NA
          results.R[myMetab, i]<-NA
          
        } else {
          
          summary<-summary(lm(formula, data=subset, na.action=na.exclude))$coefficients
          
          base.pb<-summary["(Intercept)", "Estimate"]
          R.pb<-summary["seq1", "Estimate"]
          
          results.base[myMetab, i]<-base.pb
          results.R[myMetab, i]<-R.pb
          
        }
        
      }
      
    }  
    
    
    ### eq2
    # get total abundance, TA.s, for each real sample by adding up all the log2 peak intensities for each sample
    TA.s<-apply(metabdata.t[, c(qc, a)], 2, function(x) sum(x, na.rm=TRUE))
    
    
    # median(TA*) doesn't change
    
    
    # get scaling factor, TA.scale.s, for each real sample
    TA.scale.s<-(med.TA./TA.s)
    
    
    # get the median of peak p across all qc samples, median(C.p.qc.*)
    med.C.p.qc<-apply(metabdata.t[, qc], 1, function(x) median(x, na.rm=TRUE))
    
    
    # get vector of non-normalized abundances of peak p for all real samples for the numerator
    # C<-raw.t[, c(mQC, m01, m12)] ???
    
    
    # get part 1 of the numerator, Cp,s*TA.scale.s
    part1<-data.frame(mapply(`*`, metabdata.t[, c(qc, a)], TA.scale.s, SIMPLIFY=FALSE))
    colnames(part1)<-colnames(metabdata.t[, c(qc, a)])
    rownames(part1)<-rownames(metabdata.t[, c(qc, a)])
    
    
    # get full numerator, Cp,s*TA.scale.s*median(Cp,qc,*)
    part1.t<-as.data.frame(t(part1))
    numerator<-data.frame(mapply(`*`, part1.t, med.C.p.qc, SIMPLIFY=FALSE))
    colnames(numerator)<-colnames(part1.t)
    rownames(numerator)<-rownames(part1.t)
    numerator.t<-as.data.frame(t(numerator))
    
    
    # loop through numerators and calculate normalized values
    normalized<-as.data.frame(matrix(NA, nrow=dim(numerator.t)[1], ncol=dim(numerator.t)[2]))
    rownames(normalized)<-pcs
    colnames(normalized)<-colnames(numerator.t)
    samples<-colnames(numerator.t)
    batch.pos<-CONTROL.DATA.LONG[c(qc, a), c("id", "batch", "seq1")]
    
    for (j in 1:length(samples)){
      
      mySample<-samples[j]
      
      for (i in 1:length(pcs)){
        
        myMetab<-pcs[i]
        
        numerator.p.s<-numerator.t[myMetab, mySample]
        
        batch<-batch.pos[mySample, "batch"]
        
        seq1<-as.numeric(batch.pos[mySample, "seq1"])
        
        base<-results.base[myMetab, batch]
        
        R<-results.R[myMetab, batch]
        
        denominator<-(R*seq1)+base
        
        final<-numerator.p.s/denominator
        
        normalized[myMetab, mySample]<-final
        
      }
      
    }
    
    norm2<-as.data.frame(t(normalized))
    
    return(norm2)
    
  })
  
  norm<-do.call("rbind", normDataList)
  
  norm2<-norm[rownames(CONTROL.DATA.LONG), ]
  
  norm3<-cbind(CONTROL.DATA.LONG[, c("id", "batch", "seq1")], norm2)
  rownames(norm3)<-norm3$id
  
  return(norm3)
  
}


#### metaX ####

##' @title Normalisation of peak intensity
##' @description The normalize method performs normalisation on peak 
##' intensities.
##' @rdname normalize
##' @docType methods
##' @param para A metaXpara object.
##' @param method The normalization method: sum, vsn, quantiles, 
##' quantiles.robust, sum, pqn, combat. Default is sum.
##' @param valueID The name of the column which will be normalized.
##' @param norFactor The factor that will be used for normalization. This is 
##' usually used in urine data normalization. The factor is the column name in 
##' sample list file that contains the normalization factor value, such as the 
##' value "osmolality". Default the value is NULL. Only if the value of 
##' norFactor is not NULL and the parameter "method" is NULL, this 
##' normalization will work.
##' @param ... Additional parameter
##' @return A metaXpara object.
##' @author Bo Wen \email{wenbo@@genomics.cn}
##' @exportMethod
##' @examples
##' library(faahKO)
##' xset <- group(faahko)
##' xset <- retcor(xset)
##' xset <- group(xset)
##' xset <- fillPeaks(xset)
##' peaksData <- as.data.frame(groupval(xset,"medret",value="into"))
##' peaksData$name <- row.names(peaksData)
##' para <- new("metaXpara")
##' rawPeaks(para) <- peaksData
##' sampleListFile(para) <- system.file("extdata/faahKO_sampleList.txt", 
##'     package = "metaX")
##' para <- reSetPeaksData(para)
##' para <- metaX::normalize(para)

setClass("xcmsSet",
         representation = representation(peaks = "matrix",
                                         groups = "matrix",
                                         groupidx = "list",
                                         filled="numeric",
                                         phenoData = "data.frame",
                                         rt = "list",
                                         filepaths = "character",
                                         profinfo = "list",
                                         dataCorrection="numeric",
                                         polarity = "character",
                                         progressInfo = "list",
                                         progressCallback="function",
                                         mslevel = "numeric",
                                         scanrange = "numeric",
                                         .processHistory = "list"),
         prototype = prototype(peaks = matrix(nrow = 0, ncol = 0),
                               groups = matrix(nrow = 0, ncol = 0),
                               groupidx = list(),
                               filled = integer(0),
                               phenoData = data.frame(),
                               rt = list(),
                               filepaths = character(0),
                               profinfo = vector("list"),
                               dataCorrection=integer(0),
                               polarity = character(0),
                               progressInfo = list(),
                               mslevel = numeric(0),
                               scanrange= numeric(0),
                               progressCallback = function(progress) NULL,
                               .processHistory = list()),
         validity = function(object) {
           msg <- character()
           ## Check if all slots are present.
           slNames <- slotNames(object)
           missingSlots <- character()
           for (i in 1:length(slNames)) {
             if (!.hasSlot(object, slNames[i]))
               missingSlots <- c(missingSlots, slNames[i])
           }
           if (length(missingSlots) > 0)
             msg <- c(msg, paste0("This xcmsSet lacks slot(s): ",
                                  paste(missingSlots, collapse = ","),
                                  ". Please update the object using",
                                  " the 'updateObject' method."))
           ## Check the .processHistory slot.
           if (!any(missingSlots == ".processHistory")) {
             inh <- unlist(lapply(object@.processHistory,
                                  FUN = function(z) {
                                    return(inherits(z, "ProcessHistory"))
                                  }))
             if (!all(inh))
               msg <- c(msg,
                        paste0("Slot '.processHistory' should",
                               " only contain 'ProcessHistory'",
                               " objects!"))
           }
           if (length(msg))
             return(msg)
           return(TRUE)
         })


setClass("metaXpara", slots=c(
  ## input
  dir.case = "character",
  dir.ctrl = "character",
  sampleListFile = "character",
  sampleListHead = "character",
  sampleList = "data.frame",
  ratioPairs = "character",
  missValueImputeMethod = "character",
  pairTest = "logical",
  
  ## output
  outdir = "character",
  prefix = "character",
  xcmsPeakListFile = "character",
  fig = "list",
  peaksData = "data.frame",
  VIP = "data.frame",
  #sampleList = "data.frame",
  rawPeaks = "data.frame",
  xcmsSetObj = "xcmsSet", 
  quant = "data.frame",
  
  ## identification result
  idres = "data.frame",
  
  xcmsSet = representation(
    method = "character",                      
    ppm = "numeric",                     
    peakwidth = "numeric",                      
    snthresh = "numeric", 
    prefilter = "numeric",      
    mzCenterFun = "character",       
    integrate = "numeric",
    mzdiff = "numeric",
    noise = "numeric",
    verbose.columns = "logical",
    polarity = "character", 
    profparam = "list", 
    nSlaves = "numeric" ,
    fitgauss = "logical",
    sleep = "numeric",
    
    ## findPeaks.matchedFilter-methods 
    fwhm = "numeric",
    max = "numeric",
    step = "numeric"
    
    
  ),
  group = representation(
    ## The first round
    bw0 = "numeric",                
    mzwid0 = "numeric",
    
    ## The second round
    bw = "numeric",                
    mzwid = "numeric",             
    minfrac = "numeric",          
    minsamp = "numeric",        
    max = "numeric",          
    sleep = "numeric"
  ),
  
  retcor = representation(
    method = "character",
    profStep = "numeric",
    plottype = "character"
  ),
  
  qcRlscSpan = "numeric"
  
  
  
),
prototype = prototype(xcmsSet.method = "centWave",
                      xcmsSet.ppm = 10,
                      xcmsSet.peakwidth = c(5,15),
                      xcmsSet.snthresh = 6,
                      xcmsSet.prefilter = c(1,5000),
                      xcmsSet.mzCenterFun = "wMean",
                      xcmsSet.integrate = 1,
                      xcmsSet.mzdiff = -0.001,
                      xcmsSet.noise = 5000,
                      xcmsSet.verbose.columns = TRUE,
                      xcmsSet.polarity = "positive",
                      xcmsSet.profparam = list(step=0.005),
                      xcmsSet.nSlaves = 1,
                      xcmsSet.fitgauss = FALSE,
                      xcmsSet.sleep = 0,
                      
                      ## findPeaks.matchedFilter
                      xcmsSet.fwhm = 30,
                      xcmsSet.max = 5,
                      xcmsSet.step = 0.1,
                      
                      group.bw0 = 10,
                      group.mzwid0 = 0.015,
                      
                      group.bw = 5,
                      group.mzwid = 0.015,
                      group.minfrac = 0.3,
                      group.sleep = 0,
                      group.minsamp = 1,
                      group.max = 1000,
                      
                      retcor.method = "obiwarp",
                      retcor.plottype = "deviation",
                      retcor.profStep = 0.005,
                      
                      
                      outdir = "./",
                      prefix = "metaX",
                      peaksData = NULL,
                      VIP = NULL,
                      ratioPairs = NULL,
                      qcRlscSpan = 0,
                      missValueImputeMethod = "knn",#"svdImpute",
                      quant = NULL,
                      idres = NULL,
                      pairTest = FALSE,
                      
                      sampleListHead = c(sample="sample",
                                         order="order",
                                         batch="batch",
                                         class="class",
                                         isQC="isQC")
                      
                      
                      
)
)

setGeneric("normalize",function(para,method="sum",valueID="value",
                                norFactor=NULL,...) 
  standardGeneric("normalize"))
##' @describeIn normalize
setMethod("normalize",signature(para="metaXpara"),
          function(para,method="sum",valueID="value",norFactor=NULL,...){
            message("normalize: ",valueID)
            #para@peaksData->x
            #x<-dcast(x,ID~sample,value.var = valueID)
            x <- para@peaksData %>% 
              select_("sample",valueID,"ID") %>% 
              spread_("sample",valueID)
            dat <- as.matrix(x[,-1])
            
            if(!is.null(method) && method == "vsn"){
              fdata <- ExpressionSet(assayData=dat)
              fit<-vsn2(fdata)
              fitdata = predict(fit, newdata=fdata)  ## apply fit
              nd<-exprs(fitdata)
              x[,-1] <- 2^nd
              
            }else if(!is.null(method) && method == "quantiles"){
              x[,-1] <- preprocessCore::normalize.quantiles(dat)
              
            }else if(!is.null(method) && method == "quantiles.robust"){
              x[,-1] <- preprocessCore::normalize.quantiles.robust(dat)
              
            }else if(!is.null(method) && method == "sum"){
              x[,-1] <- apply(dat,2,function(x){x/sum(x,na.rm=TRUE)}) * mean(apply(dat,2,function(x){sum(x,na.rm=TRUE)}))
              
            }else if(!is.null(method) && method == "pqn"){
              ## Normalization of data with the Probabilistic Quotient 
              ## Normalization method.
              ## "pqn": the Probabilistic Quotient Normalization is computed as 
              ## described in Dieterle, et al. (2006).
              x[,-1] <- apply(dat,2,function(y){y/sum(y,na.rm=TRUE)})
              newX <- t(apply(x[,-1],1,function(y){y/median(y,na.rm=TRUE)}))
              coe <- apply(newX, 2, median, na.rm = TRUE)
              for(i in 2:ncol(x)){
                x[,i] <- x[,i]/coe[i-1]
              }
              x[,-1] <- x[,-1] * mean(apply(dat,2,function(x){sum(x,na.rm=TRUE)}))
              save(x,file="x.rda")
            }else if(!is.null(method) && method == "combat"){
              ## The ComBat function adjusts for known batches using an 
              ## empirical Bayesian framework
              ## http://biostatistics.oxfordjournals.org/content/8/1/118.full
              message("normalization: ComBat\n")
              sample_meta <- data.frame(sample=names(x)[-1],
                                        stringsAsFactors = FALSE)
              samList <- read.delim(para@sampleListFile,stringsAsFactors=FALSE)
              m <- left_join(sample_meta,samList,by="sample")
              dat_batch <- m$batch
              modcombat = model.matrix(~1, data=m)
              x_tmp = ComBat(dat=dat, batch=dat_batch, mod=modcombat, 
                             par.prior=TRUE, 
                             prior.plot=FALSE)
              message("<=0:",sum(x_tmp<=0),"/",length(x_tmp),"\n")
              x_tmp[x_tmp<=0] <- NA
              x[,-1] <- x_tmp
              
            }else if(is.null(method) && !is.null(norFactor)){
              message("Normalization by ",norFactor)
              if(norFactor %in% names(para@peaksData)){
                para@peaksData[,valueID] <- para@peaksData[,valueID]/para@peaksData[,norFactor]
              }else{
                stop("Not found ",norFactor," in para@peaksData!")
              }
            }
            
            #else{
            #   stop("Normalization method must be in ",
            #        "'vsn,quantiles,quantiles.robust,sum'\n")
            #}
            
            if(!is.null(method) && method != "none"){
              normValue <- melt(x,id.vars="ID",value.name=valueID,
                                variable.name="sample")
              para@peaksData[,valueID] <- NULL
              para@peaksData <- plyr::join(para@peaksData,normValue,
                                           by = c("ID","sample"))
            }
            return(para)
            
          }
          
)

##' @title rawPeaks
##' @description rawPeaks
##' @rdname rawPeaks
##' @param para An object of metaXpara
##' @param value value
##' @return An object of metaXpara
##' @docType methods
##' @exportMethod 
##' @author Bo Wen \email{wenbo@@genomics.cn}
##' @examples
##' library(faahKO)
##' para <- new("metaXpara")
##' rawPeaks(para) <- data.frame()
setGeneric("rawPeaks<-",function(para,value) standardGeneric("rawPeaks<-"))
##' @describeIn metaXpara
setReplaceMethod("rawPeaks", signature(para = "metaXpara"), 
                 function(para,value){
                   para@rawPeaks <- value
                   para
                 }
)

##' @title sampleListFile
##' @description sampleListFile
##' @rdname sampleListFile
##' @param para An object of metaXpara
##' @param value value
##' @return An object of metaXpara
##' @docType methods
##' @exportMethod 
##' @author Bo Wen \email{wenbo@@genomics.cn}
##' @examples
##' para <- new("metaXpara")
##' sampleListFile(para) <- "sample.txt"
setGeneric("sampleListFile<-",function(para,value) 
  standardGeneric("sampleListFile<-"))
##' @describeIn metaXpara
setReplaceMethod("sampleListFile", signature(para = "metaXpara"), 
                 function(para,value){
                   para@sampleListFile <- value
                   para
                 }
)

##' @title reSetPeaksData
##' @description reSetPeaksData
##' @rdname reSetPeaksData
##' @docType methods
##' @param para An object of \code{metaXpara}
##' @return none
##' @exportMethod
##' @author Bo Wen \email{wenbo@@genomics.cn}
##' @examples
##' para <- new("metaXpara")
##' pfile <- system.file("extdata/MTBLS79.txt",package = "metaX")
##' sfile <- system.file("extdata/MTBLS79_sampleList.txt",package = "metaX")
##' rawPeaks(para) <- read.delim(pfile,check.names = FALSE)
##' sampleListFile(para) <- sfile
##' para <- reSetPeaksData(para)
setGeneric("reSetPeaksData",function(para) standardGeneric("reSetPeaksData"))
##' @describeIn reSetPeaksData
setMethod("reSetPeaksData", signature(para = "metaXpara"), function(para){
  
  rawPeaks <- para@rawPeaks
  samList  <- read.delim(para@sampleListFile)
  row.names(rawPeaks) <- rawPeaks$name
  rawPeaks <- rawPeaks[,names(rawPeaks) %in% samList$sample]  
  rawPeaks$ID <- row.names(rawPeaks)
  peaksData <- melt(rawPeaks,id.vars = "ID",variable.name = "sample")
  peaksData <- plyr::join(peaksData,samList,by="sample")
  para@peaksData <- peaksData
  
  return(para)
  
})

##' @title Get a data.frame which contained the peaksData in metaXpara
##' @description Get a data.frame which contained the peaksData in metaXpara
##' @rdname getPeaksTable
##' @docType methods
##' @param para An object of data
##' @param sample Sample class used
##' @param valueID The column name used
##' @return A data.frame
##' @author Bo Wen \email{wenbo@@genomics.cn}
##' @examples
##' para <- new("metaXpara")
##' pfile <- system.file("extdata/MTBLS79.txt",package = "metaX")
##' sfile <- system.file("extdata/MTBLS79_sampleList.txt",package = "metaX")
##' rawPeaks(para) <- read.delim(pfile,check.names = FALSE)
##' sampleListFile(para) <- sfile
##' para <- reSetPeaksData(para)
##' res <- getPeaksTable(para)
setGeneric("getPeaksTable",function(para,sample=NULL,valueID="value") 
  standardGeneric("getPeaksTable"))
##' @describeIn getPeaksTable
setMethod("getPeaksTable", signature(para = "metaXpara"), 
          function(para,sample=NULL,valueID="value"){
            
            message("value = ",valueID)
            
            #sampleList  <- read.delim(para@sampleListFile)
            
            peaksData <- para@peaksData
            if(!is.null(sample)){
              peaksData <- peaksData[peaksData$class %in% sample,]
            }
            peaksData$class <- as.character(peaksData$class)
            pData <- dcast(peaksData,sample+class+batch+order~ID,value.var = valueID)
            pData$class[is.na(pData$class)] <- "QC"
            
            return(pData)
            
          }
)


##' @title Batch correction using SVR normalization
##' @description Batch correction using SVR normalization
##' @param para A metaXpara object
##' @param ntop Default is top 5 correlated peaks
##' @param impute Whether do the imputation before and after normalization
##' @param cpu The number of cpu used, default is 0 which means using all cpu.
##' @author Bo Wen \email{wenbo@@genomics.cn}
##' @return A new object of metaXpara
##' @references X. Shen, X. Gong, Y. Cai, Y. Guo, J. Tu, H. Li, T.Zhang, 
##' J. Wang, F. Xue, and Z.-J. Zhu, Normalization and 
##' Integration of Large-Scale Metabolomics Data Using Support Vector 
##' Regression, Metabolomics, 2016, 12: 89
##' @export
svrNormalize=function(para,ntop=5,impute=TRUE,cpu=0){
  res <- list()           
  
  message("Using cpu: ",cpu)
  peaksData <- para@peaksData
  if(sum(is.na(peaksData$class)) <=0 ){
    message("Don't find QC sample in your data!")
    return(para)
  }
  
  message("The number of NA value in peaksData before SVR normalization: ",
          sum(is.na(peaksData$value) | peaksData$value <= 0))
  
  qcData <- peaksData[is.na(peaksData$class),]
  sampleData <- peaksData[!is.na(peaksData$class),]
  samList <- read.delim(para@sampleListFile,stringsAsFactors=FALSE)
  
  assign(x="maxOrder",value=max(samList[,para@sampleListHead['order']]),
         envir=.GlobalEnv)
  
  # require(doSMP) # will no longer work...
  cpu = ifelse(cpu==0,detectCores(),cpu)
  
  cl <- makeCluster(getOption("cl.cores", cpu))
  
  clusterExport(cl, c("svrfit","svrfit2"),envir=environment())
  
  qcData$ID_batch <- paste(qcData$ID,qcData$batch,sep="_")
  #qcData$ID_batch <- qcData$ID
  
  message(date(),"\tSVR modeling ...")
  
  if(ntop!=1){
    ## fit intensity with top n correlation peaks
    dmat <- qcData %>%
      select_("ID","sample","value") %>% 
      spread_("ID","value") %>% select(-sample)
    
    corQC<-pcor(as.matrix(dmat),use="complete.obs")
    
    ## generate 
    dsData <- list(dmat=list(),smat=list())
    for(batchID in unique(peaksData$batch)){
      
      dsData$dmat[[batchID]] <- filter(qcData,batch==batchID) %>%
        select_("ID","sample","value") %>% 
        spread_("ID","value") %>% select(-sample)
      
      dsData$smat[[batchID]] <- filter(peaksData,batch==batchID) %>%
        select_("ID","sample","value") %>% 
        spread_("ID","value")
    }
    gc()
    #intPredict <- parLapply(cl,unique(qcData$ID_batch),svrfit3,corQC=corQC,
    #                        qcData=qcData,ntop=5,peaksData=peaksData)
    intPredict <- parLapply(cl,unique(qcData$ID_batch),svrfit3,corQC=corQC,
                            ntop=5,dsData=dsData,qcData=qcData)
    intPredict <- rbindlist(intPredict)
    intPredict <- as.data.frame(intPredict)
    
  }else{
    ## fit intensity with injection order
    intPredict <- parLapply(cl,unique(qcData$ID_batch),svrfit,
                            qcData=qcData,maxOrder=maxOrder)
    intPredict <- rbindlist(intPredict)
    intPredict <- as.data.frame(intPredict)
    
  }
  
  stopCluster(cl)
  message(date(),"\tSVR modeling done.")
  if("newOrder" %in% names(intPredict)){
    intPredict <- dplyr::rename(intPredict,order=newOrder)
  }
  ## head: sample       ID     value batch class order valuePredict
  peaksData$valuePredict <- NULL
  peaksData$valueNorm <- NULL
  peaksData <- plyr::join(peaksData,intPredict,
                          by=intersect(names(peaksData),names(intPredict)))
  
  mpa <- ddply(peaksData,.(ID),summarise,mpa=median(value,na.rm = TRUE))
  peaksData <- plyr::join(peaksData,mpa,
                          by=intersect(names(peaksData),names(mpa)))
  peaksData <- dplyr::mutate(peaksData,valuePredict=valuePredict/mpa)
  peaksData$mpa <- NULL
  message("change value which ==0 to NA")
  peaksData$value[ peaksData$value<=0] <- NA
  
  peaksData$valueNorm <- peaksData$value/peaksData$valuePredict
  
  ######################################################################
  ######################################################################
  peaksData$valuePredict[ peaksData$valuePredict<=0] <- NA
  
  ## calculate CV using this value
  peaksData$valueNorm[ peaksData$valueNorm<=0] <- NA
  message("The number of NA value in peaksData after SVR normalization:",
          sum(is.na(peaksData$valueNorm)))
  ## TODO: Do we still need to do missing value imputation are QC-based 
  ## correction?
  if(impute){
    message(date(),"\tDo missing value imputation after SVR normalization...")
    paraTmp <- para
    paraTmp@peaksData <- peaksData
    paraTmp <- missingValueImpute(x = paraTmp,valueID="valueNorm")
    peaksData <- paraTmp@peaksData
  }
  ## For each batch
  ## CV plot
  cvStat <- ddply(peaksData[is.na(peaksData$class),],.(ID,batch),
                  summarise,
                  rawCV=sd(value,na.rm = TRUE)/mean(value,na.rm = TRUE),
                  normCV=sd(valueNorm,na.rm = TRUE)/mean(valueNorm,na.rm = TRUE))
  
  
  cvStatForEachBatch <- melt(cvStat,id.vars = c("ID","batch"),
                             variable.name = "CV")
  cvStatForEachBatch$batch <- as.factor(cvStatForEachBatch$batch)
  
  ## output information
  message("Summary information of the CV for QC samples:")
  cvTable <- ddply(cvStatForEachBatch,.(batch,CV),summarise,
                   lessThan30=sum(value<=0.3,na.rm = TRUE),
                   total=length(value),ratio=lessThan30/total)
  print(cvTable)
  res$cvBatch <- cvTable
  message("\n")
  para@fig$cv<- paste(para@outdir,"/",para@prefix,"-cv.pdf",sep="") 
  pdf(para@fig$cv,width = 6,height = 6)
  p<-ggplot(data=cvStatForEachBatch,aes(x=value,fill=CV,colour=CV))+
    facet_grid(batch~.)+
    geom_density(alpha = 0.5)+
    xlab(label = "CV")
  print(p)
  p<-ggplot(data=cvStatForEachBatch,aes(x=value,fill=CV,colour=CV))+
    facet_grid(batch~.)+
    geom_density(alpha = 0.5)+
    xlim(0,2)+
    xlab(label = "CV")
  print(p)
  dev.off()
  
  #######################################
  ## For all
  cvStat <- ddply(peaksData[is.na(peaksData$class),],.(ID),
                  summarise,
                  rawCV=sd(value,na.rm = TRUE)/mean(value,na.rm = TRUE),
                  normCV=sd(valueNorm,na.rm = TRUE)/mean(valueNorm,na.rm = TRUE))
  
  
  cvStatForAll <- melt(cvStat,id.vars = c("ID"),
                       variable.name = "CV")
  ## output information
  message("Summary information of the CV for QC samples:")
  cvTable <- ddply(cvStatForAll,.(CV),summarise,
                   lessThan30=sum(value<=0.3,na.rm = TRUE),
                   total=length(value),ratio=lessThan30/total)
  print(cvTable)
  res$cvAll <- cvTable
  message("\n")
  
  ########################################################
  message("Peaks with CV > ",0.3,"!")
  message(sum(cvStat$normCV > 0.3,na.rm = TRUE))
  tmpPeaksData <- merge(peaksData,cvStat,by="ID")
  if(nrow(tmpPeaksData)!=nrow(peaksData)){
    error_file <- paste(para@outdir,"/",para@prefix,"-doSVR-error.rda",
                        sep="")
    message("Please see the file: ",error_file," for detail!")
    save(peaksData,cvStat,file=error_file)
    stop("Please see detailed data in ",error_file)
  }
  peaksData <- tmpPeaksData
  peaksData <- dplyr::rename(peaksData,cv=normCV)
  para@peaksData <- peaksData
  #para@peaksData <- dplyr::filter(peaksData,!is.na(cv) & cv<=cvFilter)
  res$metaXpara <- para
  return(res)
}


