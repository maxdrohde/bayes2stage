#' Calculate V_i = Z_i D t(Z_i) + sig_e^2 `I_(n_i)`
#'
#' Calculate V_i = Z_i D t(Z_i) + sig_e^2 `I_(n_i)`
#' @param zi n_i by q design matrix for the random effects
#' @param sigma.vc vector of variance components on standard deviation scale
#' @param rho.vc vector of correlations among the random effects.  The length should be q choose 2
#' @param sigma.e std dev of the measurement error distribution
#' @return V_i
#' @export
#'
vi.calc <- function(zi, sigma.vc, rho.vc, sigma.e){
    SDMat.RE  <- diag(sigma.vc)
    ncolzi    <- ncol(zi) ## make sure this equals length(sigma.vc)
    nrowzi    <- nrow(zi)
    nERRsd    <- length(sigma.e)
    b         <- matrix(0,ncolzi,ncolzi)
    b[lower.tri(b, diag=FALSE)] <- rho.vc
    CorMat.RE <- t(b)+b+diag(rep(1,ncolzi))
    CovMat.RE <- SDMat.RE %*% CorMat.RE %*% SDMat.RE
    zi %*% CovMat.RE %*% t(zi) + diag(rep(sigma.e^2, each=nrowzi/nERRsd))
}
#' Ascertainment correction piece for univariate sampling
#'
#' Calculate the (not yet log transformed) ascertainment correction under a univariate Q_i
#'
#' @param cutpoints cutpoints defining the sampling regions. (a vector of length 2)
#' @param SampProb Sampling probabilities from within each region (vector of length 3).
#' @param mu_q a scalar for the mean value of the Q_i distribution
#' @param sigma_q a scalar for the standard deviation of the Q_i distribution
#' @return Not yet log transformed ascertainment correction
#' @export
ACi1q <- function(cutpoints, SampProb, mu_q, sigma_q){
    CDFs <- pnorm(c(-Inf, cutpoints, Inf), mu_q, sigma_q)
    sum( SampProb*(CDFs[2:length(CDFs)] - CDFs[1:(length(CDFs)-1)]) )
}

#' Ascertainment correction piece for bivariate sampling
#'
#' Calculate the (not yet log transformed) ascertainment correction under a bivariate Q_i
#'
#' @param cutpoints cutpoints defining the sampling regions. (a vector of length 4: c(xlow, xhigh, ylow, yhigh))
#' @param SampProb Sampling probabilities from within each of two sampling regions; central region and outlying region (vector of length 2).
#' @param mu_q a 2-vector for the mean value of the bivariate Q_i distribution.
#' @param sigma_q a 2 by 2 covariance matrix for the bivariate Q_i distribution.
#' @return Not yet log transformed ascertainment correction
#' @export
ACi2q <- function(cutpoints, SampProb, mu_q, sigma_q){
    (SampProb[1]-SampProb[2])*pmvnorm(lower=c(cutpoints[c(1,3)]), upper=c(cutpoints[c(2,4)]), mean=mu_q, sigma=sigma_q)[[1]] + SampProb[2]
}

#' Log of the Ascertainment correction for univariate sampling
#'
#' Calculate the log transformed ascertainment correction under a univariate Q_i.  Also return vi
#'
#' @param yi n_i-response vector
#' @param xi n_i by p design matrix for fixed effects
#' @param zi n_i by q design matric for random effects (intercept and slope)
#' @param wi the pre-multiplier of yi to generate the sampling variable q_i
#' @param beta mean model parameter vector
#' @param sigma.vc vector of variance components on standard deviation scale
#' @param rho.vc vector of correlations among the random effects.  The length should be q choose 2
#' @param sigma.e std dev of the measurement error distribution
#' @param cutpoints cutpoints defining the sampling regions. (a vector of length 2)
#' @param SampProb Sampling probabilities from within each region (vector of length 3).
#' @return log transformed ascertainment correction
#' @export
logACi1q <- function(yi, xi, zi, wi, beta, sigma.vc, rho.vc, sigma.e, cutpoints, SampProb){
    vi      <- vi.calc(zi, sigma.vc, rho.vc, sigma.e)
    mu      <- xi %*% beta
    mu_q    <- (wi %*% mu)[,1]
    sigma_q <- sqrt((wi %*% vi %*% t(wi))[1,1])
    return(list(vi=vi, logACi=log(ACi1q(cutpoints, SampProb, mu_q, sigma_q))))
}


#' Log of the Ascertainment correction piece for bivariate sampling
#'
#' Calculate the log transformed ascertainment correction under a bivariate Q_i.  Also return vi
#'
#' @param yi n_i-response vector
#' @param xi n_i by p design matrix for fixed effects
#' @param zi n_i by q design matric for random effects (intercept and slope)
#' @param wi the pre-multiplier of yi to generate the sampling variable q_i
#' @param beta mean model parameter p-vector
#' @param sigma.vc vector of variance components on standard deviation scale
#' @param rho.vc vector of correlations among the random effects.  The length should be q choose 2
#' @param sigma.e std dev of the measurement error distribution
#' @param cutpoints cutpoints defining the sampling regions. (a vector of length 4 c(xlow, xhigh, ylow, yhigh))
#' @param SampProb Sampling probabilities from within each region (vector of length 2 c(central region, outlying region)).
#' @return log transformed ascertainment correction
#' @export
logACi2q <- function(yi, xi, zi, wi, beta, sigma.vc, rho.vc, sigma.e, cutpoints, SampProb){
    vi      <- vi.calc(zi, sigma.vc, rho.vc, sigma.e)
    mu      <- xi %*% beta
    mu_q    <- as.vector(wi %*% mu)
    sigma_q <- wi %*% vi %*% t(wi)
    ## for some reason the upper and lower triangles to not always equal, so I am taking their
    ## average.  Not sure this is a problem here
    ## but doing this to be safe.  Maybe can remove later once understood.
    sigma_q <- (sigma_q + t(sigma_q))/2
    #sigma_q[upper.tri(sigma_q)]  <- t(sigma_q)[upper.tri(sigma_q)]
    #sigma_q[2,1] <- sigma_q[1,2]
    return(list(vi=vi, logACi= log( ACi2q(cutpoints=cutpoints, SampProb=SampProb, mu_q=mu_q, sigma_q=sigma_q))))
}


#' Calculate a subject-specific contribution to a log-likelihood for longitudinal normal data
#'
#' Calculate a subject-specific contribution to a log-likelihood for longitudinal normal data
#' @param yi n_i-response vector
#' @param xi n_i by p design matrix for fixed effects
#' @param beta mean model parameter vector
#' @param vi the variance covariance matrix (ZDZ+Sige2*I)
#' @return subject specific contribution to the log-likelihood
#' @export
#'
li.lme <- function(yi, xi, beta, vi){
    resid <- yi - xi %*% beta
    -(1/2) * (length(xi[,1])*log(2*pi) + log(det(vi)) + t(resid) %*% solve(vi) %*% resid )[1,1]
}

#' Calculate the conditional likelihood for the univariate and bivariate sampling cases across all subjects (Keep.liC=FALSE) or the subject specific contributions to the conditional likelihood along with the log-transformed ascertainment correction for multiple imputation (Keep.liC=TRUE).
#'
#' Calculate the conditional likelihood for the univariate and bivariate sampling cases across all subjects (Keep.liC=FALSE) or the subject specific contributions to the conditional likelihood along with the log-transformed ascertainment correction for multiple imputation (Keep.liC=TRUE).
#'
#' @param y response vector
#' @param x sum(n_i) by p design matrix for fixed effects
#' @param z sum(n_i) by q design matrix for random effects
#' @param w.function sum(n_i) vector with possible values that include "mean" "intercept" "slope" and "bivar."  There should be one unique value per subject
#' @param id sum(n_i) vector of subject ids
#' @param beta mean model parameter p-vector
#' @param sigma.vc vector of variance components on standard deviation scale
#' @param rho.vc vector of correlations among the random effects.  The length should be q choose 2
#' @param sigma.e std dev of the measurement error distribution
#' @param cutpoints A matrix with the first dimension equal to sum(n_i).  These cutpoints define the sampling regions (bivariate Q_i: each row is a vector of length 4 c(xlow, xhigh, ylow, yhigh); univariate Q_i: each row is a vector of length 2 c(k1,k2) to define the sampling regions, i.e., low, middle, high).  Each subject should have n_i rows of the same values.
#' @param SampProb A matrix with the first dimension equal to sum(n_i).   Sampling probabilities from within each region (bivariate Q_i: each row is a vector of length 2 c(central region, outlying region); univariate Q_i: each row is a vector of length 3 with sampling probabilities for each region). Each subject should have n_i rows of the same values.
#' @param Weights Subject specific sampling weights.  A vector of length sum(n_i).  Not used unless using weighted Likelihood
#' @param Keep.liC If FALSE, the function returns the conditional log likelihood across all subjects.  If TRUE, subject specific contributions and exponentiated subject specific ascertainment corrections are returned in a list.
#' @return If Keep.liC=FALSE, conditional log likelihood.  If Keep.liC=TRUE, a two-element list that contains subject specific likelihood contributions and exponentiated ascertainment corrections.
#' @export
#'
LogLikeC2 <- function(y, x, z, w.function, id, beta, sigma.vc, rho.vc, sigma.e, cutpoints, SampProb, Weights, Keep.liC=FALSE){

    subjectData    <- CreateSubjectData(id=id,y=y,x=x,z=z,Weights=Weights,SampProb=SampProb,cutpoints=cutpoints,w.function=w.function)
    liC.and.logACi <- lapply(subjectData, LogLikeiC2, beta=beta, sigma.vc=sigma.vc, rho.vc=rho.vc, sigma.e=sigma.e)

    if (Keep.liC == FALSE){out <- -1*Reduce('+', liC.and.logACi)[1]  ## sum ss contributions to liC
    }else{ out <- list(liC    = c(unlist(sapply(liC.and.logACi, function(x) x[1]))), ## ss contributions to liC
                       logACi = c(unlist(sapply(liC.and.logACi, function(x) x[2]))))} ## ss contributions to ACi
    out
}



#' Calculate the ss contributions to the conditional likelihood for the univariate and bivariate sampling cases.
#'
#' Calculate the ss contributions to the conditional likelihood for the univariate and bivariate sampling cases.
#'
#' @param subjectData a list containing: yi, xi, zi, Weights.i, w.function.i, SampProb.i, cutpoints.i
#' @param beta mean model parameter p-vector
#' @param sigma.vc vector of variance components on standard deviation scale
#' @param rho.vc vector of correlations among the random effects.  The length should be q choose 2
#' @param sigma.e std dev of the measurement error distribution
#' @return ss contributions to the conditional log likelihood.  This is an internal function used by LogLikeC2
#' @export
#'
#'
LogLikeiC2 <- function(subjectData, beta, sigma.vc, rho.vc, sigma.e){
    yi          <- subjectData[["yi"]]
    xi          <- subjectData[["xi"]]
    zi          <- subjectData[["zi"]]
    Weights.i   <- subjectData[["Weights.i"]]
    w.function  <- subjectData[["w.function.i"]]
    SampProb    <- subjectData[["SampProb.i"]]
    cutpoints   <- subjectData[["cutpoints.i"]]
    ni          <- length(yi)
    t.zi        <- t(zi)
##########
##########
##########
##########
    wi.tmp <- solve(t.zi %*% zi) %*% t.zi
    if (!(w.function %in% c("bivar", "mvints", "mvslps"))){
        if (w.function %in% c("intercept", "intercept1")){ wi <- wi.tmp[1,]
        } else if (w.function %in% c("slope", "slope1")){  wi <- wi.tmp[2,]
        } else if (w.function %in% c("intercept2")){       wi <- wi.tmp[3,]
        } else if (w.function %in% c("slope2")){           wi <- wi.tmp[4,]
        } else if (w.function=="mean"){                    wi <- t(rep(1/ni, ni))
        }
        wi         <- matrix(wi, 1, ni)
        tmp        <- logACi1q(yi, xi, zi, wi, beta, sigma.vc, rho.vc, sigma.e, cutpoints, SampProb)
        logACi     <- tmp[["logACi"]]
        liC        <- li.lme(yi, xi, beta, tmp[["vi"]])*Weights.i - logACi
    }else {
        if (w.function %in% c("bivar")){ wi <- wi.tmp[c(1,2),]
        } else if (w.function %in% c("mvints")){ wi <- wi.tmp[c(1,3),]
        } else if (w.function %in% c("mvslps")){ wi <- wi.tmp[c(2,4),]
        }
    tmp        <- logACi2q(yi, xi, zi, wi, beta, sigma.vc, rho.vc, sigma.e, cutpoints, SampProb)
    logACi     <- tmp[["logACi"]]
    liC        <- li.lme(yi, xi, beta, tmp[["vi"]])*Weights.i - logACi
    }
##########
##########
##########
##########
    # if (w.function != "bivar"){
    #     if (w.function %in% c("intercept", "intercept1")){ wi<- (solve(t.zi %*% zi) %*% t.zi)`[1,]`
    #     } else if (w.function %in% c("slope", "slope1")){     wi<- (solve(t.zi %*% zi) %*% t.zi)`[2,]`
    #     } else if (w.function %in% c("intercept2")){ wi<- (solve(t.zi %*% zi) %*% t.zi)`[3,]`
    #     } else if (w.function %in% c("slope2")){     wi<- (solve(t.zi %*% zi) %*% t.zi)`[4,]`
    #     } else if (w.function=="mean"){     wi <- t(rep(1/ni, ni))
    #     }
    #     wi         <- matrix(wi, 1, ni)
    #     tmp        <- logACi1q(yi, xi, zi, wi, beta, sigma.vc, rho.vc, sigma.e, cutpoints, SampProb)
    #     logACi     <- tmp[["logACi"]]
    #     liC        <- li.lme(yi, xi, beta, tmp[["vi"]])*Weights.i - logACi
    #
    # }else{
    #     wi         <- solve(t.zi %*% zi) %*% t.zi
    #     tmp        <- logACi2q(yi, xi, zi, wi, beta, sigma.vc, rho.vc, sigma.e, cutpoints, SampProb)
    #     logACi     <- tmp[["logACi"]]
    #     liC        <- li.lme(yi, xi, beta, tmp[["vi"]])*Weights.i - logACi
    #
    # }
    return(c(liC, logACi))
}


#' Gradient of the log of the ascertainment correction piece for sampling based on bivariate Q_i.
#'
#' Calculate the gradient of the log transformed ascertainment correction under designs that sample based on a bivariate Q_i (numerically).
#'
#' @param subjectData a list containing: yi, xi, zi, Weights.i, w.function, SampProb, cutpoints
#' @param beta mean model parameter p-vector
#' @param sigma.vc vector of variance components on standard deviation scale
#' @param rho.vc vector of correlations among the random effects.  The length should be q choose 2
#' @param sigma.e std dev of the measurement error distribution
#' @return gradient of the log transformed ascertainment correction under the bivariate sampling design
#' @importFrom numDeriv grad
#' @importFrom mvtnorm pmvnorm
#' @export
logACi2q.score2 <- function(subjectData, beta, sigma.vc, rho.vc, sigma.e){

    yi          <- subjectData[["yi"]]
    xi          <- subjectData[["xi"]]
    zi          <- subjectData[["zi"]]
    Weights.i   <- subjectData[["Weights.i"]]  ## not used
    w.function  <- subjectData[["w.function.i"]]
    SampProb    <- subjectData[["SampProb.i"]]
    cutpoints   <- subjectData[["cutpoints.i"]]

    t.zi        <- t(zi)
    #wi          <- solve(t.zi %*% zi) %*% t.zi
    ###########################
    ###########################
    ###########################
    ###########################
    wi.tmp      <- solve(t.zi %*% zi) %*% t.zi
    if (w.function %in% c("bivar")){         wi <- wi.tmp[c(1,2),]
    } else if (w.function %in% c("mvints")){ wi <- wi.tmp[c(1,3),]
    } else if (w.function %in% c("mvslps")){ wi <- wi.tmp[c(2,4),]
    }
    ###########################
    ###########################
    ###########################
    ###########################
    t.wi        <- t(wi)

    param   <- c(beta, sigma.vc, rho.vc, sigma.e)
    npar    <- length(param)

    len.beta     <- length(beta)
    len.sigma.vc <- length(sigma.vc)
    len.rho.vc   <- length(rho.vc)
    len.sigma.e  <- length(sigma.e)

    beta.index   <- c(1:len.beta)
    vc.sd.index  <- len.beta + (c(1:len.sigma.vc))
    vc.rho.index <-  len.beta + len.sigma.vc + (c(1:len.rho.vc))
    err.sd.index <- len.beta + len.sigma.vc + len.rho.vc + c(1:len.sigma.e)

    Deriv   <- sapply(1:npar,  function(rr)
    {
        grad(function(x) { new.param <- param
                           new.param[rr] <- x
                           vi      <- vi.calc(zi, new.param[vc.sd.index], new.param[vc.rho.index], new.param[err.sd.index])
                           mu_q    <- as.vector(wi %*% (xi %*% new.param[c(beta.index)]))
                           sigma_q <- wi %*% vi %*% t.wi
                           ## for some reason the upper and lower triangles to not always equal.  Not sure this is a problem here
                           ## but doing this to be safe.  Maybe can remove later once understood.
                           sigma_q <- (sigma_q + t(sigma_q))/2
                           pmvnorm(lower=c(cutpoints[c(1,3)]), upper=c(cutpoints[c(2,4)]), mean=mu_q, sigma=sigma_q)[[1]]
        },
        param[rr],
        method="simple",
        method.args=list(eps=1e-4))
    }
    )

    vi      <- vi.calc(zi, sigma.vc, rho.vc, sigma.e)
    mu_q    <- as.vector(wi %*% (xi %*% beta))
    sigma_q <- wi %*% vi %*% t.wi
    ## for some reason the upper and lower triangles to not always equal, so I am forcing
    ## the upper triangle to equal the lower triangle.
    sigma_q <- (sigma_q + t(sigma_q))/2

    (SampProb[1]-SampProb[2])*Deriv / ACi2q(cutpoints, SampProb, mu_q, sigma_q)
}



#' Gradient of the log of the ascertainment correction piece for sampling based on univariate Q_i
#'
#' Calculate the gradient of the log transformed ascertainment correction for sampling based on univariate Q_i
#'
#' @param subjectData a list containing: yi, xi, zi, Weights.i, w.function.i, SampProb.i, cutpoints.i
#' @param beta mean model parameter p-vector
#' @param sigma.vc vector of variance components on standard deviation scale
#' @param rho.vc vector of correlations among the random effects.  The length should be q choose 2
#' @param sigma.e std dev of the measurement error distribution
#' @return gradient of the log transformed ascertainment correction under univariate $Q_i$
#' @export
#' @importFrom stats dnorm
logACi1q.score2 <- function(subjectData, beta, sigma.vc, rho.vc, sigma.e){

    yi          <- subjectData[["yi"]]
    xi          <- subjectData[["xi"]]
    zi          <- subjectData[["zi"]]
    Weights.i   <- subjectData[["Weights.i"]]  ## not used
    w.function  <- subjectData[["w.function.i"]]
    SampProb    <- subjectData[["SampProb.i"]]
    cutpoints   <- subjectData[["cutpoints.i"]]

    t.zi        <- t(zi)
    ni          <- length(yi)
    #if (w.function=="mean")      wi <- t(rep(1/ni, ni))
    #if (w.function=="intercept") wi<- (solve(t.zi %*% zi) %*% t.zi)`[1,]`
    #if (w.function=="slope")     wi<- (solve(t.zi %*% zi) %*% t.zi)`[2,]`
    if (w.function %in% c("intercept", "intercept1")){ wi<- (solve(t.zi %*% zi) %*% t.zi)[1,]
    } else if (w.function %in% c("slope", "slope1")){     wi<- (solve(t.zi %*% zi) %*% t.zi)[2,]
    } else if (w.function %in% c("intercept2")){ wi<- (solve(t.zi %*% zi) %*% t.zi)[3,]
    } else if (w.function %in% c("slope2")){     wi<- (solve(t.zi %*% zi) %*% t.zi)[4,]
    } else if (w.function=="mean"){     wi <- t(rep(1/ni, ni))
    }
    wi      <- matrix(wi, 1, ni)

    vi        <- vi.calc(zi, sigma.vc, rho.vc, sigma.e)
    t.wi      <- t(wi)
    wi.zi     <- wi %*% zi
    t.wi.zi   <- t(wi.zi)

    mu        <- xi %*% beta
    mu_q      <- (wi %*% mu)[,1]
    sigma_q   <- sqrt((wi %*% vi %*% t.wi)[1,1])

    l <- ACi1q(cutpoints, SampProb, mu_q, sigma_q)
    p <- SampProb[1:(length(SampProb)-1)] - SampProb[2:(length(SampProb))]
    f <- dnorm(cutpoints, mu_q, sigma_q)

    d_li_beta <- (wi %*% xi) * sum(p*f) / l
    f_alpha_k <- sum(p*f*(cutpoints - mu_q)) / (l * 2* sigma_q^2 )

    #################################################################################
    ### now calculate d_li_dalpha
    #################################################################################
    len.sigma.vc <- length(sigma.vc)
    len.rho.vc   <- length(rho.vc)
    len.sigma.e  <- length(sigma.e)

    SDMat.RE     <- diag(sigma.vc)

    b0         <- matrix(0, len.sigma.vc, len.sigma.vc)
    b0[lower.tri(b0, diag=FALSE)] <- rho.vc
    tbb0       <- t(b0) + b0
    CorMat.RE  <- tbb0 +  diag(rep(1,len.sigma.vc))

    D             <- SDMat.RE %*% CorMat.RE %*% SDMat.RE
    m             <- matrix(0,len.sigma.vc, len.sigma.vc)

    ## derivatives w.r.t variance components SDs
    dVi.dsigma.vc <- NULL
    for (mmm in 1:len.sigma.vc){
        m1               <- m
        m1[mmm,]         <- m1[,mmm] <- 1
        m1[mmm,mmm]      <- 2
        tmp              <- sigma.vc[mmm]+ 1*(sigma.vc[mmm]==0) ## to prevent unlikely division by 0
        dViMat.dsigma.vc <- m1 * D / tmp
        dVi.dsigma.vc    <- c(dVi.dsigma.vc, (wi.zi %*% dViMat.dsigma.vc %*% t.wi.zi)[1,1])
    }

    ## derivatives w.r.t variance components rhos
    b1         <- matrix(0,len.sigma.vc,len.sigma.vc)
    b1[lower.tri(b1, diag=FALSE)] <-  c(1:len.rho.vc)
    tbb1       <- t(b1) + b1

    dVi.drho.vc <- NULL
    for (mmm in 1:len.rho.vc){
        tmp                   <- which(tbb1==mmm, arr.ind=TRUE)
        m2                    <- m
        m2[tmp[1,1],tmp[1,2]] <- 1
        m2[tmp[2,1],tmp[2,2]] <- 1
        tmp                   <- rho.vc[mmm] + 1*(rho.vc[mmm]==0) ## to prevent division by 0
        dViMat.drho.vc        <- m2 * D / tmp
        dVi.drho.vc           <- c(dVi.drho.vc, (wi.zi %*% dViMat.drho.vc %*% t.wi.zi)[1,1])
    }

    ## derivatives w.r.t error sds
    dVi.dsigma.e <- NULL
    for (mmm in 1:len.sigma.e){
        dsigma.e.vec       <- 2*sigma.e
        dsigma.e.vec[-mmm] <- 0
        dViMat.dsigma.e    <- diag(rep(dsigma.e.vec, each=ni/len.sigma.e))
        dVi.dsigma.e       <- c(dVi.dsigma.e, (wi %*% dViMat.dsigma.e %*% t.wi)[1,1])
    }

    c(d_li_beta, c(f_alpha_k * c(dVi.dsigma.vc, dVi.drho.vc, dVi.dsigma.e)))
}


#' Subject specific contribution to the lme model score (also returns marginal Vi=Cov(Y|X))
#'
#' Subject specific contribution to the lme model score (also returns marginal Vi=Cov(Y|X))
#'
#' @param subjectData a list that contains yi, xi, zi, Weights.i.  Note that Weights.i is used for inverse probability weighting only.
#' @param beta mean model parameter p-vector
#' @param sigma.vc vector of variance components on standard deviation scale
#' @param rho.vc vector of correlations among the random effects.  The length should be q choose 2
#' @param sigma.e std dev of the measurement error distribution
#' @return Subject specific contribution to the log-likelihood score (also returns marginal Vi=Cov(Y|X))
#' @export
li.lme.score2 <- function(subjectData, beta, sigma.vc, rho.vc, sigma.e){
    yi          <- subjectData[["yi"]]
    xi          <- subjectData[["xi"]]
    zi          <- subjectData[["zi"]]
    t.zi        <- t(zi)
    ni          <- length(yi)
    Weights.i   <- subjectData[["Weights.i"]]

    resid         <- yi - xi %*% beta
    vi            <- vi.calc(zi, sigma.vc, rho.vc, sigma.e)
    inv.v         <- solve(vi)
    t.resid       <- t(resid)
    t.resid.inv.v <- t.resid %*% inv.v
    inv.v.resid   <- inv.v %*% resid

    dli.dbeta <- t(xi) %*% inv.v.resid # for Beta

    len.sigma.vc <- length(sigma.vc)
    len.rho.vc   <- length(rho.vc)
    len.sigma.e  <- length(sigma.e)

    SDMat.RE     <- diag(sigma.vc)

    b0         <- matrix(0, len.sigma.vc, len.sigma.vc)
    b0[lower.tri(b0, diag=FALSE)] <- rho.vc
    tbb0       <- t(b0) + b0
    CorMat.RE  <- tbb0 +  diag(rep(1,len.sigma.vc))

    D             <- SDMat.RE %*% CorMat.RE %*% SDMat.RE
    m             <- matrix(0,len.sigma.vc, len.sigma.vc)

    ## derivatives w.r.t variance components SDs
    dli.dsigma.vc <- NULL
    for (mmm in 1:len.sigma.vc){
        m1               <- m
        m1[mmm,]         <- m1[,mmm] <- 1
        m1[mmm,mmm]      <- 2
        tmp              <- sigma.vc[mmm]+ 1*(sigma.vc[mmm]==0) ## to prevent unlikely division by 0
        dViMat.dsigma.vc <- m1 * D / tmp
        tmp              <- zi %*% dViMat.dsigma.vc %*% t.zi
        dli.dsigma.vc    <- c(dli.dsigma.vc, -0.5*(sum(diag(inv.v %*% tmp)) - t.resid.inv.v %*% tmp %*% inv.v.resid)[1,1])
    }

    ## derivatives w.r.t variance components rhos
    b1         <- matrix(0,len.sigma.vc,len.sigma.vc)
    b1[lower.tri(b1, diag=FALSE)] <-  c(1:len.rho.vc)
    tbb1       <- t(b1) + b1

    dli.drho.vc <- NULL
    for (mmm in 1:len.rho.vc){
        tmp                   <- which(tbb1==mmm, arr.ind=TRUE)
        m2                    <- m
        m2[tmp[1,1],tmp[1,2]] <- 1
        m2[tmp[2,1],tmp[2,2]] <- 1
        tmp                   <- rho.vc[mmm] + 1*(rho.vc[mmm]==0) ## to prevent division by 0
        dViMat.drho.vc        <- m2 * D / tmp
        tmp                   <- zi %*% dViMat.drho.vc %*% t.zi
        dli.drho.vc           <- c(dli.drho.vc, -0.5*(sum(diag(inv.v %*% tmp)) - t.resid.inv.v %*% tmp %*% inv.v.resid)[1,1])
    }

    ## derivatives w.r.t error sds
    dli.dsigma.e <- NULL
    for (mmm in 1:len.sigma.e){
        dsigma.e.vec       <- 2*sigma.e
        dsigma.e.vec[-mmm] <- 0
        tmp                <- diag(rep(dsigma.e.vec, each=ni/len.sigma.e))
        dli.dsigma.e       <- c(dli.dsigma.e, -0.5*(sum(diag(inv.v %*% tmp)) - t.resid.inv.v %*% tmp %*% inv.v.resid)[1,1])
    }
    #print(c(dli.dsigma.vc, dli.drho.vc, dli.dsigma.e))
    list(gr=append(dli.dbeta, c(dli.dsigma.vc, dli.drho.vc, dli.dsigma.e))*Weights.i,
         vi=vi)
}


#' Calculate the gradient of the conditional likelihood for the univariate and bivariate sampling cases across all subjects (CheeseCalc=FALSE) or the cheese part of the sandwich estimator if CheeseCalc=TRUE.
#'
#' Calculate the gradient of the conditional likelihood for the univariate and bivariate sampling cases across all subjects (CheeseCalc=FALSE) or the cheese part of the sandwich estimator if CheeseCalc=TRUE.
#'
#' @param y response vector
#' @param x sum(n_i) by p design matrix for fixed effects
#' @param z sum(n_i) by 2 design matric for random effects (intercept and slope)
#' @param w.function sum(n_i) vector with possible values that include "mean" (mean of response series), "intercept" (intercept of the regression of Yi ~ zi where zi is the design matrix for the random effects (solve(t.zi %*% zi) %*% t.zi)`[1,]`), "intercept1"  (intercept of the regression of Yi ~ zi where zi is the design matrix for the random effects (solve(t.zi %*% zi) %*% t.zi)`[1,]`). "intercept2" (second intercept of the regression of the Yi ~ zi where zi is the design matrix for the bivariate random effects (b10,b11,b20,b21) solve(t.zi %*% zi) %*% t.zi)`[3,]`), "slope" (slope of the regression of Yi ~ zi where zi is the design matrix for the random effects (solve(t.zi %*% zi) %*% t.zi)`[2,]`), "slope1" (slope of the regression of Yi ~ zi where zi is the design matrix for the random effects (solve(t.zi %*% zi) %*% t.zi)`[2,]`), "slope2" (second slope of the regression of the Yi ~ zi where zi is the design matrix for the bivariate random effects (b10,b11,b20,b21) solve(t.zi %*% zi) %*% t.zi)`[4,]`) "bivar" (intercept and slope of the regression of Yi ~ zi where zi is the design matrix for the random effects (solve(t.zi %*% zi) %*% t.zi)`[c(1,2),]`) "mvints" (first and second intercepts of the bivariate regression of the Yi ~ zi where zi is the design matrix for the bivariate random effects (b10,b11,b20,b21) solve(t.zi %*% zi) %*% t.zi)`[c(1,3),]`) "mvslps" (first and second slopes of the bivariate regression of the Yi ~ zi where zi is the design matrix for the bivariate random effects (b10,b11,b20,b21) solve(t.zi %*% zi) %*% t.zi)`[c(1,3),]`).  There should be one unique value per subject
#' @param id sum(n_i) vector of subject ids
#' @param beta mean model parameter p-vector
#' @param sigma.vc vector of variance components on standard deviation scale
#' @param rho.vc vector of correlations among the random effects.  The length should be q choose 2
#' @param sigma.e std dev of the measurement error distribution
#' @param cutpoints A matrix with the first dimension equal to sum(n_i).  These cutpoints define the sampling regions (bivariate Q_i: each row is a vector of length 4 c(xlow, xhigh, ylow, yhigh); univariate Q_i: each row is a vector of length 2 c(k1,k2) to define the sampling regions, i.e., low, middle, high).  Each subject should have n_i rows of the same values.
#' @param SampProb A matrix with the first dimension equal to sum(n_i).   Sampling probabilities from within each region (bivariate Q_i: each row is a vector of length 2 c(central region, outlying region); univariate Q_i: each row is a vector of length 3 with sampling probabilities for each region). Each subject should have n_i rows of the same values.
#' @param Weights Subject specific sampling weights.  A vector of length sum(n_i).  Not used unless using weighted Likelihood
#' @param CheeseCalc If FALSE, the function returns the gradient of the conditional log likelihood across all subjects.  If TRUE, the cheese part of the sandwich esitmator is calculated.
#' @return If CheeseCalc=FALSE, gradient of conditional log likelihood.  If CheeseCalc=TRUE, the cheese part of the sandwich estimator is calculated.
#' @export
LogLikeC.Score2 <- function(y, x, z, w.function, id, beta, sigma.vc, rho.vc, sigma.e, cutpoints, SampProb, Weights, CheeseCalc=FALSE){
    param.vec <- c(beta, log(sigma.vc),log((1+rho.vc)/(1-rho.vc)),log(sigma.e))
    #print(c("blahblah", param.vec))
    npar     <- length(param.vec)

    len.beta     <- length(beta)
    len.sigma.vc <- length(sigma.vc)
    len.rho.vc   <- length(rho.vc)
    len.sigma.e  <- length(sigma.e)

    beta.index    <- c(1:len.beta)
    vc.sd.index   <- len.beta + (c(1:len.sigma.vc))
    vc.rho.index  <-  len.beta + len.sigma.vc + (c(1:len.rho.vc))
    err.sd.index  <- len.beta + len.sigma.vc + len.rho.vc + c(1:len.sigma.e)
    notbeta.index <- c(vc.sd.index,vc.rho.index,err.sd.index)

    subjectData <- CreateSubjectData(id=id,y=y,x=x,z=z,Weights=Weights,SampProb=SampProb,cutpoints=cutpoints,w.function=w.function)

    UncorrectedScorei <- lapply(subjectData, li.lme.score2, beta=beta, sigma.vc=sigma.vc, rho.vc=rho.vc, sigma.e=sigma.e)
    Gradienti         <- lapply(UncorrectedScorei, function(x) x[['gr']]) ## create a list of ss contributions to gradient
    UncorrectedScore  <- Reduce('+', Gradienti)  ## Note if using IPW this is actually a corrected score (corrected by the IPW)
    #print(c("blah1", UncorrectedScore))
    ## NOTE HERE: I used the first element of w.function in this call.  This means, for now, we cannot mix bivar with other
    ## sampling schemes.  This also applies to the cheese calculation
    if (!(w.function[[1]] %in% c("bivar","mvints","mvslps"))){
        logACi.Score <- lapply(subjectData, logACi1q.score2, beta=beta, sigma.vc=sigma.vc, rho.vc=rho.vc, sigma.e=sigma.e)
        logAC.Score  <- Reduce('+', logACi.Score)
        CorrectedScore <- UncorrectedScore + logAC.Score
    }else{
        logACi.Score <- lapply(subjectData, logACi2q.score2, beta=beta, sigma.vc=sigma.vc, rho.vc=rho.vc, sigma.e=sigma.e)
        logAC.Score  <- Reduce('+', logACi.Score)
        CorrectedScore <- UncorrectedScore - logAC.Score  ## notice this has the opposite sign compared to above.  Remember to check
    }
    #print(c("blah2", CorrectedScore))
    if (CheeseCalc==TRUE){
        if (!(w.function[[1]] %in% c("bivar","mvints","mvslps"))){ GradiMat <- mapply("+", Gradienti, logACi.Score)
        }else{                           GradiMat <- mapply("-", Gradienti, logACi.Score)} ## notice this has the opposite sign compared to above.  Remember to check
        ## Need to use the chain rule: note that param,vec is on the unconstrained scale but Gradi was calculated on the constrained parameters
        GradiMat[notbeta.index,] <- GradiMat[notbeta.index,]*c(exp(param.vec[vc.sd.index]), 2*exp(param.vec[vc.rho.index])/((exp(param.vec[vc.rho.index])+1)^2), exp(param.vec[err.sd.index]))
        cheese <- matrix(0,  npar, npar)
        for (mm in 1:ncol(GradiMat)) cheese <- cheese + outer(GradiMat[,mm], GradiMat[,mm])
    }
    out <- -CorrectedScore
    if (CheeseCalc==TRUE) out <- cheese
    out
}

#' Calculate the ascertainment corrected log likelihood and score
#'
#' Calculate the ascertainment corrected log likelihood and score
#'
#' @param params parameter vector c(beta, log(sigma0), log(sigma1), rho, sigmae)
#' @param y response vector
#' @param x sum(n_i) by p design matrix for fixed effects
#' @param z sum(n_i) by 2 design matric for random effects (intercept and slope)
#' @param id sum(n_i) vector of subject ids
#' @param w.function sum(n_i) vector with possible values that include "mean" (mean of response series), "intercept" (intercept of the regression of Yi ~ zi where zi is the design matrix for the random effects (solve(t.zi %*% zi) %*% t.zi)`[1,]`), "intercept1"  (intercept of the regression of Yi ~ zi where zi is the design matrix for the random effects (solve(t.zi %*% zi) %*% t.zi)`[1,]`). "intercept2" (second intercept of the regression of the Yi ~ zi where zi is the design matrix for the bivariate random effects (b10,b11,b20,b21) solve(t.zi %*% zi) %*% t.zi)`[3,]`), "slope" (slope of the regression of Yi ~ zi where zi is the design matrix for the random effects (solve(t.zi %*% zi) %*% t.zi)`[2,]`), "slope1" (slope of the regression of Yi ~ zi where zi is the design matrix for the random effects (solve(t.zi %*% zi) %*% t.zi)`[2,]`), "slope2" (second slope of the regression of the Yi ~ zi where zi is the design matrix for the bivariate random effects (b10,b11,b20,b21) solve(t.zi %*% zi) %*% t.zi)`[4,]`) "bivar" (intercept and slope of the regression of Yi ~ zi where zi is the design matrix for the random effects (solve(t.zi %*% zi) %*% t.zi)`[c(1,2),]`) "mvints" (first and second intercepts of the bivariate regression of the Yi ~ zi where zi is the design matrix for the bivariate random effects (b10,b11,b20,b21) solve(t.zi %*% zi) %*% t.zi)`[c(1,3),]`) "mvslps" (first and second slopes of the bivariate regression of the Yi ~ zi where zi is the design matrix for the bivariate random effects (b10,b11,b20,b21) solve(t.zi %*% zi) %*% t.zi)`[c(1,3),]`).  There should be one unique value per subject
#' @param cutpoints A matrix with the first dimension equal to sum(n_i).  These cutpoints define the sampling regions (bivariate Q_i: each row is a vector of length 4 c(xlow, xhigh, ylow, yhigh); univariate Q_i: each row is a vector of length 2 c(k1,k2) to define the sampling regions, i.e., low, middle, high).  Each subject should have n_i rows of the same values.
#' @param SampProb A matrix with the first dimension equal to sum(n_i).   Sampling probabilities from within each region (bivariate Q_i: each row is a vector of length 2 c(central region, outlying region); univariate Q_i: each row is a vector of length 3 with sampling probabilities for each region). Each subject should have n_i rows of the same values.
#' @param Weights Subject specific sampling weights.  A vector of length sum(n_i).  Not used unless using weighted Likelihood
#' @param ProfileCol the column number(s) for which we want fixed at the value of param.  Maimizing the log likelihood for all other parameters
#'                   while fixing these columns at the values of `params[ProfileCol]`
#' @param Keep.liC If TRUE outputs subject specific conditional log lileihoods to be used for the imputation procedure described in the AOAS paper keep z sum(n_i) by 2 design matric for random effects (intercept and slope)
#' @return The conditional log likelihood with a "gradient" attribute (if Keep.liC=FALSE) and subject specific contributions to the conditional likelihood if Keep.liC=TRUE).
#' @export
LogLikeCAndScore2 <- function(params, y, x, z, id, w.function, cutpoints, SampProb, Weights, ProfileCol=NA, Keep.liC=FALSE){
    npar   <- length(params)

    nbeta <- ncol(x)
    nVCsd <- ncol(z)
    nVCrho <- choose(nVCsd,2)
    nERRsd <- npar-nbeta-nVCsd-nVCrho

    beta.index   <- c(1:nbeta)
    vc.sd.index  <- nbeta + (c(1:nVCsd))
    vc.rho.index <- nbeta + nVCsd + (c(1:nVCrho))
    err.sd.index <- nbeta + nVCsd + nVCrho + c(1:nERRsd)

    beta     <- params[beta.index]
    sigma.vc <- exp(params[vc.sd.index])
    rho.vc   <- (exp(params[vc.rho.index])-1) / (exp(params[vc.rho.index])+1)
    sigma.e  <- exp(params[err.sd.index])

    out     <- LogLikeC2( y=y, x=x, z=z, w.function=w.function, id=id, beta=beta, sigma.vc=sigma.vc, rho.vc=rho.vc, sigma.e=sigma.e, cutpoints=cutpoints,
                          SampProb=SampProb, Weights=Weights, Keep.liC=Keep.liC)
    GRAD    <- LogLikeC.Score2(y=y, x=x, z=z, w.function=w.function, id=id, beta=beta, sigma.vc=sigma.vc, rho.vc=rho.vc, sigma.e=sigma.e, cutpoints=cutpoints,
                               SampProb=SampProb, Weights=Weights)
    ## Need to use the chain rule: note that params is on the unconstrained
    ## scale but GRAD was calculated on the constrained parameters
    GRAD[vc.sd.index]  <- GRAD[vc.sd.index]*exp(params[vc.sd.index])
    GRAD[vc.rho.index] <- GRAD[vc.rho.index]*2*exp(params[vc.rho.index])/((exp(params[vc.rho.index])+1)^2)
    GRAD[err.sd.index] <- GRAD[err.sd.index]*exp(params[err.sd.index])
    ## Force the gradient of the fixed parameter to be zero, so that it does not move
    if (!is.na(ProfileCol)) GRAD[ProfileCol] <- 0
    attr(out,"gradient") <- GRAD

    out
}


## If you do not want to use the ascertainment correction term in the conditional likelihood
## set all SampProb values equal to each other.  This would be the case if you were doing
## straightforward maximum likelihood (albeit computationally inefficient) or weighted likelihood.

#' Fitting function: ACML or WL for a linear mixed effects model (random intercept and slope)
#'
#' @param formula.fixed formula for the fixed effects (of the form y~x)
#' @param formula.random formula for the random effects (of the form ~z).  Right now this model only fits random intercept and slope models.
#' @param data data frame that should contain everything in formula.fixed, formula.random, id, and Weights.  It does not include: w.function, cutpoints, SampProb
#' @param id sum(n_i) vector of subject ids (a variable contained in data)
#' @param w.function sum(n_i) vector with possible values that include "mean" (mean of response series), "intercept" (intercept of the regression of Yi ~ zi where zi is the design matrix for the random effects (solve(t.zi* zi) * t.zi)`[1,]`), "intercept1"  (intercept of the regression of Yi ~ zi where zi is the design matrix for the random effects (solve(t.zi * zi) * t.zi)`[1,]`). "intercept2" (second intercept of the regression of the Yi ~
##zi where zi is the design matrix for the bivariate random effects (b10,b11,b20,b21) solve(t.zi * zi) * t.zi)`[3,]`), "slope" (slope of the regression of Yi ~ zi where zi is the design matrix for the random effects (solve(t.zi * zi) * t.zi)`[2,]`), "slope1" (slope of the regression of Yi ~ zi where zi is the design matrix for the random effects (solve(t.zi * zi) * t.zi)`[2,]`), "slope2" (second slope of the regression of the Yi ~ zi where zi is the design matrix for the bivariate random effects (b10,b11,b20,b21) solve(t.zi * zi) * t.zi)`[4,]`) "bivar" (intercept and slope of the regression of Yi ~ zi where zi is the design matrix for the random effects (solve(t.zi * zi) * t.zi)`[c(1,2),]`) "mvints" (first and second intercepts of the bivariate regression of the Yi ~ zi where zi is the design matrix for the bivariate random effects (b10,b11,b20,b21) solve(t.zi %*% zi) * t.zi)`[c(1,3),]`) "mvslps" (first and second slopes of the bivariate regression of the Yi ~ zi where zi is the design matrix for the bivariate random effects (b10,b11,b20,b21) solve(t.zi * zi) * t.zi)`[c(1,3),]`).  There should be one unique value per subject.  There should be one unique value per subject but n_i replicates of that value.  Note that w.function should NOT be in the dat dataframe.
#' @param InitVals starting values for c(beta, log(sigma0), log(sigma1), log((1+rho)/(1-rho)), log(sigmae))
#' @param cutpoints A matrix with the first dimension equal to sum(n_i).  These cutpoints define the sampling regions for individual subjects.  If using a low, medium, high, sampling scheme, this is a sum(n_i) by 2 matrix that must be a distinct object not contained in the dat dataframe.  Each row is a vector of length 2 c(k1,k2) to define the sampling regions, i.e., low, middle, high.  If using a square doughnut design this should be sum(n_i) by 4 matrix (var1lower, var1upper, var2lower, var2upper). Each subject should have n_i rows of the same values.
#' @param SampProb A matrix with the first dimension equal to sum(n_i).   Sampling probabilities from within each region. For low medium high sampling, each row is a vector of length 3 with sampling probabilities for each region. For bivariate stratum sampling each row is a vector of length 2 with sampling probabilities for the inner and outer strata. Each subject should have n_i rows of the same values.  Not in data.
#' @param Weights Subject specific sampling weights.  A vector of length sum(n_i).  This should be a variable in the data dataframe. It should only be used if doing IPWL.  Note if doing IPWL, only use robcov (robust variances) and not covar.  If not doing IPWL, this must be a vectors of 1s.
#' @param ProfileCol the column number(s) for which we want fixed at the value of param.  Maximizing the log likelihood for all other parameters
#'                   while fixing these columns at the values of `InitVals[ProfileCol]`
#' @return Ascertainment corrected Maximum likelihood: Ests, covar, LogL, code, robcov
#' @importFrom stats model.frame
#' @importFrom stats model.matrix
#' @importFrom stats model.response
#' @importFrom stats na.omit
#' @importFrom stats nlm
#' @importFrom stats pnorm
#' @export
acml.lmem2 <- function(formula.fixed,
                       formula.random,
                       data,
                       id,
                       w.function="mean",
                       InitVals,
                       cutpoints,
                       SampProb,
                       Weights,
                       ProfileCol=NA){              ## Columns to be held fixed while doing profile likelihood.  It is fixed at its initial value.

    if(is.null(formula.random)) {stop('Specify the random effects portion of the model.  It is currently NULL.')}
    if(!is.data.frame(data)) {
        data <- as.data.frame(data)
        warning('data converted to data.frame.')
    }
    terms = unique( c(all.vars(formula.fixed), all.vars(formula.random),
                      as.character(substitute(id)), as.character(substitute(Weights))) )
    data  = data[,terms]

    if(any(is.na(data))) data = na.omit(data)

    id0   =  as.character(substitute(id))
    id    = data$id = data[ , id0 ]

    fixed.f = model.frame(formula.fixed, data)
    fixed.t = attr(fixed.f, "terms")
    y      = model.response(fixed.f,'numeric')
    uy     = unique(y)
    x      = model.matrix(formula.fixed, fixed.f)

    rand.f = model.frame(formula.random, data)
    #z      = model.matrix(formula.random, rand.f) ####### changed Aug19, 2019
    z      = model.matrix(formula.random, fixed.f) ####### changed Aug19, 2019

    #if (is.na(SampProb[1])) SampProb = c(1,1,1)
    Weights0   =  as.character(substitute(Weights))
    Weights    = data$Weights = data[ , Weights0 ]

    acml.fit <- nlm(LogLikeCAndScore2, InitVals, y=y, x=x, z=z, id=id, w.function=w.function,
                    cutpoints=cutpoints, SampProb=SampProb, Weights=Weights, ProfileCol=ProfileCol,
                    stepmax=4, iterlim=250, check.analyticals = TRUE, print.level=2)

    ## Calculate the observed information and then invert to get the covariance matrix
    npar        <- length(acml.fit$estimate)
    Hessian.eps <- 1e-7
    eps.mtx     <- diag(rep(Hessian.eps, npar))
    grad.at.max <- acml.fit$gradient
    ObsInfo.tmp <- ObsInfo <- matrix(NA, npar, npar)

    ## Observed Information
    for (j in 1:npar){
        temp            <- LogLikeCAndScore2(acml.fit$estimate+eps.mtx[j,], y=y, x=x, z=z, id=id,
                                             w.function=w.function, cutpoints=cutpoints,
                                             SampProb=SampProb,Weights=Weights, ProfileCol=ProfileCol)
        ObsInfo.tmp[j,] <- (attr(temp,"gradient")-grad.at.max)/(Hessian.eps)
    }
    for (m in 1:npar){
        for (n in 1:npar){ ObsInfo[m,n] <-  (ObsInfo.tmp[m,n]+ObsInfo.tmp[n,m])/2}}

    ## Cheese part of the sandwich estimator
    nbeta <- ncol(x)
    nVCsd <- ncol(z)
    nVCrho <- choose(nVCsd,2)
    nERRsd <- npar-nbeta-nVCsd-nVCrho

    beta.index   <- c(1:nbeta)
    vc.sd.index  <- nbeta + (c(1:nVCsd))
    vc.rho.index <- nbeta + nVCsd + (c(1:nVCrho))
    err.sd.index <- nbeta + nVCsd + nVCrho + c(1:nERRsd)

    Cheese <- LogLikeC.Score2(y=y, x=x, z=z, w.function=w.function,
                              id=id, beta=acml.fit$estimate[beta.index],
                              sigma.vc=exp(acml.fit$estimate[vc.sd.index]),
                              rho.vc=   (exp(acml.fit$estimate[vc.rho.index])-1) / (exp(acml.fit$estimate[vc.rho.index])+1),
                              sigma.e=exp(acml.fit$estimate[err.sd.index]),
                              cutpoints=cutpoints,
                              SampProb=SampProb,
                              Weights=Weights,
                              CheeseCalc=TRUE)

    if (!is.na(ProfileCol)){
        acml.fit$estimate <- acml.fit$estimate[-ProfileCol]
        ObsInfo           <- ObsInfo[-ProfileCol, -ProfileCol]
        Cheese            <- Cheese[-ProfileCol, -ProfileCol]
    }


    out              <- NULL
    out$call         <- match.call()
    out$coefficients <- acml.fit$estimate
    out$covariance   <- solve(ObsInfo)
    out$robcov       <- solve(ObsInfo)%*%Cheese%*%solve(ObsInfo)
    out$logLik       <- -acml.fit$minimum
    attr(out,'args') <- list(formula.fixed=formula.fixed,
                             formula.random=formula.random,
                             id=id,
                             w.function=w.function,
                             cutpoints = cutpoints,
                             SampProb = SampProb,
                             Weights=Weights,
                             WeightsVar = Weights0,
                             ProfileCol=ProfileCol)
    if(kappa(out$covar) > 1e5) warning("Poorly Conditioned Model")
    out
}

#' Create a list of subject-specific data
#'
#' @param id sum(n_i) vector of subject ids
#' @param y response vector
#' @param x sum(n_i) by p design matrix for fixed effects
#' @param z sum(n_i) by 2 design matric for random effects (intercept and slope)
#' @param Weights Subject specific sampling weights.  A vector of length sum(n_i).  Not used unless using weighted Likelihood
#' @param SampProb A matrix with the first dimension equal to sum(n_i).   Sampling probabilities from within each region (bivariate Q_i: each row is a vector of length 2 c(central region, outlying region); univariate Q_i: each row is a vector of length 3 with sampling probabilities for each region). Each subject should have n_i rows of the same values.
#' @param w.function sum(n_i) vector with possible values that include "mean" "intercept" "slope" and "bivar."  There should be one unique value per subject
#' @param cutpoints A matrix with the first dimension equal to sum(n_i).  These cutpoints define the sampling regions (bivariate Q_i: each row is a vector of length 4 c(xlow, xhigh, ylow, yhigh); univariate Q_i: each row is a vector of length 2 c(k1,k2) to define the sampling regions, i.e., low, middle, high).  Each subject should have n_i rows of the same values.
#' @export
CreateSubjectData <- function(id,y,x,z,Weights,SampProb,cutpoints,w.function){
    id.tmp        <- split(id,id)
    y.tmp         <- split(y,id)
    x.tmp         <- split(x,id)
    z.tmp         <- split(z,id)
    Weights.tmp <- split(Weights,id)
    SampProb.tmp  <- split(SampProb,id)
    cutpoints.tmp  <- split(cutpoints,id)
    w.function.tmp  <- split(w.function,id)

    subjectData <- vector('list', length=length(unique(id)))
    subjectData <- list()
    uid <- as.character(unique(id))
    for(j in seq(along=uid)){
        i <- uid[j]
        subjectData[[j]] <- list(idi          = as.character(unique(id.tmp[[i]])),
                                 xi           = matrix(x.tmp[[i]], ncol=ncol(x)),
                                 zi           = matrix(z.tmp[[i]], ncol=ncol(z)),
                                 yi           = y.tmp[[i]],
                                 Weights.i    = unique(Weights.tmp[[i]]),
                                 SampProb.i   = matrix(SampProb.tmp[[i]], ncol=ncol(SampProb))[1,],
                                 w.function.i = as.character(unique(w.function.tmp[[i]])),
                                 cutpoints.i  = matrix(cutpoints.tmp[[i]], ncol=ncol(cutpoints))[1,])
    }
    names(subjectData) <- uid
    subjectData
}
