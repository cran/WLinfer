#' Estimating parameters in weighed Lindley distribution
#'
#' Parameter estimate functions: \code{MME_WL}, \code{MMEm_WL}, \code{MLE_WL} and \code{MLEc_WL} can be
#' used to estimate parameters in weighed Lindley distribution. \code{MME_WL}, \code{MMEm_WL}
#' and \code{MLEc_WL} have closed form values for parameters, lambda and phi. On the other hands,
#' parameter lambda in \code{MLE_WL} is
#' based on numerical solution. It use the function \code{nleqslv} for solving one variable equation.
#'
#' @details
#' These functions implement formulas given in Hyoung-Moon Kim. et al. (2020).
#'
#' @param x a numeric vector.
#'
#' @return A numeric vector of lambda and phi estimated by each method.
#'
#' @examples
#' data <- fail_fiber
#' mme <- MME_WL(data)
#' modified_mme <- MMEm_WL(data)
#' mle <- MLE_WL(data, mme[2])
#' mlec <- MLEc_WL(data)
#'
#' rbind(mme, modified_mme, mle, mlec)
#'
#' @references Hyoung-Moon Kim. and Yu-Hyeong Jang. (2020). New Closed-Form Estimators for
#' Weighted Lindley Distribution. \emph{ }, submitted.
#'
#'
#' @importFrom LindleyR dwlindley pwlindley qwlindley
#' @importFrom pracma integral
#' @importFrom boot boot
#' @importFrom nleqslv nleqslv
#' @importFrom stats ks.test integrate qnorm quantile sd var ppoints pchisq nlminb
#' @importFrom grDevices dev.interactive devAskNewPage
#' @importFrom graphics boxplot hist legend lines par plot text abline contour points
#' @importFrom goftest ad.test cvm.test
#' @importFrom cubature hcubature
#' @importFrom bbmle mle2
#'
#' @export
MME_WL <- function(x){ #x is the variable with vector form
  x_bar <- mean(x)
  n <- length(x)
  S <- sd(x)*sqrt((n-1)/n)
  g <- S^4-x_bar*(x_bar^3+2*x_bar^2+x_bar-4*S^2)

  phihat <- (-g + sqrt( g^2+(16*S^2)*(S^2+(x_bar+1)^2)*x_bar^3   ))/(2*(S^2)*( S^2 + (x_bar+1)^2 ))
  lambdahat <- (-phihat*(x_bar-1)+sqrt( (phihat*(x_bar-1))^2 + 4*x_bar*phihat*(phihat+1) ) )/(2*x_bar)
  est <- matrix(c(lambdahat, phihat), nrow=1)
  colnames(est) <- c("lambda","phi")
  rownames(est) <- "MME"
  return(est)
}

#' @rdname MME_WL
#' @export
MMEm_WL <- function(x){ #x is the variable with vector form
  x_bar = mean(x)
  x_starbar = mean(1/(1+x))
  lambda_mm = x_starbar*(1-x_starbar)/(x_bar*x_starbar-(1-x_starbar))
  phi_mm = (1-x_starbar)^2/(x_bar*x_starbar-(1-x_starbar))

  est <- matrix(c(lambda_mm,phi_mm), nrow=1)
  colnames(est) <- c("lambda","phi")
  rownames(est) <- "modified MME"
  return(est)
}

#' @rdname MME_WL
#' @export
MLEc_WL <- function(x){ #x is the variable with vector form
  n = length(x)
  x_bar = mean(x)
  a1 = sum(x*log(x))/sum(log(x))
  a0 = ( n + sum(x*log(x)/(1+x)) )/sum(log(x))
  d <- a0 + ((1-a0)*(1+a1)-1)/(x_bar-a1)
  lambda = (  d + sqrt( d^2
                        +4*(a1+1)*a0*(a0-1)/(x_bar-a1) )  ) / (2*(a1+1))
  phi = a1*lambda - a0

  est <- matrix(c(lambda,phi), nrow=1)
  colnames(est) <- c("lambda","phi")
  rownames(est) <- "Closed form MLE"
  return(est)
}

#' MLE in weighed Lindley distribution
#'
#' \code{MLE_WL} returns maximum likelihood estimates of parameters.
#' First, the estimate of phi is obtained from one dimensional non linear equation.
#' Then, by plugging in the estimated phi, the estimate of lambda is easily obtained.
#'
#' @details
#' These functions implement formulas given in Hyoung-Moon Kim. et al. (2020).
#'
#' @param y a numeric vector of observations.
#' @param init a initial value for estimate phi.
#'
#' @return A numeric vector of MLE for lambda and phi.
#'
#' @examples
#' data <- fail_fiber
#' mme <- MME_WL(data)
#' modified_mme <- MMEm_WL(data)
#' mle <- MLE_WL(data, mme[2])
#' mlec <- MLEc_WL(data)
#'
#' rbind(mme, modified_mme, mle, mlec)
#'
#' @references Hyoung-Moon Kim. and Yu-Hyeong Jang. (2020). New Closed-Form Estimators for
#' Weighted Lindley Distribution. \emph{ }, submitted.
#'
#' @export
MLE_WL <- function(y, init){
  ml_phi = nleqslv(init,function(y,ml_phi){
    -(log((-ml_phi*(mean(y)-1)+sqrt( (ml_phi*(mean(y)-1))^2 + 4*ml_phi*(ml_phi+1)*mean(y)  ))/(2*mean(y)))-
        ((-ml_phi*(mean(y)-1)+sqrt( (ml_phi*(mean(y)-1))^2 + 4*ml_phi*(ml_phi+1)*mean(y)  ))/(2*mean(y))+ml_phi)^(-1)+
        mean(log(y))-digamma(ml_phi))
  },y=y,control=list(btol=.01))$x
  ml_lam = (-ml_phi*(mean(y)-1)+sqrt( (ml_phi*(mean(y)-1))^2 + 4*ml_phi*(ml_phi+1)*mean(y)  ))/(2*mean(y))

  est <- matrix(c(ml_lam,ml_phi), nrow=1)
  colnames(est) <- c("lambda","phi")
  rownames(est) <- "MLE"
  return(est)
}

#' Cox & Snell bias correction methods for estimators
#'
#' \code{CoxSnell_bias} and \code{CoxSnell_bias_log} provide a vector of MLE bias, and it can be used for
#' bias correction of both MLE and MLEc. The estimation method was suggested by by Cox and Snell(1968).
#'
#' @details
#' \code{CoxSnell_bias} provides the bias of original lambda and phi.
#' However, \code{CoxSnell_bias_log} provides the bias of log lambda and log phi.
#' In some cases, estimators are smaller than bias, and it means that the bias corrected estimators
#' are out of parameter space. To solve this problem, \code{CoxSnell_bias_log} is useful.
#' Correction formula is based on Fisher information and some cumulants,
#' which are given in Kim and Jang (2020).
#'
#' @param n a numeric value of data length.
#' @param lambda a numeric value of estimated lambda.
#' @param phi a numeric value of estimated phi.
#'
#' @return A numeric vector of Cox & Snell biases of lambda and phi.
#'
#' @section Background:
#' These functions implement formulas given in Hyoung-Moon Kim. et al. (2020).
#'
#' @references Hyoung-Moon Kim. and Yu-Hyeong Jang. (2020). New Closed-Form Estimators for
#' Weighted Lindley Distribution. \emph{ }, submitted.
#'
#' @examples
#' data <- fail_fiber
#' mlec <- MLEc_WL(data)
#' n <- length(data)
#' CoxSnell_bias(n, mlec[1], mlec[2])
#' CoxSnell_bias_log(n, mlec[1], mlec[2])
#'
#' @export
CoxSnell_bias <- function(n,lambda,phi){
  info_K <- matrix(0, ncol=2, nrow=2)
  info_K[1,1] <- (phi+1)/(lambda^2) - 1/(lambda+phi)^2
  info_K[1,2] = info_K[2,1] <- -1/lambda - 1/(lambda+phi)^2
  info_K[2,2] <- -1/(lambda+phi)^2 + trigamma(phi)
  K <- n*info_K

  cox_A1 = cox_A2 = matrix(0,nrow=2,ncol=2)
  cox_A1[1,1] = 1*(phi+1)/lambda^3-1/(lambda+phi)^3
  cox_A1[1,2] = cox_A1[2,1] = -1/(2*lambda^2)-1/(lambda+phi)^3
  cox_A1[2,2] = -1/(lambda+phi)^3

  cox_A2[1,1] = -1/(2*lambda^2)-1/(lambda+phi)^3
  cox_A2[1,2] = cox_A2[2,1] = -1/(lambda+phi)^3
  cox_A2[2,2] = -1/(lambda+phi)^3-1*psigamma(phi,2)/2
  A <- n*cbind(cox_A1,cox_A2)

  K_inv = solve(K)
  bias <- K_inv%*%A%*%as.vector(K_inv)
  rownames(bias) <- c("lambda","phi")
  colnames(bias) <- "bias"
  return(t(bias))
}

#' @rdname CoxSnell_bias
#' @export
CoxSnell_bias_log <- function(n,lambda, phi){  #return estimated bias of log(lambdahat) and log(phihat)

  k11 = (-(phi+1)/lambda^2 + 1/(lambda+phi)^2)*n
  k22 = ( 1/(lambda+phi)^2 - trigamma(phi))*n
  k12 = k21 = (1/lambda + 1/(lambda+phi)^2)*n

  K2 = matrix(0,ncol=2,nrow=2)
  K2[1,1] = k11*lambda^2
  K2[2,2] = k22*phi^2
  K2[1,2] = K2[2,1] = k21*lambda*phi
  CK_inv = solve(-K2)


  k111 = (2*(phi+1)/lambda^3 - 2/(lambda+phi)^3)*n
  k112 = (-1/lambda^2-2/(lambda+phi)^3)*n
  k122 = (-2/(lambda+phi)^3)*n
  k222 = (-2/(lambda+phi)^3 - psigamma(phi,2) )*n

  cox_A1 = cox_A2 = matrix(0,nrow=2,ncol=2)

  cox_A1[1,1] = lambda^3*(k111)+ lambda^2*k11
  cox_A1[1,2] = cox_A1[2,1] = lambda*phi*(lambda*k112+k12 )
  cox_A2[1,1] = lambda*phi*(lambda*k112-k12 )
  cox_A2[1,2] = cox_A2[2,1] = lambda*phi*(phi*k122+k12 )
  cox_A1[2,2] = lambda*phi*(phi*k122-k12)
  cox_A2[2,2] = phi^3*(k222)+ phi^2*k22

  CA = cbind(cox_A1,cox_A2)/2

  bias <- CK_inv %*%CA %*% matrix(as.vector(CK_inv),ncol=1)
  rownames(bias) <- c("lambda","phi")
  colnames(bias) <- "bias"
  return(t(bias))
}

#' Firth's method for bias correction of estimators
#'
#' Firth's method is bias correction option included in \code{WL}. \code{Firth_WL} and \code{Firth_WL_log} provides
#' the bias corrected lambda and phi based on Firth's method.
#'
#' @details
#' \code{Firth_WL} and \code{Firth_WL_log} returns a vector of estimates of parameters corrected by Firth's method
#' which uses the modified likelihood equations. In case of weighted Lindley distribution,
#' two non-linear equations should be solved since the solutions are not in closed form.
#' To this end, R package \code{nleqslv} is used.
#' To avoid poor local maxima, other estimators like MMEm or MLEc are recommended
#' to be used as initial values.
#'
#' @param y a numeric vector.
#' @param init a vector of initial values for iterative algorithm designed
#' to solve the modified likelihood equations.
#'
#' @return A vector of corrected estimators lambda and phi.
#'
#' @section Background:
#' Non-linear equations to be solved are derived in Kim and Jang (2020).
#'
#' @references Hyoung-Moon Kim. and Yu-Hyeong Jang. (2020). New Closed-Form Estimators for
#' Weighted Lindley Distribution. \emph{ }, submitted.
#'
#' @examples
#' data <- fail_fiber
#' Firth_WL(data,MMEm_WL(data))
#' Firth_WL_log(data,MMEm_WL(data))
#'
#' @export
Firth_WL <- function(y,init){

  temp_firth <- function(pars,y){
    n <- length(y)
    a1 <- c(n*(pars[2]+1)/pars[1]-n/(pars[1]+pars[2])-sum(y),
          n*log(pars[1])-n/(pars[1]+pars[2])-n*digamma(pars[2])+sum(log(y)) )

    cox_A1 = cox_A2 = matrix(0,nrow=2,ncol=2)
    cox_A1[1,1] = 1*(pars[2]+1)/pars[1]^3-1/(pars[1]+pars[2])^3
    cox_A1[1,2] = cox_A1[2,1] = -1/(2*pars[1]^2)-1/(pars[1]+pars[2])^3
    cox_A1[2,2] = -1/(pars[1]+pars[2])^3


    cox_A2[1,1] =  -1/(2*pars[1]^2)-1/(pars[1]+pars[2])^3
    cox_A2[1,2] = cox_A2[2,1] = -1/(pars[1]+pars[2])^3
    cox_A2[2,2] = -1/(pars[1]+pars[2])^3-1*psigamma(pars[2],2)/2
    CA <- n * cbind(cox_A1,cox_A2)

    info_K = matrix(0,ncol=2,nrow = 2)
    info_K[1,1] = 1*(pars[2]+1)/(pars[1]^2) - 1/(pars[1]+pars[2])^2
    info_K[1,2] = info_K[2,1] = -1/pars[1]-1/(pars[1]+pars[2])^2
    info_K[2,2] = -1/(pars[1]+pars[2])^2 +1*trigamma(pars[2])
    K <- n*info_K

    a2= c(CA%*%matrix(as.vector(solve(K)),ncol=1))
    return(a1-a2)
  }

  result = nleqslv(init,temp_firth,y=y,control=list(btol=.01))

  est <- matrix(result$x, nrow=1)
  colnames(est) <- c("lambda","phi")
  rownames(est) <- "MLE(Firth)"
  return(est)
}

#' @rdname Firth_WL
#' @export
Firth_WL_log = function(y,init){

  C_K2 = function(y,est){
    l = est[1]
    p = est[2]

    k11 = (-(p+1)/l^2 + 1/(l+p)^2)*length(y)
    k22 = ( 1/(l+p)^2 - trigamma(p))*length(y)
    k12 = k21 = (1/l + 1/(l+p)^2)*length(y)

    result = matrix(0,ncol=2,nrow=2)
    result[1,1] = k11*l^2
    result[2,2] = k22*p^2
    result[1,2] = result[2,1] = k21*l*p
    result = -result
    result
  }

  C_A2 = function(y,est){
    l = est[1]
    p = est[2]

    k11 = (-(p+1)/l^2 + 1/(l+p)^2)*length(y)
    k22 = ( 1/(l+p)^2 - trigamma(p))*length(y)
    k12 =  k21 = (1/l + 1/(l+p)^2)*length(y)

    k111 = (2*(p+1)/l^3 - 2/(l+p)^3)*length(y)
    k112 = (-1/l^2-2/(l+p)^3)*length(y)
    k122 = (-2/(l+p)^3)*length(y)
    k222 = (-2/(l+p)^3 - psigamma(p,2) )*length(y)


    cox_A1 = cox_A2 = matrix(0,nrow=2,ncol=2)

    cox_A1[1,1] = l^3*(k111)+ l^2*k11
    cox_A1[1,2] = cox_A1[2,1] = l*p*(l*k112+k12 )
    cox_A2[1,1] = l*p*(l*k112-k12 )
    cox_A2[1,2] = cox_A2[2,1] = l*p*(p*k122+k12 )
    cox_A1[2,2] = l*p*(p*k122-k12)
    cox_A2[2,2] = p^3*(k222)+ p^2*k22

    result = cbind(cox_A1,cox_A2)/2
    return(result)
  }

  temp_firth = function(pars,y){
    a1= exp(pars)*c(length(y)*(exp(pars)[2]+1)/exp(pars)[1]-length(y)/(exp(pars)[1]+exp(pars)[2])-sum(y),
                    length(y)*log(exp(pars)[1])-length(y)/(exp(pars)[1]+exp(pars)[2])-length(y)*digamma(exp(pars)[2])+sum(log(y)) )
    a2= c(C_A2(y,exp(pars))%*%matrix(as.vector(solve(C_K2(y,exp(pars)))),ncol=1))
    return(sum((a1-a2)^2))
  }

  result = nlminb(init,temp_firth,y=y)$par
  return(exp(result))
}


#' Asymptotic covariance matrix of estimators
#'
#' These four functions: \code{MME_var},
#' \code{MMEm_var}, \code{MLE_var} and \code{MLEc_var} provide asymptotic covariance matrixs for MME,
#' modified MME, MLE, MLEc, repectively. All of these can be calculated to a closed form value.
#'
#' @param l a numeric value.
#' @param p a numeric value.
#' @param n a numeric value.
#'
#' @details
#' These functions implement formulas given in Hyoung-Moon Kim. et al. (2020).
#'
#' @return A matrix of asymptotic covariance of lambda and phi.
#'
#'
#' @references Hyoung-Moon Kim. and Yu-Hyeong Jang. (2020). New Closed-Form Estimators for
#' Weighted Lindley Distribution. \emph{ }, submitted.
#'
#' @examples
#' data <- fail_fiber
#' n <- length(data)
#' mme <- MME_WL(data)
#' modified_mme <- MMEm_WL(data)
#' mle <- MLE_WL(data, mme[2])
#' mlec <- MLEc_WL(data)
#' MME_var(mme[1],mme[2],n)
#' MMEm_var(modified_mme[1],modified_mme[2],n)
#' MLE_var(mle[1],mle[2],n)
#' MLEc_var(mlec[1],mlec[2],n)
#'
#' @export
##### Asymptotic variance #####
MME_var <- function(l,p,n){
  H <- matrix(0,ncol=2,nrow=2)
  H[1,1] <- (l^2-(p+1)*(l+p)^2)/(l*(l+p))^2
  H[1,2] <- (l+(l+p)^2)/(l*(l+p)^2)
  H[2,1] <- -(2*p*(p+1)*(3*l+2*p+(l+p)^2))/(l^3*(l+p)^2)
  H[2,2] <- ((2*p+3)*(l+p)^2-2*l*(l-1))/(l*(l+p))^2
  W <- solve(H)

  sig.m <- matrix(0,ncol=2,nrow=2)
  sig.m[1,1] <- l^2/(2*p)*((l+p)^2-l^2/(p+1))
  sig.m[1,2] = sig.m[2,1] <- l*((l+p)^2+3*l+2*p)
  sig.m[2,2] <- 2*p^3+(4*l+9)*p^2+2*(l^2+7*l+5)*p+3*l*(l+4)
  sig.m <- sig.m*(2*p*(p+1))/l^4/(l+p)^2

  result <- W%*%sig.m%*%t(W)/n
  rownames(result)  = colnames(result) <- c("lambda","phi")
  return(result)
}

#' @rdname MME_var
#' @export
MMEm_var <- function(l,p,n){
  sig.mm <- matrix(0,ncol=2,nrow=2)
  sig.mm[1,1] <- (1+p)*(l+p)^2/l^2-1
  sig.mm[1,2] = sig.mm[2,1] <- -p
  sig.mm[2,2] <- l^(p+1)*(l+p)*hcubature(function(x) x^(p-1)*exp(-l*x)/(1+x),0,Inf)$integral/gamma(p) -l^2

  sig.mm <- sig.mm/(l+p)^2
  W <- cbind(-c(l^2,(l+p)^2+l)*(l+p)/p,-c(l^2,(l+p)^2+2*l+p)*(l+p)/l)
  result <- crossprod(W, sig.mm)%*%W /n #result <- t(W)%*%sig.mm%*%W /n

  rownames(result) = colnames(result) <- c("lambda","phi")
  return(result)
}

#' @rdname MME_var
#' @export
MLE_var <- function(l,p,n){
  cox_K <- matrix(0,ncol=2,nrow = 2)
  cox_K[1,1] <- (p+1)/(l^2) - 1/(l+p)^2
  cox_K[1,2] <- cox_K[2,1] <- -1/l-1/(l+p)^2
  cox_K[2,2] <- -1/(l+p)^2 + trigamma(p)

  result <- solve(cox_K)/n
  rownames(result) = colnames(result) <- c("lambda","phi")
  return(result)
}

#' @rdname MME_var
#' @export
MLEc_var <- function(l,p,n){
  digam_p <- digamma(p)
  digam_p1 <- digamma(p+1)
  digam_p2 <- digamma(p+2)
  logl <- log(l)

  EX = p*(l+p+1)/(l*(l+p))
  EY = digam_p-logl + 1/(l+p)
  EZ = p/l/(l+p)*( 1+(l+p+1)*(digam_p1-logl) )
  EW = p/(l+p)*(digam_p1-logl)

  EXX = p*(p+1)*(l+p+2)/(l^2*(l+p))
  EYY = (digam_p-logl)^2+2*(digam_p-logl)/(l+p)+trigamma(p)
  EZZ = p*(p+1)*(l+p+2)/(l^2*(l+p))*( (digam_p2-logl)^2+2*(digam_p2-logl)/(l+p+2)+trigamma(p+2) )
  EWW = ( l^(p+1)/(gamma(p)*(l+p)) )*integrate(function(x) x^{(p+2)-1}*(log(x)^2)*exp(-l*x)/(1+x), lower=0, upper=Inf)$value

  #EXY = EZ
  EXZ = p*(p+1)/l^2/(l+p)*( (l+p+2)*(digam_p2-logl)+1  )

  EXW = p*(p+1)/l/(l+p)*(digam_p2-logl)

  EYZ = p*(l+p+1)/l/(l+p)*( (digam_p1-logl)^2+2*(digam_p1-logl)/(l+p+1)+trigamma(p+1) )

  EYW = p/(l+p)*( (digam_p1-logl)^2+trigamma(p+1) )

  EZW = p*(p+1)/l/(l+p)*( (digam_p2-logl)^2 +trigamma(p+2))

  Cov = matrix(0,ncol=4,nrow=4)
  Cov[1,1] <- EXX - EX^2
  Cov[1,2] = Cov[2,1] <- EZ - EX*EY
  Cov[1,3] = Cov[3,1] <- EXZ-EX*EZ
  Cov[1,4] = Cov[4,1] <- EXW-EX*EW
  Cov[2,2] = Cov[2,2] <- EYY-EY*EY
  Cov[2,3] = Cov[3,2] <- EYZ-EY*EZ
  Cov[2,4] = Cov[4,2] <- EYW-EY*EW
  Cov[3,3] = Cov[3,3] <- EZZ-EZ*EZ
  Cov[3,4] = Cov[4,3] <- EZW-EZ*EW
  Cov[4,4] = Cov[4,4] <- EWW-EW*EW

  a <- EZ/EY; b <- (EW+1)/EY
  ay <- -EZ/(EY^2) ; az=1/EY ;
  by <- -(EW+1)/(EY^2) ; bw=1/EY

  g_one =  (b + ((1-b)*(1+a)-1)/(EX-a))/(2*(a+1))+
    sqrt( (b+((1-b)*(1+a)-1)/(EX-a))^2 +
            4*(a+1)*b*(b-1)/(EX-a) )/(2*(a+1))
  g_two = (g_one*EZ-EW-1)/EY


  B=b+((1-b)*(a+1)-1)/(EX-a)
  C=b*(b-1)*(a+1)/(EX-a)
  U = B+sqrt(B^2+4*C)

  Bx = -((1-b)*(a+1)-1)/(EX-a)^2 ;Bx
  By =  by+((ay-by-ay*b-a*by)*(EX-a) + (a-a*b-b)*ay )/(EX-a)^2 ; By
  Bz = (az*(1-b)*(EX-a)+(a-a*b-b)*az)/(EX-a)^2;Bz
  Bw = bw-bw*(a+1)/(EX-a) ; Bw
  Cx = -b*(b-1)*(a+1)/(EX-a)^2
  Cy = ((by*(2*b-1)*(a+1)+b*(b-1)*ay)*(EX-a) + b*(b-1)*(a+1)*ay)/(EX-a)^2
  Cz = b*(b-1)*az*(EX+1)/(EX-a)^2
  Cw = (a+1)*bw*(2*b-1)/(EX-a)

  U_x = Bx + 0.5*((B^2+4*C)^(-0.5))*(2*B*Bx+4*Cx)
  U_y = By + 0.5*((B^2+4*C)^(-0.5))*(2*B*By+4*Cy)
  U_z = Bz + 0.5*((B^2+4*C)^(-0.5))*(2*B*Bz+4*Cz)
  U_w = Bw + 0.5*((B^2+4*C)^(-0.5))*(2*B*Bw+4*Cw)

  AA=matrix(0,nrow=2,ncol=4)
  AA[1,1] = U_x/2/(a+1)
  AA[1,2] = (U_y*(a+1)-U*ay)/(2*(a+1)^2)
  AA[1,3] = (U_z*(a+1)-U*az)/(2*(a+1)^2)
  AA[1,4] = U_w/2/(a+1)
  AA[2,1] = a*AA[1,1]
  AA[2,2] = ay*g_one+a*AA[1,2]-by
  AA[2,3] = az*g_one + a*AA[1,3]
  AA[2,4] = a*AA[1,4]-bw

  result <- AA%*%Cov%*%t(AA)/n
  rownames(result) = colnames(result) <- c("lambda","phi")
  return(result)
}

#' Estimator test based on Wilks' theorem
#'
#' This is a test based on Wilks' theorem,
#' to determine which parameter space the estimated parameter is included in.
#'
#'
#' @param x a numeric vector or data frame.
#' @param estimator a numeric vector with estimated lambda and phi.
#' @param side a character string which selects the direction of wilks' theorem test
#' ("two", "less" or "greater").
#'
#' @details By Wilks' theorem, we can test the k-dimensional parameter with chi-square distribution.
#' The Wilks' theorem test can be performed by assuming the parameter space of null hypothesis
#' and setting the part not included in it as the parameter space of the alternative hypothesis.
#'
#' @return \code{wilks.test} returns a list with these components:
#' \item{side}{a character string of one of "two", "less" or "greater".}
#' \item{stat}{a numeric value the statistics from Wilks' theorem.}
#' \item{pvalue}{a numeric value the p-value of the statistics.}
#'
#' @examples
#' data <- fail_fiber
#' mme <- MME_WL(data)
#' wilks.test(data, mme, side="two")
#' wilks.test(data, mme, side="less")
#' wilks.test(data, mme, side="greater")
#'
#' @export
wilks.test <- function(x, estimator, side="two"){
  if(side!="less" & side!="greater"){
    if(side!="two"){message("Warning: wrong side input. It runs by default.")}
    side <- "two"
  }
  p <- length(estimator)
  #suppressWarnings()


  if(side=="greater"){
    est_ml <- MLE_WL(x, MME_WL(x)[2])
    LL <- function(lambda,phi){ -sum( log(dwlindley(x,lambda,phi)) ) }
    est_test <- mle2(LL, method="L-BFGS-B", upper=estimator, lower=c(1e-10, 1e-10),
                     start=list(lambda=estimator[1], phi=estimator[2]) )@coef
    nu <- sum( log( dwlindley(x,est_ml[1],est_ml[2]) ) )
    de <- sum( log( dwlindley(x,est_test[1],est_test[2]) ) )
    stat <- 2*(nu - de)
  }else if(side=="less"){
    est_ml <- MLE_WL(x, MME_WL(x)[2])
    LL <- function(lambda,phi){ -sum( log(dwlindley(x,lambda,phi)) ) }
    est_test <- mle2(LL, method="L-BFGS-B", lower=estimator,
                     start=list(lambda=estimator[1], phi=estimator[2]) )@coef
    nu <- sum( log( dwlindley(x,est_ml[1],est_ml[2]) ) )
    de <- sum( log( dwlindley(x,est_test[1],est_test[2]) ) )
    stat <- 2*(nu - de)
  }else{
    est_test <- MLE_WL(x, MME_WL(x)[2])
    nu <- sum( log( dwlindley(x,est_test[1],est_test[2]) ) )
    de <- sum( log( dwlindley(x,estimator[1],estimator[2]) ) )
    stat <- 2*(nu - de)
  }

  pvalue <- pchisq(stat,p,lower.tail=F)

  result <- list("side"=side, "stat"=stat, "pvalue"=pvalue)
  return(result)
}

#' Statistical inference in weighted Lindley distribution
#'
#' Weighted Lindley distribution is suggested by Ghitany et al. (2011). \code{WL} provides four types of
#' estimator, which are MME, modified MME, MLE, and MLEc from weighted Lindley distribution.
#' And there are four sub-options, which are bias-correction, goodness of fit test, confidence interval,
#' and Wilks' theorem test.
#'
#' @param x a numeric vector or data frame.
#' @param est_method a character string which selects the estimation method
#' ("MME", "MMEm", "MLE", "MLEc"), default is "MLEc".
#' @param bias_cor an optional character character string which selects the bias correction method
#' ("coxsnell", "boots" or "firth").
#' @param dist_test a character string or character vector which choose the test of goodness of fit
#' ("all","ks","ad","cvm").
#' @param gof_alpha a numeric value between 0 and 1 for controlling the significance level
#' of goodness of fit test; default value is 0.05.
#' @param ks_side a character string which selects the alternative hypothesis
#' ("two", "less" or "greater") for Kolmogorov-Smirnov Test, default is "two".
#' @param CI_method a character string which selects the method for calculating
#' confidence intervals ("asymp" or "boots"), default is "asymp".
#' Since the "asymp" option is not available with bias correction, only the "boots" is available
#' with bias correction.
#' @param CI_scale a character string which selects the scale of confidence intervals ("exp" or "normal")
#' @param CI_side a character string which selects the direction of confidence intervals
#' ("two", "less" or "greater").
#' @param CI_alpha a numeric value between 0 and 1 for controlling the significance level
#' of confidence intervals; default value is 0.05.
#' @param boots_iter a numeric value for iteration number of bootstrap method.
#' @param wilks_test logical. If \code{TRUE}, wilks' theorem test is performed.
#' @param wilks_alpha a numeric value between 0 and 1 for controlling the significance level
#' of wilks' theorem test; default value is 0.05.
#' @param wilks_side a character string which selects the direction of wilks' theorem test
#' ("two", "less" or "greater").
#'
#' @details
#' First, the user can determine the type of estimator from MME, modified MME, MLE, and MLEc.
#' The closed form formulas for MME, modified MME, and MLEc are given in Hyoung-Moon Kim. et al. (2020).
#' And MLE is obtained numerically. Additionally MLE and MLEc have bias correction options.
#' MLE has Cox&Snell method and Firth's method, however MLEc has Cox&Snell method and bootstrap method.\cr
#'
#' Second, it provides a goodness of fit test. There are three kinds of tests, Kolmogorov-Smirnov test,
#' Anderson Darling test, and Cramer-von Mises test. They provide statistics and also p-values.
#' If the input value \code{gof_alpha} is selected, it determines whether or not to reject the null hypothesis.
#' \cr
#'
#' Third, it provides information on the confidence interval. There are two kinds of confidence intervals,
#' one is based on bootstrap method, and the other is asymptotic variance based method.
#' Asymptotic variance based method is only available without bias correction, however bootstrap method is
#' always available. Sometimes the confidence interval is outside the parameter space. If it occers,
#' confidence interval will be calculated with log scale and show the exponential confidence interval of log
#' scaled estimators. This option can also be used separately with selecting "exp" in \code{CI_scale}.
#'
#' Lastly, through \code{wilks.test}, \code{WL} test the parameter space of estimators lambda and phi.
#' There is an option for Wilks' theorem test, and that option provides options for the side of
#' Wilks' theorem test.
#'
#'
#'
#' @return \code{WL} returns a list with these components:
#' \item{data}{a numeric vector the input values.}
#' \item{dataname}{a character string the name of input values.}
#' \item{stat_summary}{a numeric vector with min, 1st quantile, median, 3rd quantile, and max.}
#' \item{mean}{a numeric value mean of input values.}
#' \item{var}{a numeric value variance of input values.}
#' \item{est}{a numeric vector with estimated lambda and phi.}
#' \item{lambda_var}{a numeric value variance of estimated lambda.}
#' \item{phi_var}{a numeric value variance of estimationed lambda.}
#' \item{bias_cor}{a character string from bias correction method ("coxsnell","firth" or "boots").}
#' \item{est_method}{a character string from estimation method ("MME", "MMEm", "MLE" or "MLEc").}
#' \item{boots_iter}{a numeric value of bootstrap iteration.}
#' \item{test_list}{a list with results of goodness of fit test.}
#' \item{CI_list}{a list with confidence interval related outputs.}
#' \item{wilks_list}{a list with results of wilks' test.}
#'
#' @references Ghitany, M., Alqallaf, F., Al-Mutairi, D., Husain, H. (2011). A two-parameter weighted
#' Lindley distribution and its applications to survival data. \emph{Mathematics and Computers in
#' Simulation} 81: 1190-1201.
#' @references Hyoung-Moon Kim. and Yu-Hyeong Jang. (2020). New Closed-Form Estimators for
#' Weighted Lindley Distribution. \emph{ }, submitted.
#' @references Wang, M., Wang, W. (2017). Bias-Corrected maximum likelihood estimation of the parameters
#' of the weighted Lindley distribution. \emph{Communications in Statistics
#' Simulation and Computation} 46: 530-545.
#'
#' @examples
#' example <- lifetime_alum
#' result <- WL(example)
#' print(result)
#'
#' @export
WL <- function(x, est_method="MLEc", bias_cor="None",
               dist_test="ks", gof_alpha=0.05, ks_side="two",
               CI_method="asymp", CI_scale="normal", CI_side="two", CI_alpha=0.05,
               boots_iter=10^3, wilks_test=TRUE, wilks_alpha=0.05, wilks_side="two"){

  # Organize input values
  data <- x
  dataname <- deparse(substitute(x))
  stat_summary <- matrix(quantile(x), nrow=1)
  meanx <- mean(x); varx <- var(x)
  n <- length(x)
  remove(x)

  # CI related
  if(CI_method!="asymp" & CI_method!="boots"){
    message("The CI_method has wrong value. It runs by 'asymp' option.")
    CI_method <- "asymp"
  }
  if(CI_scale!="normal" & CI_scale!="exp"){
    message("The scale of CI is wrong. It runs by 'normal' option.")
    CI_scale <- "normal"
  }
  if(CI_alpha>1 | CI_alpha<0){
    message("The significance level of confidence intervals must be a number between 0 and 1.
            It runs by alpha=0.05.")
    CI_alpha <- 0.05
  }
  CI_per <- (1-CI_alpha)*100
  if(CI_side=="less"){ci <- c(0,1-CI_alpha)} else if(CI_side=="greater"){ci <- c(CI_alpha,1)
  }else{
    if(CI_side != "two"){ message("The CI_side has wrong value. It runs by default."); CI_side <- "two"}
    ci <- c(CI_alpha/2,1-CI_alpha/2)
  }

  # Distribution test related
  if(sum(dist_test %in% "all")==1){
    dist_vec <- c(1,1,1)
  }else{
    ks <- sum(dist_test=="ks")
    ad <- sum(dist_test=="ad")
    cvm <- sum(dist_test=="cvm")
    dist_vec <- c(ks,ad,cvm)
  }

  if(dist_vec[1]==1){
    #ks_side
    if(ks_side!="less" & ks_side!="greater"){
      if(ks_side!="two"){message("The side of Kolmogorov-Smirnov test is wrong. It runs by default.")}
      ks_side <- "two.sided"
    }
  }

  # Estimation
  if(est_method=="MME"){
    #estimate lambda & phi
    est <- MME_WL(data)
    # lambda <- est[1]; phi <- est[2]

    #variance of lambda & phi
    cov <- MME_var(est[1], est[2],n)
    lambda_var <- cov[1,1]
    phi_var <- cov[2,2]

    #Confidence interval(Asymptatic or Bootstrap)
    if(CI_method == "boots"){
      store.boot = meantmpt = tmpt0 <- matrix(0,nrow=boots_iter,ncol=2)

      #Confidence interval(exponential or normal)
      for(i in 1:boots_iter){
        doubleboot.data <- sample(data,replace = T)
        tmp <- boot(doubleboot.data,function(y,indices){MME_WL(y[indices])},R=1000)
        meantmpt[i,] <- colMeans(tmp$t)
        tmpt0[i,] <- tmp$t0
      }
      if(CI_scale == "normal"){
        store.boot <- 2*tmpt0-meantmpt
        CI <- apply(store.boot,2,quantile,c(CI_alpha/2,1-CI_alpha/2))
        CI_lambda <- c(CI[1,1], CI[2,1])
        CI_phi <- c(CI[1,2], CI[2,2])

        if(CI[1,1]<0 | CI[1,2]<0){
          message("Warning: Confidence interval is out of parameter space. It runs by 'exp' CI_scale option.")
          CI_scale <- "exp"
        }
      }
      if(CI_scale == "exp"){
        store.boot <- 2*log(tmpt0)-log(meantmpt)
        CI <- exp( apply(store.boot,2,quantile,c(CI_alpha/2,1-CI_alpha/2) ))
        CI_lambda <- c(CI[1,1], CI[2,1])
        CI_phi <- c(CI[1,2], CI[2,2])
      }
    }
    else{
      d <- sqrt(diag(cov))

      if(CI_scale == "normal"){
        if(CI_side == "less"){
          CI_lower <- est-d*qnorm(1-CI_alpha)
          CI_upper <- est+d*qnorm(1-CI_alpha)
          CI_lambda <- c(CI_lower[1], CI_upper[1])
          CI_phi <- c(CI_lower[2], CI_upper[2])
        }else if(CI_side == "greater"){
          CI_lower <- est-d*qnorm(1-CI_alpha)
          CI_upper <- est+d*qnorm(1-CI_alpha)
          CI_lambda <- c(CI_lower[1], CI_upper[1])
          CI_phi <- c(CI_lower[2], CI_upper[2])
        }else{
          CI_lower <- est-d*qnorm(1-CI_alpha/2)
          CI_upper <- est+d*qnorm(1-CI_alpha/2)
          CI_lambda <- c(CI_lower[1], CI_upper[1])
          CI_phi <- c(CI_lower[2], CI_upper[2])
        }
        if(CI_lower[1]<0 | CI_lower[2]<0){
          message("Warning: Confidence interval is out of parameter space. It runs by 'exp' CI_scale option.")
          CI_scale <- "exp"
        }
      }
      if(CI_scale == "exp"){
        if(CI_side == "less"){
          CIexp_lower <- est*exp(-qnorm(1-CI_alpha)*d/est)
          CIexp_upper <- est*exp(+qnorm(1-CI_alpha)*d/est)
          CI_lambda <- c(CIexp_lower[1], CIexp_upper[1])
          CI_phi <- c(CIexp_lower[2], CIexp_upper[2])
        }else if(CI_side == "greater"){
          CIexp_lower <- est*exp(-qnorm(1-CI_alpha)*d/est)
          CIexp_upper <- est*exp(+qnorm(1-CI_alpha)*d/est)
          CI_lambda <- c(CIexp_lower[1], CIexp_upper[1])
          CI_phi <- c(CIexp_lower[2], CIexp_upper[2])
        }else{
          CIexp_lower <- est*exp(-qnorm(1-CI_alpha/2)*d/est)
          CIexp_upper <- est*exp(+qnorm(1-CI_alpha/2)*d/est)
          CI_lambda <- c(CIexp_lower[1], CIexp_upper[1])
          CI_phi <- c(CIexp_lower[2], CIexp_upper[2])
        }
      }
    }
  }
  else if(est_method=="MMEm"){
    #estimate lambda & phi
    est <- MMEm_WL(data) #lambda <- est[1]; phi <- est[2]

    #variance of lambda & phi
    cov <- MMEm_var(est[1], est[2],n)
    lambda_var <- cov[1,1]
    phi_var <- cov[2,2]

    #Confidence interval(Asymptatic or Bootstrap)
    if(CI_method == "boots"){
      store.boot = meantmpt = tmpt0 <- matrix(0,nrow=boots_iter,ncol=2)

      #Confidence interval(exponential or normal)
      for(i in 1:boots_iter){
        doubleboot.data <- sample(data,replace = T)
        tmp <- boot(doubleboot.data,function(y,indices){MMEm_WL(y[indices])},R=1000)
        meantmpt[i,] <- colMeans(tmp$t)
        tmpt0[i,] <- tmp$t0
      }
      if(CI_scale == "normal"){
        store.boot <- 2*tmpt0-meantmpt
        CI <- apply(store.boot,2,quantile,c(CI_alpha/2,1-CI_alpha/2))
        CI_lambda <- c(CI[1,1], CI[2,1])
        CI_phi <- c(CI[1,2], CI[2,2])

        if(CI[1,1]<0 | CI[1,2]<0){
          message("Warning: Confidence interval is out of parameter space. It runs by 'exp' CI_scale option.")
          CI_scale <- "exp"
        }
      }
      if(CI_scale == "exp"){
        store.boot <- 2*log(tmpt0)-log(meantmpt)
        CI <- exp( apply(store.boot,2,quantile,c(CI_alpha/2,1-CI_alpha/2) ))
        CI_lambda <- c(CI[1,1], CI[2,1])
        CI_phi <- c(CI[1,2], CI[2,2])
      }
    }
    else{
      d <- sqrt(diag(cov))

      #Confidence interval(exponential or normal)
      if(CI_scale == "normal"){
        if(CI_side == "less"){
          CI_lower <- est-d*qnorm(1-CI_alpha)
          CI_upper <- est+d*qnorm(1-CI_alpha)
          CI_lambda <- c(CI_lower[1], CI_upper[1])
          CI_phi <- c(CI_lower[2], CI_upper[2])
        }else if(CI_side == "greater"){
          CI_lower <- est-d*qnorm(1-CI_alpha)
          CI_upper <- est+d*qnorm(1-CI_alpha)
          CI_lambda <- c(CI_lower[1], CI_upper[1])
          CI_phi <- c(CI_lower[2], CI_upper[2])
        }else{
          CI_lower <- est-d*qnorm(1-CI_alpha/2)
          CI_upper <- est+d*qnorm(1-CI_alpha/2)
          CI_lambda <- c(CI_lower[1], CI_upper[1])
          CI_phi <- c(CI_lower[2], CI_upper[2])
        }

        if(CI_lower[1]<0 | CI_lower[2]<0){
          message("Warning: Confidence interval is out of parameter space. It runs by 'exp' CI_scale option.")
          CI_scale <- "exp"
        }
      }
      if(CI_scale == "exp"){
        if(CI_side == "less"){
          CIexp_lower <- est*exp(-qnorm(1-CI_alpha)*d/est)
          CIexp_upper <- est*exp(+qnorm(1-CI_alpha)*d/est)
          CI_lambda <- c(CIexp_lower[1], CIexp_upper[1])
          CI_phi <- c(CIexp_lower[2], CIexp_upper[2])
        }else if(CI_side == "greater"){
          CIexp_lower <- est*exp(-qnorm(1-CI_alpha)*d/est)
          CIexp_upper <- est*exp(+qnorm(1-CI_alpha)*d/est)
          CI_lambda <- c(CIexp_lower[1], CIexp_upper[1])
          CI_phi <- c(CIexp_lower[2], CIexp_upper[2])
        }else{
          CIexp_lower <- est*exp(-qnorm(1-CI_alpha/2)*d/est)
          CIexp_upper <- est*exp(+qnorm(1-CI_alpha/2)*d/est)
          CI_lambda <- c(CIexp_lower[1], CIexp_upper[1])
          CI_phi <- c(CIexp_lower[2], CIexp_upper[2])
        }
      }
    }
  }
  else if(est_method=="MLE"){
    #bias correction
    if(bias_cor=="coxsnell"){
      est <- MLE_WL(data, MME_WL(data)[2])
      est_temp <- est - CoxSnell_bias(n,est[1], est[2])
      if(est_temp[1]<0 | est_temp[2]<0){
        message("Warning: The estimators are out of parameter space. It runs by log scaled Cox &
                Snell bias correction option.")
        est_temp <- exp( log(est) - CoxSnell_bias_log(n,est[1], est[2]) )
      }
      est <- est_temp
      remove(est_temp)

      boot_tmp <- boot(data,function(data,indices){
        tmp = MLE_WL(data[indices], MME_WL(data[indices])[2])
        exp( log(tmp) - CoxSnell_bias_log(length(data[indices]),tmp[1], tmp[2]) ) },
        R=boots_iter)
      est_var <- apply(boot_tmp$t,2,var)
      lambda_var <- est_var[1]
      phi_var <- est_var[2]

      CI_method <- "boots"
      #CI_scale
      if(CI_method == "boots"){
        store.boot = meantmpt = tmpt0 <- matrix(0,nrow=boots_iter,ncol=2)

        #Confidence interval(exponential or normal)
        for(i in 1:boots_iter){
          doubleboot.data <- sample(data,replace = T)
          tmp <- boot(doubleboot.data,function(data,indices){
            tmp = MLE_WL(data[indices], MME_WL(data[indices])[2])
            tmp - CoxSnell_bias(length(data[indices]),tmp[1], tmp[2]) },R=1000)
          meantmpt[i,] <- colMeans(tmp$t)
          tmpt0[i,] <- tmp$t0
        }
        if(CI_scale == "normal"){
          store.boot <- 2*tmpt0-meantmpt
          CI <- apply(store.boot,2,quantile,c(CI_alpha/2,1-CI_alpha/2))
          CI_lambda <- c(CI[1,1], CI[2,1])
          CI_phi <- c(CI[1,2], CI[2,2])

          if(CI[1,1]<0 | CI[1,2]<0){
            message("Warning: Confidence interval is out of parameter space. It runs by 'exp' CI_scale option.")
            CI_scale <- "exp"
          }
        }
        if(CI_scale == "exp"){
          store.boot <- 2*log(tmpt0)-log(meantmpt)
          CI <- exp( apply(store.boot,2,quantile,c(CI_alpha/2,1-CI_alpha/2) ))
          CI_lambda <- c(CI[1,1], CI[2,1])
          CI_phi <- c(CI[1,2], CI[2,2])
        }
      }
    }
    else if(bias_cor=="firth"){
      est <- Firth_WL(data, MME_WL(data))
      if(est[1]<0 | est[2]<0){
        message("Warning: The estimators are out of parameter space. It runs by log scaled Firth method
                bias correction option.")
        est <- Firth_WL_log(data, MME_WL(data))
      }
      lambda <- est[1]; phi <- est[2]

      boot_tmp <- boot(data,function(data,indices){ Firth_WL(data[indices], MME_WL(data[indices])) },
                       R=boots_iter)
      est_var <- apply(boot_tmp$t,2,var)
      lambda_var <- est_var[1]
      phi_var <- est_var[2]

      CI_method <- "boots"
      #CI_scale
      if(CI_method == "boots"){
        store.boot = meantmpt = tmpt0 <- matrix(0,nrow=boots_iter,ncol=2)

        #Confidence interval(exponential or normal)
        for(i in 1:boots_iter){
          doubleboot.data <- sample(data,replace = T)
          tmp <- boot(doubleboot.data,function(data,indices){ Firth_WL(data[indices], MME_WL(data[indices])) },R=1000)
          meantmpt[i,] <- colMeans(tmp$t)
          tmpt0[i,] <- tmp$t0
        }
        if(CI_scale == "normal"){
          store.boot <- 2*tmpt0-meantmpt
          CI <- apply(store.boot,2,quantile,c(CI_alpha/2,1-CI_alpha/2))
          CI_lambda <- c(CI[1,1], CI[2,1])
          CI_phi <- c(CI[1,2], CI[2,2])

          if(CI[1,1]<0 | CI[1,2]<0){
            message("Warning: Confidence interval is out of parameter space. It runs by 'exp' CI_scale option.")
            CI_scale <- "exp"
          }
        }
        if(CI_scale == "exp"){
          store.boot <- 2*log(tmpt0)-log(meantmpt)
          CI <- exp( apply(store.boot,2,quantile,c(CI_alpha/2,1-CI_alpha/2) ))
          CI_lambda <- c(CI[1,1], CI[2,1])
          CI_phi <- c(CI[1,2], CI[2,2])
        }
      }
    }
    else{
      if(bias_cor!="None") message("The bias correction method is wrong. Proceed without bias correction.")
      bias_cor <- "None"

      est <- MLE_WL(data, MME_WL(data)[2])
      #lambda <- est[1]; phi <- est[2]

      #variance of lambda & phi
      cov <- MLE_var(est[1],est[2],n)
      lambda_var <- cov[1,1]
      phi_var <- cov[2,2]

      #Confidence interval(Asymptatic or Bootstrap)
      if(CI_method == "boots"){
        store.boot = meantmpt = tmpt0 <- matrix(0,nrow=boots_iter,ncol=2)

        #Confidence interval(exponential or normal)
        for(i in 1:boots_iter){
          doubleboot.data <- sample(data,replace = T)
          tmp <- boot(doubleboot.data,function(data,indices){
            MLE_WL(data[indices], MME_WL(data[indices])[2])},R=1000)
          meantmpt[i,] <- colMeans(tmp$t)
          tmpt0[i,] <- tmp$t0
        }
        if(CI_scale == "normal"){
          store.boot <- 2*tmpt0-meantmpt
          CI <- apply(store.boot,2,quantile,c(CI_alpha/2,1-CI_alpha/2))
          CI_lambda <- c(CI[1,1], CI[2,1])
          CI_phi <- c(CI[1,2], CI[2,2])

          if(CI[1,1]<0 | CI[1,2]<0){
            message("Warning: Confidence interval is out of parameter space. It runs by 'exp' CI_scale option.")
            CI_scale <- "exp"
          }
        }
        if(CI_scale == "exp"){
          store.boot <- 2*log(tmpt0)-log(meantmpt)
          CI <- exp( apply(store.boot,2,quantile,c(CI_alpha/2,1-CI_alpha/2) ))
          CI_lambda <- c(CI[1,1], CI[2,1])
          CI_phi <- c(CI[1,2], CI[2,2])
        }
      }
      else{
        d <- c(sqrt(lambda_var), sqrt(phi_var))

        #Confidence interval(exponential or normal)
        if(CI_scale == "normal"){
          if(CI_side == "less"){
            CI_lower <- est-d*qnorm(1-CI_alpha)
            CI_upper <- est+d*qnorm(1-CI_alpha)
            CI_lambda <- c(CI_lower[1], CI_upper[1])
            CI_phi <- c(CI_lower[2], CI_upper[2])
          }else if(CI_side == "greater"){
            CI_lower <- est-d*qnorm(1-CI_alpha)
            CI_upper <- est+d*qnorm(1-CI_alpha)
            CI_lambda <- c(CI_lower[1], CI_upper[1])
            CI_phi <- c(CI_lower[2], CI_upper[2])
          }else{
            CI_lower <- est-d*qnorm(1-CI_alpha/2)
            CI_upper <- est+d*qnorm(1-CI_alpha/2)
            CI_lambda <- c(CI_lower[1], CI_upper[1])
            CI_phi <- c(CI_lower[2], CI_upper[2])
          }

          if(CI_lower[1]<0 | CI_lower[2]<0){
            message("Warning: Confidence interval is out of parameter space. It runs by 'exp' CI_scale option.")
            CI_scale <- "exp"
          }
        }
        if(CI_scale == "exp"){
          if(CI_side == "less"){
            CIexp_lower <- est*exp(-qnorm(1-CI_alpha)*d/est)
            CIexp_upper <- est*exp(+qnorm(1-CI_alpha)*d/est)
            CI_lambda <- c(CIexp_lower[1], CIexp_upper[1])
            CI_phi <- c(CIexp_lower[2], CIexp_upper[2])
          }else if(CI_side == "greater"){
            CIexp_lower <- est*exp(-qnorm(1-CI_alpha)*d/est)
            CIexp_upper <- est*exp(+qnorm(1-CI_alpha)*d/est)
            CI_lambda <- c(CIexp_lower[1], CIexp_upper[1])
            CI_phi <- c(CIexp_lower[2], CIexp_upper[2])
          }else{
            CIexp_lower <- est*exp(-qnorm(1-CI_alpha/2)*d/est)
            CIexp_upper <- est*exp(+qnorm(1-CI_alpha/2)*d/est)
            CI_lambda <- c(CIexp_lower[1], CIexp_upper[1])
            CI_phi <- c(CIexp_lower[2], CIexp_upper[2])
          }
        }
      }
    }
  }
  else{
    if(est_method!="MLEc"){message("The estimation method is wrong. It runs by default.") }
    est_method <- "MLEc"
    #estimate lambda & phi
    if(bias_cor=="coxsnell"){
      est <- MLEc_WL(data)
      est_temp <- est - CoxSnell_bias(n,est[1], est[2])
      if(est_temp[1]<0 | est_temp[2]<0){
        message("Warning: The estimators are out of parameter space. It runs by log scaled Cox &
                Snell bias correction option.")
        est_temp <- exp( log(est) - CoxSnell_bias_log(n,est[1], est[2]) )
      }
      est <- est_temp

      #variance
      boot_tmp <- boot(data,function(data,indices){
        tmp = MLEc_WL(data[indices])
        exp( log(tmp) - CoxSnell_bias_log(length(data[indices]),tmp[1], tmp[2]) ) },R=boots_iter)
      est_var <- apply(boot_tmp$t,2,var)
      lambda_var <- est_var[1]
      phi_var <- est_var[2]

      CI_method <- "boots"
      #CI_scale
      if(CI_method == "boots"){
        store.boot = meantmpt = tmpt0 <- matrix(0,nrow=boots_iter,ncol=2)

        #Confidence interval(exponential or normal)
        for(i in 1:boots_iter){
          doubleboot.data <- sample(data,replace = T)
          tmp <- boot(doubleboot.data,function(data,indices){
            tmp = MLEc_WL(data[indices])
            tmp - CoxSnell_bias(length(data[indices]),tmp[1], tmp[2]) },R=1000)
          meantmpt[i,] <- colMeans(tmp$t)
          tmpt0[i,] <- tmp$t0
        }
        if(CI_scale == "normal"){
          store.boot <- 2*tmpt0-meantmpt
          CI <- apply(store.boot,2,quantile,c(CI_alpha/2,1-CI_alpha/2))
          CI_lambda <- c(CI[1,1], CI[2,1])
          CI_phi <- c(CI[1,2], CI[2,2])

          if(CI[1,1]<0 | CI[1,2]<0){
            message("Warning: Confidence interval is out of parameter space. It runs by 'exp' CI_scale option.")
            CI_scale <- "exp"
          }
        }
        if(CI_scale == "exp"){
          store.boot <- 2*log(tmpt0)-log(meantmpt)
          CI <- exp( apply(store.boot,2,quantile,c(CI_alpha/2,1-CI_alpha/2) ))
          CI_lambda <- c(CI[1,1], CI[2,1])
          CI_phi <- c(CI[1,2], CI[2,2])
        }
      }
    }
    else if(bias_cor=="boots"){

      boot_tmp <- boot(data,function(data,indices){
        MLEc_WL(data[indices]) }, R=boots_iter*10)
      est <- 2*boot_tmp$t0 - colMeans(boot_tmp$t) #2*original - boot_hat
      est <- matrix(est, nrow=1)

      est_var <- apply(boot_tmp$t,2,var)
      lambda_var <- est_var[1]
      phi_var <- est_var[2]

      CI_method <- "boots"
      #CI_scale
      if(CI_method == "boots"){
        store.boot = meantmpt = tmpt0 <- matrix(0,nrow=boots_iter,ncol=2)

        #Confidence interval(exponential or normal)
        for(i in 1:boots_iter){
          doubleboot.data <- sample(data,replace = T)
          tmp <- boot(doubleboot.data,function(data,indices){
            MLEc_WL(data[indices]) },R=1000)
          meantmpt[i,] <- colMeans(tmp$t)
          tmpt0[i,] <- tmp$t0
        }
        if(CI_scale == "normal"){
          store.boot <- 2*tmpt0-meantmpt
          CI <- apply(store.boot,2,quantile,c(CI_alpha/2,1-CI_alpha/2))
          CI_lambda <- c(CI[1,1], CI[2,1])
          CI_phi <- c(CI[1,2], CI[2,2])

          if(CI[1,1]<0 | CI[1,2]<0){
            message("Warning: Confidence interval is out of parameter space. It runs by 'exp' CI_scale option.")
            CI_scale <- "exp"
          }
        }
        if(CI_scale == "exp"){
          store.boot <- 2*log(tmpt0)-log(meantmpt)
          CI <- exp( apply(store.boot,2,quantile,c(CI_alpha/2,1-CI_alpha/2) ))
          CI_lambda <- c(CI[1,1], CI[2,1])
          CI_phi <- c(CI[1,2], CI[2,2])
        }
      }
    }
    else{
      if(bias_cor!="None") message("The bias correction method is wrong. Proceed without bias correction.")
      bias_cor <- "None"

      est <- MLEc_WL(data) #lambda <- est[1]; phi <- est[2]

      #variance of lambda & phi
      cov <- MLEc_var(est[1], est[2],n)
      lambda_var <- cov[1,1]
      phi_var <- cov[2,2]

      #Confidence interval(Asymptatic or Bootstrap)
      if(CI_method == "boots"){
        store.boot = meantmpt = tmpt0 <- matrix(0,nrow=boots_iter,ncol=2)

        #Confidence interval(exponential or normal)
        for(i in 1:boots_iter){
          doubleboot.data <- sample(data,replace = T)
          tmp <- boot(doubleboot.data,function(data,indices){MLEc_WL(data[indices])},R=1000)
          meantmpt[i,] <- colMeans(tmp$t)
          tmpt0[i,] <- tmp$t0
        }
        if(CI_scale == "normal"){
          store.boot <- 2*tmpt0-meantmpt
          CI <- apply(store.boot,2,quantile,c(CI_alpha/2,1-CI_alpha/2))
          CI_lambda <- c(CI[1,1], CI[2,1])
          CI_phi <- c(CI[1,2], CI[2,2])

          if(CI[1,1]<0 | CI[1,2]<0){
            message("Warning: Confidence interval is out of parameter space. It runs by 'exp' CI_scale option.")
            CI_scale <- "exp"
          }
        }
        if(CI_scale == "exp"){
          store.boot <- 2*log(tmpt0)-log(meantmpt)
          CI <- exp( apply(store.boot,2,quantile,c(CI_alpha/2,1-CI_alpha/2) ))
          CI_lambda <- c(CI[1,1], CI[2,1])
          CI_phi <- c(CI[1,2], CI[2,2])
        }
      }
      else{
        d <- sqrt(diag(cov))
        #d <- c(sqrt(lambda_var/length(x)), sqrt(phi_var/length(x)))

        #Confidence interval(exponential or normal)
        if(CI_scale == "normal"){
          if(CI_side == "less"){
            CI_lower <- est-d*qnorm(1-CI_alpha)
            CI_upper <- est+d*qnorm(1-CI_alpha)
            CI_lambda <- c(CI_lower[1], CI_upper[1])
            CI_phi <- c(CI_lower[2], CI_upper[2])
          }else if(CI_side == "greater"){
            CI_lower <- est-d*qnorm(1-CI_alpha)
            CI_upper <- est+d*qnorm(1-CI_alpha)
            CI_lambda <- c(CI_lower[1], CI_upper[1])
            CI_phi <- c(CI_lower[2], CI_upper[2])
          }else{
            if(CI_side != "two"){ message("The CI has wrong side. It runs by default.") }
            CI_side <- "two"
            CI_lower <- est-d*qnorm(1-CI_alpha/2)
            CI_upper <- est+d*qnorm(1-CI_alpha/2)
            CI_lambda <- c(CI_lower[1], CI_upper[1])
            CI_phi <- c(CI_lower[2], CI_upper[2])
          }

          if(CI_lower[1]<0 | CI_lower[2]<0){
            message("Warning: Confidence interval is out of parameter space. It runs by 'exp' CI_scale option.")
            CI_scale <- "exp"
          }
        }
        if(CI_scale == "exp"){
          if(CI_side == "less"){
            CIexp_lower <- est*exp(-qnorm(1-CI_alpha)*d/est)
            CIexp_upper <- est*exp(+qnorm(1-CI_alpha)*d/est)
            CI_lambda <- c(CIexp_lower[1], CIexp_upper[1])
            CI_phi <- c(CIexp_lower[2], CIexp_upper[2])
          }else if(CI_side == "greater"){
            CIexp_lower <- est*exp(-qnorm(1-CI_alpha)*d/est)
            CIexp_upper <- est*exp(+qnorm(1-CI_alpha)*d/est)
            CI_lambda <- c(CIexp_lower[1], CIexp_upper[1])
            CI_phi <- c(CIexp_lower[2], CIexp_upper[2])
          }else{
            CIexp_lower <- est*exp(-qnorm(1-CI_alpha/2)*d/est)
            CIexp_upper <- est*exp(+qnorm(1-CI_alpha/2)*d/est)
            CI_lambda <- c(CIexp_lower[1], CIexp_upper[1])
            CI_phi <- c(CIexp_lower[2], CIexp_upper[2])
          }
        }
      }
    }
  }

  # Goodness of fit test
  test_list <- list("dist_vec"=dist_vec, "gof_alpha"=gof_alpha)
  ks_list = ad_list = cvm_list <- NA
  if(dist_vec[1]==1){
    ks <- ks.test(data,pwlindley,theta=est[1],alpha=est[2],alternative=ks_side)
    ks_list <- list("ks_side"=ks_side, "ks_stat"=as.numeric(ks$statistic), "ks_pvalue"=ks$p.value)
  }
  if(dist_vec[2]==1){
    ad <- ad.test(data,null = "pwlindley",theta = est[1],alpha =est[2])
    ad_list <- list("ad_stat"=as.numeric(ad$statistic), "ad_pvalue"=ad$p.value)
  }
  if(dist_vec[3]==1){
    cvm <- cvm.test(data,null = "pwlindley",theta = est[1],alpha =est[2])
    cvm_list <- list("cvm_stat"=as.numeric(cvm$statistic), "cvm_pvalue"=cvm$p.value)
  }

  # wilks test
  if(wilks_test == TRUE){
    # wilks tset side
    if(wilks_side!="less" & wilks_side!="greater"){
      if(wilks_side!="two"){message("The side of Wilks' theorem test is wrong. It runs by default.")}
      wilks_side <- "two"
    }

    wilks_temp <- wilks.test(data, est, side=wilks_side)
    wilks_list <- list("wilks_test"=wilks_test, "alpha"=wilks_alpha)
    wilks_list <- c(wilks_list, wilks_temp)
  }else{ wilks_list <- list("wilks_test"=wilks_test) }

  # making list
  test_list <- c(test_list,ks_list,ad_list,cvm_list)
  test_list <- test_list[!sapply(test_list, function(x) all(is.na(x)))]
  CI_list <- list("CI_lambda"=CI_lambda, "CI_phi"=CI_phi, "CI_side"=CI_side,
                  "CI_per"=CI_per, "CI_method"=CI_method, "CI_scale"=CI_scale)

  result_list <- list("data"=data, "dataname"=dataname,"stat_summary"=stat_summary, "mean"=meanx, "var"=varx,
                      "est"=est, "lambda_var"=lambda_var, "phi_var"=phi_var,
                      "bias_cor"=bias_cor, "est_method"=est_method, "boots_iter"=boots_iter,
                      "test_list"=test_list, "CI_list"=CI_list, "wilks_list"=wilks_list)

  class(result_list)<- "WL"
  return(result_list)
}

#' @method print WL
#' @export
print.WL <- function(x, digits = max(3, getOption("digits") - 3), ...){
  dataname <- x$dataname
  a0 <- paste0("\nData: ",dataname,"\n")
  cat(a0)

  cat("\nEstimate:\n")
  lambda <- x$est[1]; phi <- x$est[2]
  lambda <- round(lambda, digits=digits); phi <- round(phi, digits=digits)
  a1 <- paste0("lambda= ",lambda,", ","phi= ",phi)
  cat(a1)
  cat("\n")
}

#' Summarizing WL function
#'
#' \code{summary} method for a class "WL".
#'
#' @param object an object of class "WL" made by the function \code{WL}.
#' @param ... not used, but exists because of the compatibility.
#' @param x an object of class "summary.WL".
#' @param digits a numeric number of significant digits.
#'
#' @method summary WL
#'
#' @examples
#' example <- fail_fiber
#' result <- WL(example, est_method="MLEc")
#' summary(result)
#'
#' @export
summary.WL <- function(object,...){
  stat_summary <- object$stat_summary
  colnames(stat_summary) <- c("Min","1st Qu","Median","3rd Qu","Max")

  est_method <- object$est_method
  if(est_method == "MMEm"){est_method <- "Modified MME"}
  if(est_method=="MLE" | est_method=="MLEc"){
    bias_cor <- object$bias_cor
    if(bias_cor=="coxsnell"){bias_cor <- "Cox & Snell"}
    else if(bias_cor=="boots"){bias_cor <- "Bootstrap"}
    else if(bias_cor=="firth"){bias_cor <- "Firth's"}
    else(bias_cor <- "None")
  }else{bias_cor <- "None"}

  dist_vec <- object$test_list$dist_vec
  ks_side <- NA; ks_stat <- NA; ks_pvalue <- NA
  ad_stat <- NA; ad_pvalue <- NA
  cvm_stat <- NA; cvm_pvalue <- NA

  if(dist_vec[1]==1){
    ks_side <- object$test_list$ks_side
    if(ks_side=="two.sided"){ks_side <- "Two-sided"}
    else if(ks_side=="less"){ks_side <- "One-sided(less)"}
    else if(ks_side=="greater"){ks_side <- "One-sided(greater)"}

    ks_side <- ks_side
    ks_stat <- object$test_list$ks_stat
    ks_pvalue <- object$test_list$ks_pvalue
  }
  if(dist_vec[2]==1){
    ad_stat <- object$test_list$ad_stat
    ad_pvalue <- object$test_list$ad_pvalue
  }
  if(dist_vec[3]==1){
    cvm_stat <- object$test_list$ad_stat
    cvm_pvalue <- object$test_list$ad_pvalue
  }

  CI_side <- object$CI_list$CI_side
  if(CI_side=="two"){CI_side <- "Two-sided"}
  else if(CI_side=="less"){CI_side <- "One-sided(less)"}
  else if(CI_side=="greater"){CI_side <- "One-sided(greater)"}

  CI_method <- object$CI_list$CI_method
  if(CI_method=="asymp"){CI_method <- "Asymptotic"}
  else if(CI_method=="boots"){CI_method <- "Bootstrap"}

  CI_scale <- object$CI_list$CI_scale
  if(CI_scale=="exp"){CI_scale <- "Exponential"}
  else if(CI_scale=="normal"){CI_scale <- "Normal"}

  wilks_stat=wilks_alpha=wilks_pvalue=wilks_result <- NA
  if(object$wilks_list$wilks_test == TRUE){
    wilks_side <- object$wilks_list$side
    wilks_stat <- object$wilks_list$stat
    wilks_alpha <- object$wilks_list$alpha
    wilks_pvalue <- object$wilks_list$pvalue

    if(wilks_side=="two"){wilks_side <- "Two-sided"}
    else if(wilks_side=="less"){wilks_side <- "One-sided(less)"}
    else if(wilks_side=="greater"){wilks_side <- "One-sided(greater)"}

    if(wilks_pvalue > wilks_alpha){
      wilks_result <- paste0("The null hypothesis cannot be rejected.")
    }else{
      wilks_result <- paste0("The null hypothesis can be rejected.")
    }
  }


  result <- list("stat_summary"=stat_summary, "mean"=object$mean, "var"=object$var, "lambda"=object$est[1],
                 "phi"=object$est[2], "lambda_var"=object$lambda_var, "phi_var"=object$phi_var,
                 "dataname"=object$dataname, "est_method"=est_method, "bias_cor"=bias_cor,
                 "boots_iter"=object$boots_iter,
                 "dist_vec"=dist_vec, "gof_alpha"=object$test_list$gof_alpha,
                 "ks_side"=ks_side, "ks_stat"=ks_stat, "ks_pvalue"=ks_pvalue,
                 "ad_stat"=ad_stat, "ad_pvalue"=ad_pvalue, "cvm_stat"=cvm_stat, "cvm_pvalue"=cvm_pvalue,
                 "CI_lambda"=object$CI_list$CI_lambda, "CI_phi"=object$CI_list$CI_phi, "CI_side"=CI_side,
                 "CI_per"=object$CI_list$CI_per, "CI_method"=CI_method, "CI_scale"=CI_scale,
                 "wilks_test"=object$wilks_list$wilks_test, "wilks_stat"=wilks_stat, "wilks_side"=wilks_side,
                 "wilks_pvalue"=wilks_pvalue, "wilks_alpha"=wilks_alpha, "wilks_result"=wilks_result)
  class(result) <- "summary.WL"
  result
}

#' @rdname summary.WL
#' @method print summary.WL
#' @export
print.summary.WL <- function(x, digits = max(3, getOption("digits") - 3), ...){
  dataname <- x$dataname
  cat("\nData:", dataname, "\n")
  cat("\n")
  cat("Data summary:\n")

  mean <- round(x$mean, digits=digits)
  var <- round(x$var, digits=digits)
  a1 <- paste0("  Mean: ",mean,", Variance: ",var)
  cat(a1, "\n")

  a2 <- x$stat_summary
  rownames(a2) <- rep("", nrow(a2))
  print(a2, quote=FALSE, right=TRUE, digits = digits)

  alpha <- x$gof_alpha
  if(x$dist_vec[1]==1){
    a3 <- paste0("\n",x$test_side,"Kolmogorov-Smirnov test for ","alpha=",alpha)
    cat(a3,"\n")
    s <- round(x$ks_stat, digits=digits);
    p <- round(x$ks_pvalue, digits=digits);
    a4 <- paste0("  D = ",s,", p-value: ",p)
    cat(a4)
    cat("\n")
    if(p > alpha){
      a5 <- paste0("  >> This data follows the weighted lindley distribution with estimated parameters.")
      cat(a5)
    }else{
      a5 <- paste0("  >> This data does not follow the weighted lindley distribution
                   with estimated parameters.")
      cat(a5)
    }
    cat("\n")
  }
  if(x$dist_vec[2]==1){
    a3 <- paste0("\n",x$test_side,"Anderson-Darling test for ","alpha=",alpha)
    cat(a3,"\n")
    s <- round(x$ad_stat, digits=digits);
    p <- round(x$ad_pvalue, digits=digits);
    a4 <- paste0("  A = ",s,", p-value: ",p)
    cat(a4)
    cat("\n")
    if(p > alpha){
      a5 <- paste0("  >> This data follows the weighted lindley distribution with estimated parameters.")
      cat(a5)
    }else{
      a5 <- paste0("  >> This data does not follow the weighted lindley distribution
                   with estimated parameters.")
      cat(a5)
    }
    cat("\n")
  }
  if(x$dist_vec[3]==1){
    a3 <- paste0("\n",x$test_side,"Cramer-von Mises test for ","alpha=",alpha)
    cat(a3,"\n")
    s <- round(x$cvm_stat, digits=digits);
    p <- round(x$cvm_pvalue, digits=digits);
    a4 <- paste0("  T = ",s,", p-value: ",p)
    cat(a4)
    cat("\n")
    if(p > alpha){
      a5 <- paste0("  >> This data follows the weighted lindley distribution with estimated parameters.")
      cat(a5)
    }else{
      a5 <- paste0("  >> This data does not follow the weighted lindley distribution
                   with estimated parameters.")
      cat(a5)
    }
    cat("\n")
  }

  if(x$bias_cor=="Bootstrap"){
    cat("\nEstimation method(bias correction): ")
    a6 <- paste0(x$est_method,"(",x$bias_cor,"with ",x$boots_iter," times)")
    cat(a6,"\n")
    lambda <- round(x$lambda, digits=digits); phi <- round(x$phi, digits=digits)
    a7 <- paste0("  lambda: ",lambda,", phi: ",phi)
    cat(a7,"\n")
  }else{
    cat("\nEstimation method(bias correction): ")
    a6 <- paste0(x$est_method,"(",x$bias_cor,")")
    cat(a6,"\n")
    lambda <- round(x$lambda, digits=digits); phi <- round(x$phi, digits=digits)
    a7 <- paste0("  lambda: ",lambda,", phi: ",phi)
    cat(a7,"\n")
  }

  cat("\nVariance of lambda & phi:\n")
  var1 <- round(x$lambda_var, digits=digits)
  var2 <- round(x$phi_var, digits=digits)
  a8 <- paste0("  Var(lambda) = ", var1,", Var(phi) = ", var2)
  cat(a8, "\n")
  cat("\n")

  if(x$CI_method == "Bootstrap"){
    a9 <- paste0(x$CI_side," confidence intervals for ",x$CI_per,"%")
    cat(a9,"\n")
    a10 <- paste0("(Bootstrap method with ", x$boots_iter*10, " times & ", x$CI_scale," scaled )")
    cat(a10,"\n")
  }else{
    a9 <- paste0(x$CI_side," confidence interval for ",x$CI_per,"%")
    cat(a9,"\n")
    a10 <- paste0("(Asymptotic method & ", x$CI_scale," scaled )")
    cat(a10,"\n")
  }
  CI_lambda <- round(x$CI_lambda, digits=digits)
  a11 <- paste0("  CI for lambda: ","(",CI_lambda[1],", ",CI_lambda[2],")")
  cat(a11,"\n")
  CI_phi <- round(x$CI_phi, digits=digits)
  a12 <- paste0("  CI for phi: ","(",CI_phi[1],", ",CI_phi[2],")")
  cat(a12,"\n")

  if(x$wilks_test==TRUE){
    a13 <- paste0("\n",x$wilks_side," Wilks' theorem test for estimated parameters")
    cat(a13,"\n")
    s <- round(x$wilks_stat, digits=digits);
    p <- round(x$wilks_pvalue, digits=digits);
    a14 <- paste0("  X = ",s,", p-value: ",p)
    cat(a14)
    cat("\n")
    a15 <- paste0("  >>",x$wilks_result)
    cat(a15)
    cat("\n")
  }
  cat("\n")
}

#' Some plots for \code{WL}
#'
#' \code{plot} method for a class "WL".
#'
#' @method plot WL
#'
#' @param x an object of class "WL" made by the function \code{WL}.
#' @param which  if a subset of the plots is required, specify a subset of 1:4; see 'Details' for a
#' description of the plots.
#' @param ask logical; if TRUE, the user is asked before each plot.
#' @param ... other parameters to be passed through to plotting functions.
#'
#' @details The first plot (\code{which=1}) is histogram of input values and pdf from estimated lambda and phi.
#' This is a good plot to judge how well the estimated values match the actual data.
#' \cr
#' \cr
#' The second plot (\code{which=2}) is the boxplot of the input data. It shows the simple information about median and outlier
#' points. Additionally, you can easily see which the "Outlier Points" are with box on the right side.\cr
#' \cr
#' The third plot (\code{which=3}) is Q-Q plot. This provides the same information as the normal Q-Q plot.\cr
#' \cr
#' The last plot (\code{which=4}) is contour plot of lambda and phi. You can observe the location of estimator, and also
#' compare estimator to Cox&Snell bias corrected MLE and MLEc.\cr
#' \cr
#' Plots are basically displayed one by one in the order mentioned above,
#' and it is also possible to view only one plot or four at a time.
#' How to do this is shown in the example below.
#'
#'
#' @seealso See also \code{\link{WL}} and \code{\link{plot}}.
#'
#' @examples
#' example <- lifetime_alum
#' result <- WL(example)
#'
#' plot(result)
#' plot(result, which=1)
#' plot(result, which=3)
#' par(mfrow = c(2, 2)) ; plot(result)
#' par(mfrow = c(1, 1))
#'
#' @export
plot.WL <- function(x, which=c(1,2,3,4),
                    ask = prod(par("mfcol")) < length(which) && dev.interactive(), ...){
  #Histogram and density function
  if(class(x) != "WL") stop("object not of class 'WL'")
  show <- rep(FALSE, 4)
  show[which] <- TRUE
  dots <- list(...)
  nmdots <- names(dots)

  one.fig <- prod(par("mfcol")) == 1
  if (ask) {
    oask <- devAskNewPage(TRUE)
    on.exit(devAskNewPage(oask))
  }

  #first plot
  if(show[1]==TRUE){
    h <- hist(x$data, plot=FALSE)
    plot_main <- paste0("Histogram of ", x$dataname, " and density function")
    xlines <-seq(min(h$breaks),max(h$breaks),length.out=100)
    lambda <- x$est[1]
    phi <- x$est[2]
    y <- dwlindley(xlines,theta=lambda, alpha=phi)
    max_y <- max(y*length(x$data)*diff(h$breaks)[1])
    hist(x$data, col="grey", main=plot_main, xlab=x$dataname, ylim=c(0,max_y))
    lines(x=xlines, y*length(x$data)*diff(h$breaks)[1])
    component1 <- x$est_method
    component2 <- paste0(intToUtf8(955),"=",round(lambda,2))
    component3 <- paste0(intToUtf8(966),"=",round(phi,2))
    component <- c(component1,component2,component3)
    legend("topright", legend=component)
  }

  #second plot
  if(show[2]==TRUE){
    box_main <- paste0("Boxplot of ", x$dataname, " and its outliers")
    quant_vec <- quantile(x$data)
    upper_line <- quant_vec[4] + 1.5*(quant_vec[4] - quant_vec[2])
    lower_line <- quant_vec[2] - 1.5*(quant_vec[4] - quant_vec[2])
    outliers <- c(which(x$data > upper_line), which(x$data < lower_line))
    if(length(outliers)==0){
      label_out <- paste0("outlier does not exist.")
      boxplot(x$data, main=box_main ,xlab=x$dataname)
      legend("topright", legend=label_out, text.font=14)
    }else{
      label_out <- paste0("   ", outliers)
      boxplot(x$data, main=box_main ,xlab=x$dataname)
      text(x=1.22, y=x$data[outliers], labels=outliers)
      legend("topright", legend=c("Outlier points",label_out), text.font=14)
    }
  }

  #third plot
  if(show[3]==TRUE){
    data_sort <- sort(x$data)
    x0 <- qwlindley(ppoints(length(data_sort)),x$est[1],x$est[2])
    plot(x = x0, y = data_sort, main="Weigthed Lindley Q-Q Plot",
         xlab = "Theoretical quantiles", ylab = "Observed quantiles");
    abline(a = 0, b = 1, col = "red");

    component1 <- x$est_method
    component2 <- paste0(intToUtf8(955),"=",round(x$est[1],2))
    component3 <- paste0(intToUtf8(966),"=",round(x$est[2],2))
    component <- c(component1,component2,component3)
    legend("topleft", legend=component)
  }

  #fourth plot
  if(show[4]==TRUE){
    n <- length(x$data)
    est <- MLE_WL(x$data, MME_WL(x$data)[1])
    cmle <- est - CoxSnell_bias(n,est[1], est[2])
    if(cmle[1]<0 | cmle[2]<0){
      cmle <- exp( log(est) - CoxSnell_bias_log(n,est[1], est[2]) )
    }

    est <- MLEc_WL(x$data)
    cmlec <- est - CoxSnell_bias(n,est[1], est[2])
    if(cmlec[1]<0 | cmlec[2]<0){
      cmlec <- exp( log(est) - CoxSnell_bias_log(n,est[1], est[2]) )
    }

    n1 <- 100
    lambda_vec <- seq(x$CI_list$CI_lambda[1], x$CI_list$CI_lambda[2], length.out=n1)
    phi_vec <- seq(x$CI_list$CI_phi[1], x$CI_list$CI_phi[2], length.out=n1)

    z <- matrix(0, nrow=n1, ncol=n1)
    for(i in 1:n1){
      lambda <- lambda_vec[i]
      for(j in 1:n1){
        phi <- phi_vec[j]
        z[i,j] <- sum( log( dwlindley(x$data,lambda,phi) ) )
      }
    }
    contour_main <- paste0("Contour plot of log-likelihood function")
    contour(lambda_vec,phi_vec,z, xlab="lambda", ylab="phi", main=contour_main)
    points(cmle[1], cmle[2], pch=5,col="black");
    points(cmlec[1], cmlec[2], pch=4,col="black");
    points(x$est[1], x$est[2], pch=19,col="red")
    legend("topright", legend=c("CMLE","CMLEc",x$est_method), col=c("black","black","red"),
           pch=c(5,4,19))#, text.font=14)
  }
  invisible()
}

#' The failure stresses (in GPA) of single carbon fibers
#'
#' The data represents the strength measured in GPA for single carbon fibers and impregnated
#' 1000-carbon fiber tows.
#' Single fibers were tested under tension at gauge lengths of 50 mm,
#' with sample sizes n = 65 (Bader et al., 1982).
#'
#' @usage data(fail_fiber, package = "WL")
#' @format A dataframe with 65 observations.
#'
#' @references Bader, M. G., Priest, A. M. (1982). Statistical aspects of fibre and bundle strength in
#' hybrid composites. In: Hayashi, T., Kawata, K., Umekawa, S., eds. Progress in Science and
#' Engineering Composites: Proceedings of the Fourth International Conference on Composite Materials,
#' ICCM-IV. Tokyo: Japan Society for Composite Materials.
#' @references Hyoung-Moon Kim. and Yu-Hyeong Jang. (2020). New Closed-Form Estimators for Weighted Lindley
#' Distribution. \emph{ }, submitted.
"fail_fiber"

#' The lifetimes of aluminum specimens
#'
#' This data set represents the lifetimes of 102 aluminum specimens exposed to
#' a maximum stress/cycle of 26000 PSI (PSI is a unit of pressure expressed in
#' pounds of force per square inch of area; it stands for Pounds per Square Inch.
#' 1 PSI = 6894 Pascals = 0.070 atmospheres = 51.715 torr) (Owen et al., 2000).
#'
#' @usage data(lifetime_alum, package = "WL")
#' @format A dataframe with 102 observations.
#'
#' @references Owen W.J., Padget W.J. (2000). A Birnbaum-Saunders accelerated life model.
#' \emph{IEEE Trans Reliability}, 49(2):224-229.
#' @references Hyoung-Moon Kim. and Yu-Hyeong Jang. (2020). New Closed-Form Estimators for Weighted Lindley
#' Distribution. \emph{ }, submitted.
"lifetime_alum"

#' The lifetime failure data of an electronic device
#'
#' It represents the lifetime failure data of an electronic device (Wang, 2000).
#'
#' @usage data(fail_time, package = "WL")
#' @format A dataframe with 18 observations.
#'
#' @references Wang, F. (2000). A new model with bathtub-shaped failure rate using an additive Burr XII
#' distribution. \emph{Reliability Engineering and System Safety}, 70:305-312.
#' @references Hyoung-Moon Kim. and Yu-Hyeong Jang. (2020). New Closed-Form Estimators for Weighted Lindley
#' Distribution. \emph{ }, submitted.
"fail_time"

