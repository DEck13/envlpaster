#####################################################
###                                               ###
###  Functions for the 1D algortihm               ###
###                                               ###
#####################################################


###################################################
#                projection                       #
###################################################
projection <- function(a){
  d <- dim(a)[2]
  if(sum(t(a)%*%a)==0){
    return(0)
  }
  pa <- a%*%matpower(t(a)%*%a,-1)%*%t(a)
  return(pa)
}

##################################################
#                  matrix power                  #
##################################################
matpower <- function(a,alpha){
  small <- 0.000001
  p1<-nrow(a)
  eva<-eigen(a)$values
  eve<-as.matrix(eigen(a)$vectors)
  eve<-eve/t(matrix((diag(t(eve)%*%eve)^0.5),p1,p1))
  index<-(1:p1)[eva>small]
  evai<-eva
  evai[index]<-(eva[index])^(alpha)
  foo <- NULL
  if(length(evai) == 1) foo <- diag(evai, nrow = 1)
  else foo <- diag(evai)
  ai<-as.matrix(eve)%*%foo%*%t(as.matrix(eve))
  return(ai)
}

##################################################
#         1D objective function                  #
##################################################
get1Dobj <- function(w,A,B){
  small <- 0.000001
  p <- dim(A)[1]
  foo <- eigen((A+B), symmetric = TRUE)
  if(p == 1) B.int <- foo$vec %*% 1/foo$val %*% t(foo$vec)
  else B.int <- foo$vec %*% diag(1/foo$val) %*% t(foo$vec)
  Fw <- log(t(w)%*%A%*%w + small) + log(t(w)%*%B.int%*%w + small) - 2*log(t(w)%*%w)
  return(Fw)
}

##################################################
#    get initial value for 1D algorithm          #
##################################################
get1Dini <- function(A,B){
  p <- dim(A)[1]
  vecs <- cbind(eigen(A, symmetric = TRUE)$vectors,
                eigen(A+B, symmetric = TRUE)$vectors)
  idx <- order(apply(vecs,2,get1Dobj,A,B))[1]
  w <- vecs[,idx]
  return(w)
}

##################################################
#         1D objective function gradient         #
##################################################
get1Dderiv <- function(w,A,B){
  p <- dim(A)[1]
  foo <- eigen((A + B), symmetric = TRUE)
  if(p == 1) B.int <- foo$vec %*% 1/foo$val %*% t(foo$vec)
  else B.int <- foo$vec %*% diag(1/foo$val) %*% t(foo$vec)
  dF <- c(2/(t(w)%*%A%*%w))*A%*%w + c(2/(t(w)%*%B.int%*%w))*B.int%*%w - c(4/(t(w)%*%w))*w
  return(dF)
}

##################################################
# 1D manifold algorithm for u-dim envelope       #
# the algorithm based on M and inv(M+U)          #
# and Polak-Ribiere conjugate gradient (PRCG)    #
##################################################
#' manifold1Dplus
#' @description
#' \loadmathjax
#' The 1D algorithm
#' @details This function calls \code{get1Dobj}, \code{get1Dini}, and \code{get1Dderiv}
#' in order to find \mjdeqn{ \max_{w} \left[ \log(w^TMw) + \log(w^T(M+U)w) - 2\log(w^Tw) \right] }{ascii}
#' using Polak-Ribiere conjugate gradient in \code{optim}.
#' This maximization is conducted a total of \code{u} times and at each iteration
#' a vector belonging to the envelope space is returned.
#' The vector returned at a specific iteration is orthogonal to the vectors
#' returned at previous iterations. When finished, a basis matrix for the envelope space is returned.
#' @param M A \mjeqn{\sqrt{n}}{ascii} estimate of an estimator's asymptotic covariance matrix.
#' @param U A \mjeqn{\sqrt{n}}{ascii} estimate of the parameter associated with the space we are enveloping. 
#' For our purposes this quantity is either the outer product of the MLE of the mean-value 
#' submodel parameter vector with itself or the outer product of the 
#' MLE of the canonical submodel parameter vector with itself.
#' @param u The dimension of the envelope space assumed.
#' @return \item{G}{A \mjeqn{\sqrt{n}}{ascii} estimator of the basis matrix for the
#' envelope subspace. This matrix has \code{u} columns}
#' @references 
#' 
#' Cook, R.D. and Zhang, X. (2015 a). Foundations for Envelope Models and Methods. 
#' \emph{Journal of the American Statistical Association}, \strong{110}, 599-611. \cr
#' \cr
#' Cook, R.D. and Zhang, X. (2015 b). Algorithms for Envelope Estimation. 
#' \emph{Journal of Computational and Graphical Statistics}, Published online. \doi{10.1080/10618600.2015.1029577}
#' @examples \dontrun{library(envlpaster)
#'  data(simdata30nodes)
#'  data <- simdata30nodes.asterdata
#'  nnode <- length(vars)
#'  xnew <- as.matrix(simdata30nodes[,c(1:nnode)])
#'  m1 <- aster(xnew, root, pred, fam, modmat)
#'  avar <- m1$fisher
#'  beta <- m1$coef
#'  U <- beta \%o\% beta
#'  manifold1Dplus(M = avar, U = U, u = 1)}
#' @export
manifold1Dplus <- function(M,U,u){
  p <- dim(M)[1]
  Mnew <- M
  Unew <- U
  G <- matrix(0,p,u)
  G0 <- diag(1,p)
  for(i in 1:u){  # used to be A = solve(Mnew + Unew), B = Mnew
    ans <- optim(get1Dini(Mnew,Unew),get1Dobj,get1Dderiv,
                 A=Mnew,B=Unew,method="CG",
                 control=list(maxit=500,type=2))
    w <- c(ans$par)
    gk <- c(1/sqrt(sum(w^2)))*w
    if(p == 1) G[,i] <- G0 * gk
    else G[,i] <- G0%*%gk
    G0 <- qr.Q(qr(G[,1:i]),complete=T)
    G0 <- G0[,(i+1):p]
    Mnew <- t(G0)%*%M%*%G0
    Unew <- t(G0)%*%U%*%G0
  }
  return(G)
}

