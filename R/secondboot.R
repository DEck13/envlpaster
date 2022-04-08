#' secondboot
#' @description
#' \loadmathjax
#'   A parametric bootstrap procedure evaluated at an envelope estimator of the submodel mean-value parameter vector \mjeqn{\tau}{ascii} that was obtained using eigenstructures or the 1d algorithm.
#' @details
#' This function implements the second level of the parametric bootstrap
#' procedure given by either Algorithm 1 or Algorithm 2 in Eck (2015) with
#' respect to the mean-value parameterization. This is detailed in Steps 4
#' through 5c in the algorithm below. At iteration \mjeqn{b}{ascii}, this parametric
#' bootstrap generates resamples from the distribution evaluated at the
#' envelope estimator (\mjeqn{\hat{\tau}_{env}^{(b)}}{ascii}) of \mjeqn{\tau}{ascii}.
#' In this case, the selected indices producing the eigenstructure which was used to
#' construct the envelope estimator \mjeqn{\hat{\tau}_{env}^{(b)}}{ascii} are used to
#' construct envelope estimators for the generated data. These resamples
#' are used to estimate the variability of \mjeqn{\hat{\tau}_{env}^{(b)}}{ascii}.
#' The algorithm using eigenstructures is as follows:
#'
#'   \enumerate{
#'     \item Fit the aster model to the data and obtain
#'     \mjeqn{\hat{\tau} = (\hat{\gamma}^T, \hat{\upsilon}^T)}{ascii} and \mjeqn{\hat{\Sigma}}{ascii}
#'     from the aster model fit.
#'     \item Compute the envelope estimator of \mjeqn{\upsilon}{ascii} in the original sample, given as
#'     \mjeqn{\hat{\upsilon}_{env} = P_{\hat{G}}\hat{\upsilon}}{ascii} where \mjeqn{P_{\hat{G}}}{ascii} is computed using eigenstructures
#'     and selected via a model selection criterion of choice.
#'     \item Perform a parametric bootstrap by generating resamples from the distribution of the
#'     aster submodel evaluated at \mjeqn{\hat{\tau}_{env} = (\hat{\gamma}^T,\hat{\upsilon}_{env}^T)^T}{ascii}. For iteration
#'     \mjeqn{b=1,...,B}{ascii} of the procedure:
#'       \enumerate{
#'         \item Compute \mjeqn{\hat{\tau}^{(b)}}{ascii} and \mjeqn{\widehat{\Sigma}_{\upsilon,\upsilon}^{(b)}}{ascii} from the aster
#'         model fit to the resampled data.
#'         \item Build \mjeqn{P_{\hat{G}}^{(b)}}{ascii} using the indices of \mjeqn{\hat{\Sigma}_{\upsilon,\upsilon}^{(b)}}{ascii} that
#'         are selected using the same model selection criterion as Step 2 to build \mjeqn{\hat{G}}{ascii}.
#'         \item Compute \mjeqn{\hat{\upsilon}_{env}^{(b)} = P_{\hat{\mathcal{E}}}^{(b)}\hat{\upsilon}^{(b)}}{ascii}
#'         and
#'         \item Store \mjeqn{\hat{\tau}_{env}^{(b)}}{ascii} and \mjeqn{g\left(\hat{\tau}_{env}^{(b)}\right)}{ascii}
#'         where \mjeqn{g}{ascii} maps \mjeqn{\tau}{ascii} to the parameterization of Darwinian fitness.
#'       }
#'     \item After \mjeqn{B}{ascii} steps, the bootstrap estimator of expected Darwinian fitness is the
#'     average of the envelope estimators stored in Step 3d. This completes the first part of the
#'     bootstrap procedure.
#'     \item We now proceed with the second level of bootstrapping at the \mjeqn{b^{th}}{ascii} stored envelope estimator
#'     \mjeqn{\hat{\tau}_{env}^{(b)}}{ascii}. For iteration \mjeqn{k=1,...,K}{ascii} of the procedure:
#'       \enumerate{
#'         \item Generate data from the distribution of the aster submodel evaluated at \mjeqn{\hat{\tau}_{env}^{(b)}}{ascii}.
#'         \item Perform Steps 3a through 3d with respect to the dataset obtained in Step 5a.
#'         \item Store \mjeqn{\hat{\tau}_{env}^{(b)^{(k)}}}{ascii} and \mjeqn{g\left(\hat{\tau}_{env}^{(b)^{(k)}}\right)}{ascii}.
#'       }
#'   }
#'
#' When the second level of bootstrapping is completed for all \mjeqn{b = 1,...,B}{ascii} then this
#' function reports the standard deviation of the bootstrapped envelope estimator of
#' expected Darwinian fitness. In this case, the bootstrap procedure accounts for model
#' selection volatility. The bootstrapped envelope estimator is
#'
#' \mjdeqn{\hat{\mu}_g = \frac{1}{B} \sum_{b=1}^B g(\hat{\tau}_{env}^{(b)})}{ascii}
#'
#' where \mjeqn{g(\hat{\tau}_{env}^{(b)})}{ascii} are the stored envelope estimators of expected Darwinian
#' fitness in the \code{env.boot.out} matrix included in the output of \code{fit.boot.Efron}.
#' The standard deviation of the bootstrapped envelope estimator of expected Darwinian fitness is
#'
#' \mjdeqn{ \sum_{b=1}^B\left[\widehat{cov}^{(b)^T}\hat{V}^{-1}\widehat{cov}^{(b)}\right] / B }{ascii}
#'
#' where \mjeqn{\widehat{cov}^{(b)} = \textbf{B}^{(b)^T} C^{(b)} / K}{ascii} and \mjeqn{\hat{V} = \textbf{B}^{(b)^T}\textbf{B}^{(b)}/K}{ascii}. The matrix \mjeqn{\textbf{B}^{(b)} \in R^{K\times p}}{ascii} has rows given by
#'
#' \mjdeqn{\hat{\tau}_{env}^{(b)^{(k)}} - \sum_{k=1}^K\hat{\tau}_{env}^{(b)^{(k)}}/K}{ascii}
#'
#' and the matrix \eqn{C^{(b)} \in R^{K \times d}} has columns given by
#'
#' \mjdeqn{g\left(\tau_{env}^{(b)^{(k)}}\right) - g\left(\tau_{env}^{(b)}\right)}{ascii}.
#'
#' For more details, see Efron (2014) and Eck (2015). The parametric bootstrap
#' procedure which uses the 1d algorithm to construct envelope estimators is
#' analogous to the above algorithm. To use the 1d algorithm, the user
#' specifies \code{method = "1d"} instead of \code{method = "eigen"}.
#' @param k The index of the top level parametric bootstrap procedure conducted by fit.boot.Efron that the second level of bootstrapping is being applied to.
#' @param nboot2 The bootstrap sample size for the second level of parametric bootstrapping.
#' @param out The output of fit.boot.Efron.
#' @param model An aster model object.
#' @param index The indices denoting which components of the canonical parameter vector are parameters of interest.
#' @param data An asterdata object corresponding to the original data.
#' @param amat This object can either be an array or a matrix. It specifies a linear combination of mean-value parameters that correspond to expected Darwinian fitness. See the aster function help page in the original aster package for more details.
#' @param newdata A dataframe corresponding to hypothetical individuals in which expected Darwinian fitness is to be estimated.
#' @param method The procedure used to obtain envelope estimators.
#' @returns
#'   \item{sd.Efron}{The estimated standard deviation (sd) for estimated expected Darwinian fitness where is estimation is conducted using envelope methodology. This sd accounts for model selection volatility. An eigenvalue decomposition using eigen is used internally to calculate this quantity.}
#'   \item{cov}{A components needed to construct sd.Efron if other numerical methods are desired.}
#'   \item{V}{A components needed to construct sd.Efron if other numerical methods are desired.}
#'   \item{MLE.tau.boot.subsample}{A components needed to construct sd.Efron if other numerical methods are desired.}
#'   \item{est.env.subsample}{A components needed to construct sd.Efron if other numerical methods are desired.}
#' @references Cook, R.D. and Zhang, X. (2014). Foundations for Envelope Models and Methods. \emph{JASA}, In Press.\cr
#' \cr
#' Cook, R.D. and Zhang, X. (2015). Algorithms for Envelope Estimation. \emph{Journal of Computational and Graphical Statistics}, Published online. \doi{10.1080/10618600.2015.1029577}.\cr
#' \cr
#' Eck, D. J., Geyer, C. J., and Cook, R. D. (2016). Enveloping the aster model. \emph{in prep}. \cr
#' \cr
#' Eck, D.~J., Geyer, C.~J., and Cook, R.~D. (2016). Web-based Supplementary Materials for ``Enveloping the aster model.'' \emph{in prep}. \cr
#' \cr
#' Efron, B. (2014). Estimation and Accuracy After Model Selection. \emph{JASA}, \strong{109:507}, 991-1007.\cr
#' @examples ### Web-based Supplementary Materials for ``Enveloping the aster model.'' ###
#' @export
#'

secondboot <- function(k, nboot2, out, model, index, data, amat,
                       newdata, method = c("eigen","1d"))
{

  # extract necessary components from the top-level of
  # bootstrapping
  env.boot.out <- out$env.boot.out
  env.1d.boot.out <- out$env.1d.boot.out
  MLE.boot.out <- out$MLE.boot.out
  MLE.tau.boot <- out$MLE.tau.boot
  env.tau.boot <- out$env.tau.boot
  env.1d.tau.boot <- out$env.1d.tau.boot
  P.list <- out$P.list
  P.1d.list <- out$P.1d.list
  vectors.list <- out$vectors.list
  u.1d.list <- out$u.1d.list
  nboot <- ncol(env.boot.out)
  npop <- nrow(env.boot.out)


  # specify necessary quantities for secondboot function
  # not extracted in the above
  aout4star2 <- b2 <- model
  beta <- model$coef
  p <- length(beta)
  n <- nrow(model$x)
  nnode <- ncol(model$x)
  modmat.mat <- matrix(model$modmat, nrow = n * nnode)
  mu <- predict(model, parm.type = "mean.value",
                model.type = "unconditional")
  #tau <- crossprod(modmat.mat, mu)
  tau <- MLE.tau.boot[, k]
  offset <- as.vector(model$origin)
  code <- data$code
  families <- data$families
  vars <- colnames(model$x)
  fam <- model$fam
  pred <- model$pred
  root <- model$root


  # necessary for trial run
  cov <- V <- var.Efron <- NULL
  MLE.tau.boot.subsample <- NULL
  est.env.subsample <- NULL


  # the array (matrix) that specifies Darwinian fitness
  amat.mat <- NULL
  if(class(amat) == "array"){
    amat.mat <- matrix(amat, nrow = npop,
                       byrow = TRUE)
  }
  if(class(amat) == "matrix"){
    amat.mat <- amat
    amat <- array(amat.mat, dim = c(npop, nnode, p))
  }

  # initialize important quantities
  P.foo <- M.foo <- fit.foo <- foo <- NULL
  vectors.foo <- u <- NULL

  if(method == "eigen"){
    foo <- env.tau.boot[, k]
    vectors.foo <- vectors.list[[k]]
    u <- length(vectors.foo)
    P.foo <- P.list[[k]]
    M.foo <- modmat.mat; M.foo[, index] <- M.foo[, index] %*% P.foo
    fit.foo <- rowSums(out$env.boot.out) / nboot
  }

  if(method == "1d"){
    u <- u.1d.list[[k]]
    fit.foo <- rowSums(out$env.1d.boot.out) / nboot
    if(u == length(index)){
      foo <- MLE.tau.boot[, k]
      P.foo <- P.1d.list[[k]]
      M.foo <- modmat.mat
    }
    if(u < length(index)){
      foo <- env.1d.tau.boot[, k]
      P.foo <- P.1d.list[[k]]
      M.foo <- modmat.mat; M.foo[, index] <- M.foo[, index] %*% P.foo
    }
  }

  # theta values
  theta.samp <- transformUnconditional(foo, M.foo, data,
                                       from = "tau", to = "theta", offset = offset)
  theta.samp <- matrix(theta.samp, nrow = n, ncol = nnode)

  # setup of model matrices and asterdata object
  # obtained from top level of bootstrapping
  cond <- !(colnames(modmat.renew) %in% model$dropped)
  modmat.renew <- modmat.env.renew <- modmat.renew[, cond]
  modmat.env.renew[, index] <-
    modmat.env.renew[, index] %*% P.foo
  data.renew <- asterdata(newdata, vars = vars, pred = pred,
                          group = rep(0, length(vars)), code = code,
                          families = families)
  origin.renew <- model$origin[1:npop,]
  offset.renew <- as.vector(origin.renew)

  # initial quantities
  #M.foo.array <- array(M.foo, dim = dim(model$modmat))
  #print(dim(M.foo.array))
  #print(dim(model$modmat))
  r <- length(index)
  U.2nd <- P.2nd <- Sigma.uu <- matrix(0, r, r)
  Gamma.2nd <- matrix(0, nrow = r, ncol = u)
  M.2nd <- M.foo
  est.env.subsample <- matrix(0, nrow = npop, ncol = nboot2)
  MLE.tau.boot.subsample <- matrix(0, nrow = p, ncol = nboot2)
  sd.Efron <- 0

  # generate many new resamples
  for(j in 1:nboot2){

    xstar.samp <- raster(theta.samp, pred, fam, root)

    # fit the aster model using the secondary generated data
    class(b2) <- class(try(aout4star2 <- aster(xstar.samp, root,
                                               pred, fam, model$modmat, parm = beta), silent = TRUE))[1]

    if(class(b2) == "try-error"){
      class(b2) <- class(try(aout4star2 <- aster(xstar.samp, root, pred,
                                                 fam, model$modmat, parm = beta,
                                                 method = "nlm"), silent = TRUE))[1]
    }

    if(class(b2) == "try-error"){
      class(b2) <- class(try( aout4star2 <- aster(xstar.samp, root,
                                                  pred, fam, model$modmat), silent = TRUE))[1]
    }

    if(class(b2) == "try-error"){
      class(b2) <- class(try(aout4star2 <- aster(xstar.samp, root, pred,
                                                 fam, model$modmat, method = "nlm"), silent = TRUE))[1]
    }


    # get tau and mu from this fit
    mu.star <- predict(aout4star2, parm.type = "mean.value",
                       model.type = "unconditional")

    #tau.renew <- transformUnconditional(aout4star2$coef, data,
    #  modmat = modmat.mat, from = "beta", to = "tau",
    #  offset = offset)

    # get beta and tau using the model matrix containing
    # the projection into the envelope
    Sigma.uu <- aout4star2$fisher[index, index]
    if(method == "eigen"){
      if(length(vectors.foo) < length(index)){
        Gamma.2nd <- eigen(Sigma.uu, symmetric = TRUE)$vec[, vectors.foo]
        P.2nd <- projection(Gamma.2nd)
        M.2nd[, index] <- M.2nd[, index] %*% P.2nd
      }
      if(length(vectors.foo) == length(index)){
        M.2nd <- M.foo
      }
    }
    if(method == "1d"){
      if(u < length(index)){
        tau.star <- crossprod(M.foo, mu.star)
        U.2nd <- tau.star[index] %*% t(tau.star[index])
        Gamma.2nd <- manifold1Dplus(Sigma.uu, U = U.2nd, u = u)
        P.2nd <- tcrossprod(Gamma.2nd)
      }
      if(u == length(index)){
        P.2nd <- diag(length(index))
        M.2nd[, index] <- M.foo[, index] %*% P.2nd
      }
    }


    tau.env.star <- crossprod(M.2nd, mu.star)
    beta.env.star <- transformUnconditional(parm = tau.env.star,
                                            M.2nd, data, from = "tau", to = "beta",
                                            offset = offset)


    # compute the envelope estimator of expected
    # Darwinian fitness
    modmat.env.renew[, index] <- modmat.renew[,index] %*% P.2nd
    mu.env.renew <- transformUnconditional(parm = beta.env.star,
                                           modmat.env.renew, data.renew, from = "beta", to = "mu",
                                           offset = offset.renew)
    est.env.subsample[, j] <- (amat.mat %*% mu.env.renew)
    M.2nd <- M.foo

    # store the canonical statistic value
    MLE.tau.boot.subsample[, j] <- tau.env.star

    # compute the Efron sd estimator99i5
    if(j == nboot2){

      B <- t(MLE.tau.boot.subsample - rowSums(MLE.tau.boot.subsample)/nboot2)
      V <- (t(B) %*% B) / nboot2;  eig.V <- eigen(V)
      env.centered <- t(est.env.subsample - fit.foo)
      cov <- t(B) %*% env.centered / nboot2

      var.Efron <- t(cov) %*% eig.V$vec %*% diag(1/eig.V$val) %*%
        t(eig.V$vec) %*% cov
      sd.Efron <- sqrt(diag(var.Efron))
    }
  }

  ### output
  out <- list(sd.Efron = sd.Efron, cov = cov, V = V,
              MLE.tau.boot.subsample = MLE.tau.boot.subsample,
              est.env.subsample = est.env.subsample)
  return(out)
}
