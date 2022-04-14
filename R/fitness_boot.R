#' fitness_boot
#' @description
#' \loadmathjax
#' This function implements the first level of the parametric bootstrap procedure given by Algorithm 1 in Eck et al. (2020) with 
#' respect to the submodel mean-value parameterization (parameterization closest to that of expected Darwinian fitness). 
#' This is detailed in Steps 1 through 3d in the algorithm below. This parametric bootstrap generates resamples from the 
#' distribution evaluated at an envelope estimator of \mjeqn{\tau}{ascii} adjusting for model selection volatility.\cr
#' @details
#' This function implements the first level of the parametric bootstrap procedure given 
#' by Algorithm 1 in Eck (2020) with respect to either the 1d (1d Algorithm in Cook and Zhang 
#' (2015 a,b)) or eigen (see Section 4 in Eck et al. (2020)) approaches. 
#' This is detailed in Steps 1 through 3d in the algorithm below.  This parametric bootstrap generates 
#' resamples from the distribution evaluated at an envelope estimator of \mjeqn{\tau}{ascii}
#' adjusting for model selection volatility.
#' The user specifies a model selection criterion which selects vectors that
#' construct envelope estimators using the reducing subspace approach. The user also
#' specifies which method is to be used in order to calculate envelope
#' estimators. When one is using a partial envelope, then this function
#' constructs envelope estimators of \mjeqn{\upsilon}{ascii} where we write \mjeqn{\tau}{ascii} = \mjeqn{(\gamma^T,\upsilon^T)^T}{ascii}
#' and \mjeqn{\upsilon}{ascii} corresponds to aster model parameters of interest.
#'
#' In applications, candidate reducing subspaces are indices of eigenvectors of \mjeqn{ \widehat{\Sigma}_{\upsilon,\upsilon} }{ascii}
#'  where \mjeqn{ \widehat{\Sigma}_{\upsilon,\upsilon}}{ascii} is the part of \mjeqn{\hat{\Sigma}}{ascii}
#' corresponding to our parameters of interest. These indices are specified
#' by \code{vectors}. When all of the components of \mjeqn{\tau}{ascii} are components
#' of interest, then we write \mjeqn{\widehat{\Sigma}_{\upsilon,\upsilon} = \widehat{\Sigma} }{ascii}. When data
#' is generated via the parametric bootstrap, it is the indices (not the
#' original reducing subspaces) that are used to construct envelope estimators
#' constructed using the generated data. The algorithm using reducing subspaces
#' is as follows:
#'
#'   \enumerate{
#'     \item Fit the aster model to the data and obtain
#'     \mjeqn{ \hat{\tau} = (\hat{\gamma}^T, \hat{\upsilon}^T) }{ascii} and \mjeqn{\hat{\Sigma}}{ascii}
#'    from the aster model fit.
#'     \item Compute the envelope estimator of \mjeqn{\upsilon}{ascii} in the original sample, given as
#'     \mjeqn{ \hat{\upsilon_{env}} = P_{\hat{G}}\hat{\upsilon}}{ascii} where \mjeqn{P_{\hat{G}}}{ascii} is computed using reducing subspaces
#'     and selected via a model selection criterion of choice.
#'     \item Perform a parametric bootstrap by generating resamples from the distribution of the
#'     aster submodel evaluated at \mjeqn{\hat{\tau}_{env} = (\hat{\gamma}^T,\hat{\upsilon_{env}}^T)^T}{ascii}. For iteration
#'     \mjeqn{b=1,...,B}{ascii} of the procedure:
#'       \enumerate{
#'         \item Compute \mjeqn{\hat{\tau}^{(b)}}{ascii} and \mjeqn{\hat{\Sigma}_{\upsilon,\upsilon}^{(b)}}{ascii} from the aster
#'         model fit to the resampled data.
#'         \item Build \mjeqn{P_{\hat{G}}^{(b)}}{ascii} using the indices of \mjeqn{\widehat{\Sigma}_{\upsilon,\upsilon}^{(b)}}{ascii} that
#'         are selected using the same model selection criterion as Step 2 to build \mjeqn{\hat{G}}{ascii}.
#'         \item Compute \mjeqn{\hat{\upsilon}_{env}^{(b)} = P_{\hat{\mathcal{E}}}^{(b)}\hat{\upsilon}^{(b)}}{ascii} and
#'         \mjeqn{\hat{\tau}_{env}^{(b)} = \left(\hat{\gamma}^{(b)^T},\hat{\upsilon}_{env}^{(b)^T}\right)^T}{ascii}.
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
#'#' For more details, see Eck et al. (2020) and Efron (2014). The parametric bootstrap
#' procedure which uses the 1d algorithm to construct envelope estimators is
#' analogous to the above algorithm. To use the 1d algorithm, the user
#' specifies \code{method = "1d"} instead of \code{method = "eigen"}.
#' @param model An aster model object.
#' @param nboot The number of bootstrap iterations desired.
#' @param index The indices denoting which components of the canonical parameter vector are parameters of interest.
#' @param vectors The indices denoting which reducing subspace of Fisher information is desired to construct envelope estimators. \cr
#' Must be specified if \code{method = "eigen"}.
#' @param dim The dimension of the envelope space used to construct envelope estimators.\cr
#'  Must be specified if \code{method = "1d"}.
#' @param data An asterdata object corresponding to the original data.
#' @param amat This object can either be an array or a matrix.\cr
#'  It specifies a linear combination of mean-value parameters that correspond to expected Darwinian fitness.\cr
#'  See the \code{aster} function help page in the original \code{aster} package for more details.
#' @param newdata A dataframe corresponding to hypothetical individuals in which expected Darwinian fitness is to be estimated.
#' @param modmat.new A model matrix corresponding to hypothetical individuals in which expected Darwinian fitness is to be estimated.
#' @param renewdata A dataframe in long format corresponding to hypothetical individuals in which expected Darwinian fitness is to be estimated.
#' @param criterion A model selection criterion of choice.
#' @param alpha The type 1 error rate desired for the LRT.
#' @param fit.name An expression that appears in the name of the nodes that correspond to Darwinian fitness.\cr
#'  This is only necessary if \code{renewdata} is not provided.
#' @param method The procedure used to obtain envelope estimators.
#' @param quiet A logical argument. If FALSE, the function displays how much time it takes to run \code{m} iterations.
#' @param corenum The number of cores specified for speeding up the bootstrap process.
#' @return a list containing the following elements in order:
#' \item{env.boot.out}{Estimated expected Darwinian fitness from generated data obtained from Steps 3a-3d in the bootstrap procedure using the envelope
#' estimator constructed using reducing subspaces. }
#' \item{MLE.boot.out}{Estimated expected Darwinian fitness from generated data obtained from Steps 3a-3d in the bootstrap procedure using the MLE. }
#' \item{MLE.tau.boot}{Estimated mean-value parameter vectors from generated data obtained from Steps 3a-3d in the bootstrap procedure using the MLE. }
#' \item{env.tau.boot}{Estimated mean-value parameter vectors from generated data obtained from Steps 3a-3d in the bootstrap procedure using the envelope estimator constructed using the 1d algorithm. }
#' \item{P.list}{A list of all estimated projections into the envelope space constructed from reducing subspaces for Steps 3a-3d in the bootstrap
#' procedure. }
#' \item{vectors.list}{A list of indices of eigenvectors used to build the projections in P.list. These indices are selected using the user specified model selection criterion as indicated in Steps 3a-3d in the bootstrap procedure. }
#' @references
#'
#' Eck, D. J., Geyer, C. J., and Cook, R. D. (2020). Combining envelope methodology and aster models for variance reduction in life 
#' history analyses. \emph{Journal of Statistical Planning and Inference}, \strong{205}, 283-292. \cr 
#' \cr
#' Eck, D.~J., Geyer, C.~J., and Cook, R.~D. (2018). Supporting Data Analysis for 
#' "Combining Envelope Methodology and Aster Models for Variance Reduction in Life History Analyses." \cr
#' \cr
#' Cook, R.D. and Zhang, X. (2015 a). Foundations for Envelope Models and Methods. 
#' \emph{Journal of the American Statistical Association}, \strong{110}, 599-611. \cr
#' \cr
#' Cook, R.D. and Zhang, X. (2015 b). Algorithms for Envelope Estimation. 
#' \emph{Journal of Computational and Graphical Statistics}, Published online. \doi{10.1080/10618600.2015.1029577}.\cr
#' \cr
#' Efron, B. (2014). Estimation and Accuracy After Model Selection. \emph{Journal of the American Statistical Association}, 
#' \strong{109:507}, 991-1007.
#' @noMd
#' @export
#' @examples \dontrun{# see Supporting Data Analysis for 
#' "Combining Envelope Methodology and Aster Models for Variance Reduction in Life History Analyses."}
#' @import aster
#' @import doSNOW
#' @import foreach
#' @import parallel
#' @import stats
#' @importFrom aster2 transformUnconditional transformSaturated asterdata
#' @importFrom  mathjaxr preview_rd

fitness_boot <- function(model, nboot, index, vectors = NULL, dim = NULL,
    data, amat, newdata, modmat.new = NULL, renewdata = NULL,
    criterion = c("AIC","BIC","LRT"), alpha = 0.05, fit.name = NULL,
    method = c("eigen","1d"), quiet = FALSE,corenum=NULL)
{

  # stopping condition
  #stopifnot(!(is.null(vectors) & is.null(u)))
  #stopifnot(!(is.null(fit.name) & is.null(renewdata)))
  #timer <- proc.time()

  # extract important envelope initial quantities
  # obtain the target and nuisance parameters
  # build U w.r.t. the target parameters
  beta <- model$coef
  x <- model$x
  vars <- colnames(x)
  code <- data$code
  families <- data$families
  fam <- model$fam
  pred <- model$pred
  root <- model$root
  n <- nrow(x)
  nnode <- ncol(x)
  npop <- nrow(newdata)
  p <- length(beta)
  formula <- model$formula
  modmat.model <- model$modmat
  modmat.mat <- matrix(modmat.model, nrow = n * nnode)
  dimensions <- dim(modmat.model)
  offset <- as.vector(model$origin)

  # obtain tau
  mu <- predict(model, parm.type = "mean.value",
    model.type = "unconditional")
  tau <- crossprod(modmat.mat, mu)
  nuis.ind <- c(1:p)[!index]
  k <- length(index)
  target <- tau[index]

  # obtain the projection used to construct the envelope
  # estimator
  avar <- (model$fisher)[index,index]
  U <- P <- NULL; u <- length(index)
  if(method == "eigen"){
    u <- length(vectors)
    eig <- eigen(avar, symmetric = TRUE)
    G <- eig$vec[,c(vectors)]
    P <- tcrossprod(G)
  }
  if(method == "1d"){
    u <- dim
    U <- target %o% target
    G <- manifold1Dplus(M = avar, U = U, u = u)
    P <- tcrossprod(G)
  }

  tau.env <- crossprod(P,target)
  fulltau <- tau
  fulltau[index] <- tau.env

  # change the model matrix for the envelope estimator
  modelmatrix.int <- modmat.mat
  modelmatrix.int[, index] <- modelmatrix.int[, index] %*% P
  modmat.model.int <- array(modmat.mat, dimensions)

  # convert from tau to beta changed from fulltau to tau
  beta.foo <- transformUnconditional(parm = fulltau,
    modelmatrix.int, data, from = "tau", to = "beta",
    offset = offset)

  # set up for the bootstrap for the MLE
  theta.hat2 <- predict(model, model.type = "cond",
    parm.type = "canon")
  theta.hat2 <- matrix(theta.hat2, nrow = n, ncol = nnode)
  MLE.tau.boot <- matrix(nrow = npop, ncol = nboot) # changed from nrow = p
  b2 <- "try-error"; class(b2) <- "try-error"

  # set up for the bootstrap for the envelope estimator
  # (this envelope estimator corresponds to the selected method)
  theta.hat <- transformUnconditional(parm = fulltau,
    modelmatrix.int, data, from = "tau", to = "theta",
    offset = offset)
  theta.hat <- matrix(theta.hat, nrow = n, ncol = nnode)
  #theta.hat <- theta.hat2
  b <- "try-error"; class(b) <- "try-error"


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


  # initial quantities
  est <- est2 <- matrix(0, ncol = npop)
  env.tau.boot <- matrix(0, ncol = p)
  P.list <-  NULL
  class(P.list) <- class(P.list) <- "list"
  vectors.list <- u.1d.list <- NULL
  class(vectors.list) <- class(u.1d.list) <- "list"
  out <- NULL
  table <- NULL


  #modmat.renew <- model.matrix(m1$formula, data = renewdata)
  cond <- !(colnames(modmat.renew) %in% model$dropped)
  modmat.renew <- modmat.renew[, cond]
  M1.renew <- t(modmat.renew[,-index])
  M2.renew <- P %*% t(modmat.renew[,index])
  modmat.env.renew <- t(rbind(M1.renew,M2.renew))
  data.renew <- asterdata(newdata, vars = vars, pred = pred,
    group = rep(0, length(vars)), code = code,
    families = families)
  origin.renew <- model$origin[1:npop,]
  offset.renew <- as.vector(origin.renew)

  # fit the aster model using different methods

  check_error = function(class, modmat_status=modmat.model) {
    if(class == "try-error"){
      class <- class(try(aout4star2 <- aster(xstar2, root, pred,
        fam, modmat_status, parm = beta,
        method = "nlm"), silent = TRUE))[1]
    }

    if(class == "try-error"){
      class <- class(try( aout4star2 <- aster(xstar2, root,
        pred, fam, modmat_status), silent = TRUE))[1]
    }

    if(class == "try-error"){
      class<- class(try(aout4star2 <- aster(xstar2, root, pred,
        fam, modmat_status, method = "nlm"), silent = TRUE))[1]
    }
    return(class)
  }
  method_check = function(method, tau, index, model, data, alpha) {
    result = 0
    if (method == "1d") {
      barbaz <- envlpaster::selection(tau, index, model, data,
        alpha = alpha, type = "mean-value", method = "1d")
      if(criterion == "AIC") result <- barbaz$aic
      if(criterion == "BIC") result <- barbaz$bic
      if(criterion == "LRT") result <- barbaz$LRT
    }
    if(method == "eigen"){
      fubar <- envlpaster::selection(tau, index, model, data,
        alpha = alpha, type = "mean-value", method = "eigen")
      if(criterion == "AIC") result <- fubar$aic
      if(criterion == "BIC") result <- fubar$bic
      if(criterion == "LRT") result <- fubar$LRT
    }
    return(result)
  }


  #list of vars and functs that need to be move into the parallel from global
  cl <- makeCluster(corenum*2,type="SOCK")
  registerDoSNOW(cl)
  clusterExport(cl, c('check_error','method_check','theta.hat2',
    'pred', 'fam', 'root','modmat.model','beta',
    'modmat.mat','offset','amat.mat','vectors.list',
    'MLE.tau.boot','est2','b2','offset.renew','index'),
    envir=environment())
  clusterEvalQ(cl,library("aster2","envlpaster","aster"))
  resultlist<-foreach(k=1:nboot) %dopar%{

    # fit the aster model to the resampled data obtained from
    # the MLE of theta
    xstar2 <-aster::raster(theta.hat2, pred, fam, root)
  # fit the aster model to the resampled data obtained from
    # the MLE of theta
    class(b2) <- class(try(aout4star2 <-aster::aster(xstar2, root, pred,
      fam, modmat.model, parm = beta), silent = TRUE))[1]
    class(b2) <- check_error(class(b2))
    # MLE of expected Darwinian fitness
    tau.renew <- transformUnconditional(aout4star2$coef, data,
      modmat = modmat.mat, from = "beta", to = "tau",
      offset = offset)
    phi.MLE.renew <- offset.renew + modmat.renew %*% aout4star2$coef
    mu.MLE.renew <- transformSaturated(parm = phi.MLE.renew,
      data.renew, from = "phi", to = "mu")
    MLE.tau.boot <- tau.renew
    est2 <- amat.mat %*% mu.MLE.renew #Dar.fit.MLE
    # selection of the envelope eigenstructure
    vectors.list = method_check(method, tau = tau.renew, index=index, model=aout4star2, data=data, alpha=alpha)
    dim.1d=0
    if (method=="1d") {
      dim.1d=vectors.list
    }

    # get the bootstrapped envelope estimators
    M <- matrix(modmat.model.int, nrow = n * nnode)
    avar.star <- (aout4star2$fisher)[index,index]
    beta.env.star <- aout4star2$coef
    mu.star <- predict(aout4star2, parm.type = "mean.value",
      model.type = "unconditional")
    modelmatrix.star <- M

    tau.env.star <- crossprod(modelmatrix.star, mu.star)
    env.tau.boot <- tau.env.star

    # envelope estimator of expected Darwinian fitness using
    # eigenstructures
    if(method == "eigen"){
      G.star <- eigen(avar.star, symmetric = TRUE)$vec[,c(vectors)]
      P.star <- tcrossprod(G.star)
      P.list <- P.star
      modelmatrix.star[, index] <- modelmatrix.star[, index] %*% P.star
      M1.renew <- t(modmat.renew[,-index])
      M2.renew <- P.star %*% t(modmat.renew[,index])
      modmat.env.renew <- t(rbind(M1.renew,M2.renew))
    }

    # envelope estimator of expected Darwinian fitness using the
    # 1d algorithm
    if(method == "1d"){
      if(vectors.list < length(index)){
        up.env.star <- crossprod(M, mu.star)[index]
        U.star <- up.env.star %o% up.env.star
        G.star <- manifold1Dplus(avar.star, U = U.star, u = dim.1d)
        P.star <- tcrossprod(G.star)
        P.list[k] <- P.star

        modelmatrix.star[, index] <- modelmatrix.star[, index] %*% P.star
        beta.env.star <- transformUnconditional(parm = tau.env.star,
          modelmatrix.star, data, from = "tau", to = "beta",
          offset = offset)
        modmat.env.renew <- modmat.renew

      }
      if(dim.1d == length(index)){
        P.list <- diag(length(index))
        env.tau.boot<- tau.renew
      }
      if(dim.1d == length(index)){
        env.tau.boot<- tau.renew
      }
    }
    # P.list
    # env.tau.boot
    beta.env.star <- transformUnconditional(parm = tau.env.star,
        modelmatrix.star, data, from = "tau", to = "beta",
        offset = offset)
    modmat.env.renew[, index] <- modmat.env.renew[, index] %*% P.star
        mu.env.renew <- transformUnconditional(parm = beta.env.star,
    modmat.env.renew, data.renew, from = "beta", to = "mu",
          offset = offset.renew)

    # the stored expected Darwinian fitness estimates
    est<- amat.mat %*% mu.env.renew
    list(est,
         est2,
         MLE.tau.boot,
        env.tau.boot,
        P.list,
        vectors.list)
  }
  stopCluster(cl)
  env.boot.out=do.call(cbind, lapply(1:nboot, function(j){
    unlist(resultlist[[j]][1])
  }))
   MLE.boot.out=do.call(cbind, lapply(1:nboot, function(j){
    unlist(resultlist[[j]][2])
  }))
  MLE.tau.boot=do.call(cbind, lapply(1:nboot, function(j){
    unlist(resultlist[[j]][3])
  }))
  env.tau.boot=do.call(cbind, lapply(1:nboot, function(j){
    unlist(resultlist[[j]][4])
  }))
   P.list=do.call('c', lapply(1:nboot, function(j){
    resultlist[[j]][5]
  }))
  vectors.list=do.call('c', lapply(1:nboot, function(j){
    resultlist[[j]][6]
  }))
  out<-list( env.boot.out = env.boot.out, MLE.boot.out = MLE.boot.out, MLE.tau.boot = MLE.tau.boot,
        env.tau.boot = env.tau.boot,
        P.list = P.list,
        vectors.list = vectors.list)

  #secondboot

  return(out)
}





