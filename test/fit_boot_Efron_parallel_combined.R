
library(foreach)
library(doSNOW)

fit.boot.Efron.par <- function(model, nboot,nboot2, index, vectors = NULL, dim = NULL,
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
  # necessary for trial run
  cov <- V <- var.Efron <- NULL
  MLE.tau.boot.subsample <- NULL
  est.env.subsample <- NULL

  # initialize important quantities
  P.foo <- M.foo <- fit.foo <- foo <- NULL
  vectors.foo <- u <- NULL

  if(method == "eigen"){
    foo <- env.tau.boot[, k] 
    vectors.foo <- vectors.list[[k]]
    u <- length(vectors.foo)
    P.foo <- P.list[[k]] 
    M.foo <- modmat.mat; M.foo[, index] <- M.foo[, index] %*% P.foo
    fit.foo <- rowSums(env.boot.out) / nboot
  }

  if(method == "1d"){
    u <- u.1d.list[[k]]
    fit.foo <- rowSums(env.boot.out) / nboot
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
  secondstrap=function(x) {
    for(j in 1:nboot2){
        xstar.samp <- raster(theta.samp, pred, fam, root)
        
        # fit the aster model using the secondary generated data
        class(b2) <- class(try(aout4star2 <- aster(xstar.samp, root, 
          pred, fam,  modmat.model, parm = beta), silent = TRUE))[1] 
        class(b2) <- check_error(class(b2))
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
    result <- list(sd.Efron = sd.Efron, cov = cov, V = V, 
        MLE.tau.boot.subsample = MLE.tau.boot.subsample, 
        est.env.subsample = est.env.subsample)
    return(result)
  }
  out2=lapply(1:nboot, secondstrap)
  ### output 
  
  return(list(out2))
}





