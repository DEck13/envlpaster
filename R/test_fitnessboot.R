#install.packages("aster2")
#install.packages("envlpaster")
library(aster2)
library(envlpaster)
data(simdata30nodes)
data <- simdata30nodes.asterdata
set.seed(13)
vars <- variables
pred <- data$pred
families <- data$families
code <- data$code
fam <- code
fam[fam == 2] <- 3
nnode <- length(vars)
xnew <- as.matrix(simdata30nodes[,c(1:nnode)])
m1 <- aster(xnew, root, pred, fam, modmat)
m1$formula <- formula
m1$xlevels <- xlevels
m1$terms <- terms
summary(m1)

beta <- m1$coef; p <- length(beta)
target <- 5:9; m <- length(target)
x <- m1$x
nind <- n <- nrow(x)

mu <- predict(m1, parm.type = "mean.value",
              model.type = "unconditional")
modmat.mat <- matrix(modmat, ncol = p)
tau <- crossprod(modmat.mat, mu)

vars <- as.vector(outer(c("u", "v", "w"), 1:10, paste,
                        sep = ""))

nx <- 12; npop <- nx^2
z1.cand <- seq(from = -6, to = 6, length = nx)
z2.cand <- seq(from = -6, to = 6, length = nx)
fred <- data.frame( cbind(rep(z1.cand, each = nx),
                          rep(z2.cand, nx), matrix(1,nrow = nx^2, ncol = 30)))
colnames(fred) <- c("z1","z2",vars)
fred$root <- 1

amat <- array(0, c(npop, nnode, npop))
foot <- grepl("w", vars)
for (k in 1:npop) amat[k, foot, k] <- 1

renewdata <- reshape(fred, varying = list(vars),
                     direction = "long", timevar = "varb",
                     times = as.factor(vars), v.names = "resp")
layer2 <- gsub("[0-9]", "",
               as.character(renewdata$varb))
renewdata <- data.frame(renewdata, layer = layer2)
renewdata$vtype <- as.factor(substr(as.character(
  renewdata$varb),1, 1))
renewdata$year <- as.numeric(substring(as.character(
  renewdata$varb), 2))
renewdata$uyear <- renewdata$year *
  as.numeric(as.character(renewdata$vtype) == "u")

modmat.renew <- model.matrix(m1$formula, data = renewdata)
modmat.renew.array <- array(modmat.renew, dim = c(npop, nnode, p))

p1 <- predict(m1, modmat = modmat.renew.array, varvar = varb,
              idvar = id, root = renewdata$root, se.fit = TRUE)
p1.mat <- matrix(p1$fit, ncol = nnode)
p1.mat.fit <- rowSums(
  p1.mat[, which(grepl("w", unique(renewdata$varb)))]
)

# for foreach
library(parallel)
library(foreach)
library(doParallel)

fit.boot.Efron <- function(model, nboot, index, vectors = NULL, dim = NULL,
                           data, amat, newdata, modmat.new = NULL, renewdata = NULL, 
                           criterion = c("AIC","BIC","LRT"), alpha = 0.05, fit.name = NULL, 
                           method = c("eigen","1d"), quiet = FALSE)
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
  est  <- matrix(nrow = npop, ncol = nboot) # changed from nrow = p
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
  est <- est2 <- est.1d <- matrix(0, nrow = npop, ncol = nboot)
  MLE.tau.boot <- env.tau.boot <- env.1d.tau.boot <- matrix(0, nrow = p, ncol = nboot)
  P.list <- P.1d.list <-  NULL
  class(P.list) <- class(P.1d.list) <- "list"
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
  
  c1 <- makeCluster(numCores)
  registerDoParallel(c1)
  
  
  
  # the parametric bootstrap which takes Efron's procedure into account
  r1 <- foreach(iboot= 1:nboot, .combine = "cbind", .packages = c("envlpaster", "aster2")) %dopar% {
  #for(k in 1:nboot){
    
    # generate a resample of responses from the MLE
    xstar2 <- raster(theta.hat2, pred, fam, root)
    
    # fit the aster model to the resampled data obtained from 
    # the MLE of theta
    class(b2) <- class(try(aout4star2 <- aster(xstar2, root, pred,
                                               fam, modmat.model, parm = beta), silent = TRUE))[1] 
    
    if(class(b2) == "try-error"){
      class(b2) <- class(try(aout4star2 <- aster(xstar2, root, pred,
                                                 fam, modmat.model, parm = beta, 
                                                 method = "nlm"), silent = TRUE))[1] 
    }
    
    if(class(b2) == "try-error"){
      class(b2) <- class(try( aout4star2 <- aster(xstar2, root, 
                                                  pred, fam, modmat.model), silent = TRUE))[1]
    }
    
    if(class(b2) == "try-error"){
      class(b2) <- class(try(aout4star2 <- aster(xstar2, root, pred,
                                                 fam, modmat.model, method = "nlm"), silent = TRUE))[1] 
    }
    
    
    # MLE of expected Darwinian fitness
    tau.renew <- transformUnconditional(aout4star2$coef, data,
                                        modmat = modmat.mat, from = "beta", to = "tau",
                                        offset = offset)
    phi.MLE.renew <- offset.renew + modmat.renew %*% aout4star2$coef
    mu.MLE.renew <- transformSaturated(parm = phi.MLE.renew, 
                                       data.renew, from = "phi", to = "mu")
    MLE.tau.boot[, k] <- tau.renew
    Dar.fit.MLE <- amat.mat %*% mu.MLE.renew
    
    
    # generate a resample of responses from the 
    # envelope estimator
    xstar <- raster(theta.hat, pred, fam, root)
    
    # fit the aster model to the resampled data obtained from 
    # the envelope estimator of theta
    class(b2) <- class(try(aout4star2 <- aster(xstar, root, pred,
                                               fam, modmat.model.int, parm = beta), silent = TRUE))[1] 
    
    if(class(b2) == "try-error"){
      class(b2) <- class(try(aout4star2 <- aster(xstar, root, pred,
                                                 fam, modmat.model.int, parm = beta, 
                                                 method = "nlm"), silent = TRUE))[1] 
    }
    
    if(class(b2) == "try-error"){
      class(b2) <- class(try( aout4star2 <- aster(xstar, root, 
                                                  pred, fam, modmat.model.int), silent = TRUE))[1]
    }
    
    if(class(b2) == "try-error"){
      class(b2) <- class(try(aout4star2 <- aster(xstar, root, pred,
                                                 fam, modmat.model.int, method = "nlm"), silent = TRUE))[1] 
    }
    
    
    # MLE of expected Darwinian fitness
    tau.renew <- transformUnconditional(aout4star2$coef, data,
                                        modmat = modelmatrix.int, from = "beta", to = "tau",
                                        offset = offset)
    #phi.MLE.renew <- offset.renew + modmat.renew %*% aout4star2$coef
    #mu.MLE.renew <- transformSaturated(parm = phi.MLE.renew, 
    #  data.renew, from = "phi", to = "mu")
    
    # selection of the envelope dimension using the 1d algorithm
    barbaz <- NULL
    dim.1d <- 0
    if(method == "1d"){
      barbaz <- selection(tau.renew, index, model = aout4star2, data = data, 
                          alpha = alpha, type = "mean-value", method = "1d")
      if(criterion == "AIC") dim.1d <- barbaz$aic
      if(criterion == "BIC") dim.1d <- barbaz$bic
      if(criterion == "LRT") dim.1d <- barbaz$LRT
    }
    
    u.1d.list[[k]] <- dim.1d
    
    # selection of the envelope eigenstructure
    fubar <- NULL
    vectors <- NULL
    if(method == "eigen"){
      fubar <- selection(tau.renew, index, model = aout4star2, data = data, 
                         alpha = alpha, type = "mean-value", method = "eigen")
      if(criterion == "AIC") vectors <- fubar$aic
      if(criterion == "BIC") vectors <- fubar$bic
      if(criterion == "LRT") vectors <- fubar$LRT    
    }
    
    vectors.list[[k]] <- vectors
    
    if(!quiet){
      cat("iteration: " , k, " indices: ", vectors, " dim.1d: ", dim.1d, "\n")
    }
    
    # get the bootstrapped envelope estimators    
    M <- matrix(modmat.model.int, nrow = n * nnode)
    avar.star <- (aout4star2$fisher)[index,index]
    beta.env.star <- aout4star2$coef
    mu.star <- predict(aout4star2, parm.type = "mean.value", 
                       model.type = "unconditional")
    
    # envelope estimator of expected Darwinian fitness using
    # eigenstructures
    Dar.fit.env <- Dar.fit.MLE
    if(method == "eigen"){    
      G.star <- eigen(avar.star, symmetric = TRUE)$vec[,c(vectors)]  
      P.star <- tcrossprod(G.star)
      P.list[[k]] <- P.star
      modelmatrix.star <- M
      modelmatrix.star[, index] <- modelmatrix.star[, index] %*% P.star
      
      tau.env.star <- crossprod(modelmatrix.star, mu.star)
      env.tau.boot[, k] <- tau.env.star
      beta.env.star <- transformUnconditional(parm = tau.env.star, 
                                              modelmatrix.star, data, from = "tau", to = "beta", 
                                              offset = offset)
      M1.renew <- t(modmat.renew[,-index])
      M2.renew <- P.star %*% t(modmat.renew[,index])
      modmat.env.renew <- t(rbind(M1.renew,M2.renew))      
      mu.env.renew <- transformUnconditional(parm = beta.env.star, 
                                             modmat.env.renew, data.renew, from = "beta", to = "mu", 
                                             offset = offset.renew)
      Dar.fit.env <- amat.mat %*% mu.env.renew  
    }
    
    # envelope estimator of expected Darwinian fitness using the
    # 1d algorithm
    Dar.fit.env.1d <- Dar.fit.MLE  
    if(method == "1d"){
      if(dim.1d < length(index)){
        up.env.star <- crossprod(M, mu.star)[index]
        U.1d.star <- up.env.star %o% up.env.star
        G.1d.star <- manifold1Dplus(avar.star, U = U.1d.star, u = dim.1d)
        P.1d.star <- tcrossprod(G.1d.star)
        P.1d.list[[k]] <- P.1d.star
        
        modelmatrix.star <- M
        modelmatrix.star[, index] <- modelmatrix.star[, index] %*% P.1d.star
        
        tau.env.star <- crossprod(modelmatrix.star, mu.star)
        env.1d.tau.boot[, k] <- tau.env.star
        beta.env.star <- transformUnconditional(parm = tau.env.star, 
                                                modelmatrix.star, data, from = "tau", to = "beta", 
                                                offset = offset)
        modmat.env.renew <- modmat.renew
        modmat.env.renew[, index] <- modmat.env.renew[, index] %*% P.1d.star
        mu.env.renew <- transformUnconditional(parm = beta.env.star, 
                                               modmat.env.renew, data.renew, from = "beta", to = "mu", 
                                               offset = offset.renew)
        Dar.fit.env.1d <- amat.mat %*% mu.env.renew 
      }
      
      if(dim.1d == length(index)){ 
        P.1d.list[[k]] <- diag(length(index))
        env.1d.tau.boot[, k] <- tau.renew
      }
    }
    
    
    # the stored expected Darwinian fitness estimates
    est[, k] <- Dar.fit.env
    est2[, k] <- Dar.fit.MLE
    est.1d[, k] <- Dar.fit.env.1d
    
    
    if(k == nboot){
      #means <- apply(est, FUN = mean, MARGIN = 1)
      #S <- var(t(est)); S2 <- var(t(est2)); S.1d <- var(t(est.1d))
      #ratio <- sqrt( diag(S2) / diag(S) )
      #table <- cbind(Dar.fit.env, sqrt(diag(S)), Dar.fit.MLE, 
      #  sqrt(diag(S2)), ratio)
      #colnames(table) <- c("env","se(env)","MLE","se(MLE)","ratio")  
      
      #ratio.1d <- sqrt( diag(S.1d) / diag(S) )
      #table.1d <- cbind(Dar.fit.env, sqrt(diag(S)), Dar.fit.env.1d, 
      #  sqrt(diag(S.1d)), ratio.1d)
      #colnames(table.1d) <- c("env","se(env)","env.1d","se(env.1d)","ratio")  
      out <- list( env.boot.out = est, MLE.boot.out = est2, 
                   env.1d.boot.out = est.1d, MLE.tau.boot = MLE.tau.boot, 
                   env.tau.boot = env.tau.boot, env.1d.tau.boot = env.1d.tau.boot, 
                   P.list = P.list, P.1d.list = P.1d.list,
                   vectors.list = vectors.list, u.1d.list = u.1d.list)
    }
  }
  
  return(r1)
}

set.seed(13)
nboot <- 5

t1 <- Sys.time()
blah <- fit.boot.Efron(model = m1, nboot = nboot, index = target,
                       dim = foo2$bic, data = data, amat = amat, newdata = fred,
                       modmat.new = modmat.renew, renewdata = renewdata, criterion = "BIC", alpha = 0.05, fit.name = NULL,
                       method = c("eigen"), quiet = FALSE)


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
  
  ######
  c1 <- makeCluster(numCores)
  registerDoParallel(c1)

  r1 <- foreach(iboot= 1:nboot, .combine = "cbind", .packages = c("envlpaster", "aster2")) %dopar% {
  # generate many new resamples    
  #for(j in 1:nboot2){
    
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
  out <- list(sd.Efron = r1, cov = cov, V = V, 
              MLE.tau.boot.subsample = MLE.tau.boot.subsample, 
              est.env.subsample = est.env.subsample)
  return(out)
}

nboot2 <- 500
set.seed(13)
internal <- function(k){
  set.seed(13)
  bar <- secondboot(k, out = blah, model = m1, data = data,
                    nboot2 = nboot2, index = target, newdata = fred,
                    amat = amat, method = "eigen")
  return(bar$sd.Efron)
}