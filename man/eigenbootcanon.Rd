\name{eigenbootcanon}
\alias{eigenbootcanon}
\title{eigenbootcanon}
\usage{
  eigenbootcanon(model, nboot, index, vectors, 
    code, families, quiet = FALSE, m = 100)
}

\description{
  A parametric bootstrap procedure evaluated at an envelope estimator 
  of the submodel mean-value parameter vector \eqn{\beta} that was 
  obtained using reducing subspaces.
}

\arguments{
  \item{model}{An aster model object.}

  \item{nboot}{The number of bootstrap iterations desired.}

  \item{index}{The indices denoting which components of the 
    canonical parameter vector are parameters of interest.}

  \item{vectors}{The indices denoting which reducing subspace 
    of inverse Fisher information is desired to construct 
    envelope estimators.}

  \item{code}{A vector of dimension equal to the number 
    of nodes in the aster graph. This vector specifies 
    which exponential family in the \code{families} list
    labels each arrow of the aster graphical structure.}

  \item{families}{A list of family specifications.}  

  \item{quiet}{A logical argument. If FALSE, the function 
    displays how much time it takes to run \code{m} iterations.}

  \item{m}{The length of the output interval.}
}

\value{
  \item{u}{The dimension of the envelope space assumed.}

  \item{table}{A table of output. The first two columns display 
    the envelope estimator and its bootstrapped standard error. 
    The next two columns display the MLE and its bootstrapped 
    standard error. The last column displays the ratio of the 
    standard error for the bootstrapped envelope estimator 
    to the standard error for the bootstrapped MLE.}

  \item{S}{The bootstrap estimator of the variability of 
    the partial envelope estimator.}

  \item{S2}{The bootstrap estimator of the variability of 
    the MLE.}

  \item{env.boot.out}{The realizations from the bootstrap procedure 
    using envelope methodology.}
    
  \item{MLE.boot.out}{The realizations from the bootstrap procedure 
    using maximum likelihood estimation.}
}

\details{
  This function implements the parametric bootstrap procedure given by 
  the algorithm below with respect to the canonical 
  parameterization. This parametric bootstrap generates resamples from the 
  distribution evaluated at an envelope estimator of \eqn{\beta}. The user 
  specifies which vectors are used in order to construct envelope estimators 
  using the reducing subspace approach. When one is using a partial envelope then 
  this function constructs envelope estimators of \eqn{\upsilon} where we 
  write \eqn{\beta} = \eqn{(\gamma^T,\upsilon^T)^T} and \eqn{\upsilon} 
  corresponds to aster model parameters of interest. In applications, candidate 
  reducing subspaces are indices of eigenvectors of \eqn{\widehat{\Sigma}_{\upsilon,\upsilon}^{-1}} 
  where \eqn{\widehat{\Sigma}_{\upsilon,\upsilon}^{-1}} is the part of \eqn{\widehat{\Sigma}^{-1}} 
  corresponding to our parameters of interest. When all of the components 
  of \eqn{\beta} are components of interest, then we write 
  \eqn{\widehat{\Sigma}_{\upsilon,\upsilon}^{-1} = \widehat{\Sigma}^{-1}}. The algorithm is as follows: 

  \enumerate{
  \item[1.] Fit the aster model to the data and obtain 
    \eqn{\hat{\beta} = (\hat{\gamma}^T, \hat{\upsilon}^T)} and 
    \eqn{\widehat{\Sigma}^{-1}} from the aster model fit.
  \item[2.] Compute the envelope estimator of \eqn{\upsilon} in the original 
    sample, given as \eqn{\hat{\upsilon}_{env} = P_{\hat{G}}\hat{\upsilon}} 
    where \eqn{P_{\hat{G}}} is the projection into the reducing subspace of 
    \eqn{\widehat{\Sigma}_{\upsilon,\upsilon}^{-1}} specified by \code{vectors}. 
  \item[3.] Perform a parametric bootstrap by generating resamples from the 
  distribution evaluated at \eqn{\hat{\upsilon}_{env}}. For iteration 
  \eqn{b=1,...,B} of the procedure: 
    \enumerate{
      \item[(3a)] Compute \eqn{\hat{\upsilon}^{(b)}} and 
        \eqn{\widehat{\Sigma}_{\upsilon,\upsilon}^{(b)^{-1}}} from the aster model fit to 
        the resampled data.
      \item[(3b)] Obtain \eqn{P_{\hat{G}}^{(b)}} as done in Step 2. 
      \item[(3c)] Store \eqn{\hat{\upsilon}_{env}^{(b)} = P_{\hat{G}}^{(b)}\hat{\upsilon}^{(b)}}.
    }
  }

  A parametric bootstrap generating resamples from the distribution 
  evaluated at the aster model MLE is also conducted by this function.  
}

\references{
  Eck, D. J., Geyer, C. J., and Cook, R. D. (2016).
  Enveloping the aster model.
  \eqn{\emph{in prep}}.
}

\examples{
\dontrun{set.seed(13)
library(envlpaster)
library(aster2)
data(generateddata)
m1 <- aster(resp ~ 0 + varb + mass + timing, 
  fam = fam, pred = pred, varvar = varb, idvar = id,
  root = root, data = redata)
target <- c(9:10)
nboot <- 2000; timer <- nboot/2
bar <- eigenbootcanon(m1, nboot = nboot, index = target,
  vectors = c(2), u = 1, m = timer)
bar}
}
