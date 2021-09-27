\name{targetboot}
\alias{targetboot}
\title{targetboot}
\usage{
  targetboot(model, nboot, index, u, data, 
    quiet = FALSE, m = 100)
}

\description{
  A parametric bootstrap procedure evaluated at an envelope estimator 
  of the submodel mean-value parameter vector \eqn{\tau} that was 
  obtained using the 1D algorithm.
}

\arguments{
  \item{model}{An aster model object}

  \item{nboot}{The number of bootstrap iterations desired}

  \item{index}{The indices denoting which components of the 
    canonical parameter vector are parameters of interest}

  \item{u}{The dimension of the envelope space assumed}

  \item{data}{An asterdata object}

  \item{quiet}{A logical argument. If FALSE, the function 
    displays how much time it takes to run \code{m} iterations}

  \item{m}{The length of the output interval}
}

\value{
  \item{u}{The dimension of the envelope space assumed}

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
  the algorithm below with respect to the mean-value parameterization. 
  This parametric bootstrap generates resamples from the distribution 
  evaluated at an envelope estimator of \eqn{\tau}. Envelope estimators are 
  constructed using the 1D algorithm at a user-specified envelope model 
  dimension \code{u}. When one is using a partial envelope then this function 
  constructs envelope estimators of \eqn{\upsilon} where we write 
  \eqn{\tau} = \eqn{(\gamma^T,\upsilon^T)^T} and \eqn{\upsilon} 
  corresponds to aster model parameters of interest. In the sample, the 
  1D algorithm uses \eqn{M = \widehat{\Sigma}_{\upsilon,\upsilon}} and 
  \eqn{U = \hat{\tau}\hat{\tau}^T} as inputs where \eqn{\widehat{\Sigma}_{\upsilon,\upsilon}} 
  is the part of \eqn{\widehat{\Sigma}} corresponding to our parameters 
  of interest. When all of the components of \eqn{\tau} are components 
  of interest, then we write \eqn{\widehat{\Sigma}_{\upsilon,\upsilon} = \widehat{\Sigma}}. 
  The algorithm is as follows: 


  \enumerate{
  \item[1.] Fit the aster model to the data and obtain 
    \eqn{\hat{\tau} = (\hat{\gamma}^T, \hat{\upsilon}^T)} and \eqn{\widehat{\Sigma}} 
    from the aster model fit.
  \item[2.] Compute the envelope estimator of \eqn{\upsilon} in the original 
    sample, given as \eqn{\hat{\upsilon}_{env} = P_{\hat{\mathcal{E}}}\hat{\upsilon}} 
    where \eqn{P_{\hat{\mathcal{E}}}} is obtained from the 1D algorithm.
  \item[3.] Perform a parametric bootstrap by generating resamples from the 
  distribution evaluated at \eqn{\hat{\upsilon}_{env}}. For iteration 
  \eqn{b=1,...,B} of the procedure: 
    \enumerate{
      \item[(3a)] Compute \eqn{\hat{\upsilon}^{(b)}} and 
        \eqn{\widehat{\Sigma}_{\upsilon,\upsilon}^{(b)}} from the aster model fit to 
        the resampled data.
      \item[(3b)] Obtain \eqn{P_{\hat{\mathcal{E}}}^{(b)}} as done in Step 2. 
      \item[(3c)] Store \eqn{\hat{\upsilon}_{env}^{(b)} = P_{\hat{\mathcal{E}}}^{(b)}\hat{\upsilon}^{(b)}}.
    }
  }

  A parametric bootstrap generating resamples from the distribution 
  evaluated at the aster model MLE is also conducted by this function.    
}

\references{
  Cook, R.D. and Zhang, X. (2014).
  Foundations for Envelope Models and Methods.
  \emph{JASA}, In Press.

  Cook, R.D. and Zhang, X. (2015).
  Algorithms for Envelope Estimation.
  \emph{Journal of Computational and Graphical Statistics}, 
  Published online. \doi{10.1080/10618600.2015.1029577}.

  Eck, D. J., Geyer, C. J., and Cook, R. D. (2016).
  Enveloping the aster model.
  \emph{in prep}.
}

\examples{
\dontrun{set.seed(13)
library(envlpaster)
library(aster2)
data(simdata30nodes)
data <- simdata30nodes.asterdata
nnode <- length(vars)
xnew <- as.matrix(simdata30nodes[,c(1:nnode)])
m1 <- aster(xnew, root, pred, fam, modmat)
target <- 5:9
indices <- c(1,2,4,5)
u <- length(indices)
nboot <- 2000; timer <- nboot/2
bar <- targetboot(m1, nboot = nboot, index = target, 
  u = u, data = data, m = timer)
bar}
}
