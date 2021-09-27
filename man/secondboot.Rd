\name{secondboot}
\alias{secondboot}
\title{secondboot}
\usage{
  secondboot(k, nboot2, out, model, index, data, amat, 
    newdata, method = c("eigen","1d"))
}

\description{
  A parametric bootstrap procedure evaluated at an envelope estimator 
  of the submodel mean-value parameter vector \eqn{\tau} that was 
  obtained using eigenstructures or the 1d algorithm.  
}

\arguments{

  \item{k}{The index of the top level parametric bootstrap procedure conducted by \code{fit.boot.Efron}
    that the second level of bootstrapping is being applied to. }

  \item{nboot2}{The bootstrap sample size for the second level of parametric bootstrapping.}

  \item{out}{The output of \code{fit.boot.Efron}.}

  \item{model}{An aster model object.}

  \item{index}{The indices denoting which components of the 
    canonical parameter vector are parameters of interest.}

  \item{data}{An asterdata object corresponding to the 
    original data.}

  \item{amat}{This object can either be an array or a matrix. 
    It specifies a linear combination of mean-value parameters 
    that correspond to expected Darwinian fitness. See the 
    \code{aster} function help page in the original \code{aster} 
    package for more details.}

  \item{newdata}{A dataframe corresponding to hypothetical 
    individuals in which expected Darwinian fitness is to 
    be estimated.}    

  \item{method}{The procedure used to obtain envelope estimators.}    

}

\value{

  \item{sd.Efron}{The estimated standard deviation (sd) for estimated expected
    Darwinian fitness where is estimation is conducted using envelope
    methodology. This sd accounts for model selection volatility. An eigenvalue 
    decomposition using \code{eigen} is used internally to calculate this quantity. }

  \item{cov}{A components needed to construct \code{sd.Efron} if other numerical 
    methods are desired. }

  \item{V}{A components needed to construct \code{sd.Efron} if other numerical 
    methods are desired. }

  \item{MLE.tau.boot.subsample}{A components needed to construct \code{sd.Efron} if 
    other numerical methods are desired. }

  \item{est.env.subsample}{A components needed to construct \code{sd.Efron} if other 
    numerical methods are desired. }  

}

\details{

  This function implements the second level of the parametric bootstrap 
  procedure given by either Algorithm 1 or Algorithm 2 in Eck (2015) with 
  respect to the mean-value parameterization. This is detailed in Steps 4 
  through 5c in the algorithm below. At iteration \eqn{b}, this parametric 
  bootstrap generates resamples from the distribution evaluated at the 
  envelope estimator (\eqn{\hat{\tau}_{env}^{(b)}}) of \eqn{\tau}. 
  In this case, the selected indices producing the eigenstructure which was used to 
  construct the envelope estimator \eqn{\hat{\tau}_{env}^{(b)}} are used to
  construct envelope estimators for the generated data. These resamples 
  are used to estimate the variability of \eqn{\hat{\tau}_{env}^{(b)}}. 
  The algorithm using eigenstructures is as follows:

  \enumerate{
  \item[1.] Fit the aster model to the data and obtain 
    \eqn{\hat{\tau} = (\hat{\gamma}^T, \hat{\upsilon}^T)} and \eqn{\hat{\Sigma}} 
    from the aster model fit.
  \item[2.] Compute the envelope estimator of \eqn{\upsilon} in the original sample, given as 
    \eqn{\hat{\upsilon}_{env} = P_{\hat{G}}\hat{\upsilon}} where \eqn{P_{\hat{G}}} is computed using eigenstructures 
    and selected via a model selection criterion of choice.
  \item[3.] Perform a parametric bootstrap by generating resamples from the distribution of the 
    aster submodel evaluated at \eqn{\hat{\tau}_{env} = (\hat{\gamma}^T,\hat{\upsilon}_{env}^T)^T}. For iteration 
    \eqn{b=1,...,B} of the procedure: 
    \enumerate{
      \item[(3a)] Compute \eqn{\hat{\tau}^{(b)}} and \eqn{\widehat{\Sigma}_{\upsilon,\upsilon}^{(b)}} from the aster 
        model fit to the resampled data.
      \item[(3b)] Build \eqn{P_{\hat{G}}^{(b)}} using the indices of \eqn{\hat{\Sigma}_{\upsilon,\upsilon}^{(b)}} that 
        are selected using the same model selection criterion as Step 2 to build \eqn{\hat{G}}.
      \item[(3c)] Compute \eqn{\hat{\upsilon}_{env}^{(b)} = P_{\hat{\mathcal{E}}}^{(b)}\hat{\upsilon}^{(b)}} 
        and  
        \eqn{\hat{\tau}_{env}^{(b)} = \left(\hat{\gamma}^{(b)^T},\hat{\upsilon}_{env}^{(b)^T}\right)^T}.
      \item[(3d)] Store \eqn{\hat{\tau}_{env}^{(b)}} and \eqn{g\left(\hat{\tau}_{env}^{(b)}\right)} 
        where \eqn{g} maps \eqn{\tau} to the parameterization of Darwinian fitness.
    }
  \item[4.] After \eqn{B} steps, the bootstrap estimator of expected Darwinian fitness is the 
    average of the envelope estimators stored in Step 3d. This completes the first part of the 
    bootstrap procedure. 
  \item[5.] We now proceed with the second level of bootstrapping at the \eqn{b^{th}} stored envelope estimator 
    \eqn{\hat{\tau}_{env}^{(b)}}. For iteration \eqn{k=1,...,K} of the procedure:
    \enumerate{
      \item[(5a)] Generate data from the distribution of the aster submodel evaluated at \eqn{\hat{\tau}_{env}^{(b)}}.
      \item[(5b)] Perform Steps 3a through 3d with respect to the dataset obtained in Step 5a.
      \item[(5c)] Store \eqn{\hat{\tau}_{env}^{(b)^{(k)}}} and \eqn{g\left(\hat{\tau}_{env}^{(b)^{(k)}}\right)}.
    }
  }

  When the second level of bootstrapping is completed for all \eqn{b = 1,...,B} then this 
  function reports the standard deviation of the bootstrapped envelope estimator of 
  expected Darwinian fitness. In this case, the bootstrap procedure accounts for model 
  selection volatility. The bootstrapped envelope estimator is 

  \deqn{\hat{\mu}_g = \frac{1}{B} \sum_{b=1}^B g(\hat{\tau}_{env}^{(b)})}
  
  where \eqn{g(\hat{\tau}_{env}^{(b)})} are the stored envelope estimators of expected Darwinian 
  fitness in the \code{env.boot.out} matrix included in the output of \code{fit.boot.Efron}. 
  The standard deviation of the bootstrapped envelope estimator of expected Darwinian fitness is 

  \deqn{ \sum_{b=1}^B\left[\widehat{cov}^{(b)^T}\hat{V}^{-1}\widehat{cov}^{(b)}\right] / B } 
  
  where \eqn{\widehat{cov}^{(b)} = \textbf{B}^{(b)^T} C^{(b)} / K} and \eqn{\hat{V} = \textbf{B}^{(b)^T}\textbf{B}^{(b)}/K}. The matrix \eqn{\textbf{B}^{(b)} \in R^{K\times p}} has rows given by

  \deqn{\hat{\tau}_{env}^{(b)^{(k)}} - \sum_{k=1}^K\hat{\tau}_{env}^{(b)^{(k)}}/K} 
  
  and the matrix \eqn{C^{(b)} \in R^{K \times d}} has columns given by 

  \deqn{g\left(\tau_{env}^{(b)^{(k)}}\right) - g\left(\tau_{env}^{(b)}\right)}. 
  
  For more details, see Efron (2014) and Eck (2015). The parametric bootstrap 
  procedure which uses the 1d algorithm to construct envelope estimators is 
  analogous to the above algorithm. To use the 1d algorithm, the user 
  specifies \code{method = "1d"} instead of \code{method = "eigen"}.

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
  \eqn{\emph{in prep}}.

  Eck, D.~J., Geyer, C.~J., and Cook, R.~D. (2016).
  Web-based Supplementary Materials for ``Enveloping the aster model.''
  \eqn{\emph{in prep}}.

  Efron, B. (2014).
  Estimation and Accuracy After Model Selection.
  \eqn{\emph{JASA}}, \eqn{\textbf{109:507}}, 991-1007.
}

\examples{
### Web-based Supplementary Materials for ``Enveloping the aster model.'' ###
}
