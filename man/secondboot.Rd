% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/secondboot.R
\name{secondboot}
\alias{secondboot}
\title{secondboot}
\usage{
secondboot(
  k,
  nboot2,
  out,
  model,
  index,
  data,
  amat,
  newdata,
  method = c("eigen", "1d")
)
}
\arguments{
\item{k}{The index of the top level parametric bootstrap procedure conducted by fit.boot.Efron that the second level of bootstrapping is being applied to.}

\item{nboot2}{The bootstrap sample size for the second level of parametric bootstrapping.}

\item{out}{The output of fit.boot.Efron.}

\item{model}{An aster model object.}

\item{index}{The indices denoting which components of the canonical parameter vector are parameters of interest.}

\item{data}{An asterdata object corresponding to the original data.}

\item{amat}{This object can either be an array or a matrix. It specifies a linear combination of mean-value parameters that correspond to expected Darwinian fitness. See the aster function help page in the original aster package for more details.}

\item{newdata}{A dataframe corresponding to hypothetical individuals in which expected Darwinian fitness is to be estimated.}

\item{method}{The procedure used to obtain envelope estimators.}
}
\value{
\item{sd.Efron}{The estimated standard deviation (sd) for estimated expected Darwinian fitness where is estimation is conducted using envelope methodology. This sd accounts for model selection volatility. An eigenvalue decomposition using eigen is used internally to calculate this quantity.}
  \item{cov}{A components needed to construct sd.Efron if other numerical methods are desired.}
  \item{V}{A components needed to construct sd.Efron if other numerical methods are desired.}
  \item{MLE.tau.boot.subsample}{A components needed to construct sd.Efron if other numerical methods are desired.}
  \item{est.env.subsample}{A components needed to construct sd.Efron if other numerical methods are desired.}
}
\description{
\loadmathjax
  A parametric bootstrap procedure evaluated at an envelope estimator of the submodel mean-value parameter vector \mjeqn{\tau}{ascii} that was obtained using eigenstructures or the 1d algorithm.
}
\details{
This function implements the second level of the parametric bootstrap
procedure given by either Algorithm 1 or Algorithm 2 in Eck (2015) with
respect to the mean-value parameterization. This is detailed in Steps 4
through 5c in the algorithm below. At iteration \mjeqn{b}{ascii}, this parametric
bootstrap generates resamples from the distribution evaluated at the
envelope estimator (\mjeqn{\hat{\tau}_{env}^{(b)}}{ascii}) of \mjeqn{\tau}{ascii}.
In this case, the selected indices producing the eigenstructure which was used to
construct the envelope estimator \mjeqn{\hat{\tau}_{env}^{(b)}}{ascii} are used to
construct envelope estimators for the generated data. These resamples
are used to estimate the variability of \mjeqn{\hat{\tau}_{env}^{(b)}}{ascii}.
The algorithm using eigenstructures is as follows:

  \enumerate{
    \item Fit the aster model to the data and obtain
    \mjeqn{\hat{\tau} = (\hat{\gamma}^T, \hat{\upsilon}^T)}{ascii} and \mjeqn{\hat{\Sigma}}{ascii}
    from the aster model fit.
    \item Compute the envelope estimator of \mjeqn{\upsilon}{ascii} in the original sample, given as
    \mjeqn{\hat{\upsilon}_{env} = P_{\hat{G}}\hat{\upsilon}}{ascii} where \mjeqn{P_{\hat{G}}}{ascii} is computed using eigenstructures
    and selected via a model selection criterion of choice.
    \item Perform a parametric bootstrap by generating resamples from the distribution of the
    aster submodel evaluated at \mjeqn{\hat{\tau}_{env} = (\hat{\gamma}^T,\hat{\upsilon}_{env}^T)^T}{ascii}. For iteration
    \mjeqn{b=1,...,B}{ascii} of the procedure:
      \enumerate{
        \item Compute \mjeqn{\hat{\tau}^{(b)}}{ascii} and \mjeqn{\widehat{\Sigma}_{\upsilon,\upsilon}^{(b)}}{ascii} from the aster
        model fit to the resampled data.
        \item Build \mjeqn{P_{\hat{G}}^{(b)}}{ascii} using the indices of \mjeqn{\hat{\Sigma}_{\upsilon,\upsilon}^{(b)}}{ascii} that
        are selected using the same model selection criterion as Step 2 to build \mjeqn{\hat{G}}{ascii}.
        \item Compute \mjeqn{\hat{\upsilon}_{env}^{(b)} = P_{\hat{\mathcal{E}}}^{(b)}\hat{\upsilon}^{(b)}}{ascii}
        and
        \item Store \mjeqn{\hat{\tau}_{env}^{(b)}}{ascii} and \mjeqn{g\left(\hat{\tau}_{env}^{(b)}\right)}{ascii}
        where \mjeqn{g}{ascii} maps \mjeqn{\tau}{ascii} to the parameterization of Darwinian fitness.
      }
    \item After \mjeqn{B}{ascii} steps, the bootstrap estimator of expected Darwinian fitness is the
    average of the envelope estimators stored in Step 3d. This completes the first part of the
    bootstrap procedure.
    \item We now proceed with the second level of bootstrapping at the \mjeqn{b^{th}}{ascii} stored envelope estimator
    \mjeqn{\hat{\tau}_{env}^{(b)}}{ascii}. For iteration \mjeqn{k=1,...,K}{ascii} of the procedure:
      \enumerate{
        \item Generate data from the distribution of the aster submodel evaluated at \mjeqn{\hat{\tau}_{env}^{(b)}}{ascii}.
        \item Perform Steps 3a through 3d with respect to the dataset obtained in Step 5a.
        \item Store \mjeqn{\hat{\tau}_{env}^{(b)^{(k)}}}{ascii} and \mjeqn{g\left(\hat{\tau}_{env}^{(b)^{(k)}}\right)}{ascii}.
      }
  }

When the second level of bootstrapping is completed for all \mjeqn{b = 1,...,B}{ascii} then this
function reports the standard deviation of the bootstrapped envelope estimator of
expected Darwinian fitness. In this case, the bootstrap procedure accounts for model
selection volatility. The bootstrapped envelope estimator is

\mjdeqn{\hat{\mu}_g = \frac{1}{B} \sum_{b=1}^B g(\hat{\tau}_{env}^{(b)})}{ascii}

where \mjeqn{g(\hat{\tau}_{env}^{(b)})}{ascii} are the stored envelope estimators of expected Darwinian
fitness in the \code{env.boot.out} matrix included in the output of \code{fit.boot.Efron}.
The standard deviation of the bootstrapped envelope estimator of expected Darwinian fitness is

\mjdeqn{ \sum_{b=1}^B\left[\widehat{cov}^{(b)^T}\hat{V}^{-1}\widehat{cov}^{(b)}\right] / B }{ascii}

where \mjeqn{\widehat{cov}^{(b)} = \textbf{B}^{(b)^T} C^{(b)} / K}{ascii} and \mjeqn{\hat{V} = \textbf{B}^{(b)^T}\textbf{B}^{(b)}/K}{ascii}. The matrix \mjeqn{\textbf{B}^{(b)} \in R^{K\times p}}{ascii} has rows given by

\mjdeqn{\hat{\tau}_{env}^{(b)^{(k)}} - \sum_{k=1}^K\hat{\tau}_{env}^{(b)^{(k)}}/K}{ascii}

and the matrix \eqn{C^{(b)} \in R^{K \times d}} has columns given by

\mjdeqn{g\left(\tau_{env}^{(b)^{(k)}}\right) - g\left(\tau_{env}^{(b)}\right)}{ascii}.

For more details, see Efron (2014) and Eck (2015). The parametric bootstrap
procedure which uses the 1d algorithm to construct envelope estimators is
analogous to the above algorithm. To use the 1d algorithm, the user
specifies \code{method = "1d"} instead of \code{method = "eigen"}.
}
\examples{
### Web-based Supplementary Materials for ``Enveloping the aster model.'' ###
}
\references{
Cook, R.D. and Zhang, X. (2014). Foundations for Envelope Models and Methods. \emph{JASA}, In Press.\cr
\cr
Cook, R.D. and Zhang, X. (2015). Algorithms for Envelope Estimation. \emph{Journal of Computational and Graphical Statistics}, Published online. \doi{10.1080/10618600.2015.1029577}.\cr
\cr
Eck, D. J., Geyer, C. J., and Cook, R. D. (2016). Enveloping the aster model. \emph{in prep}. \cr
\cr
Eck, D.~J., Geyer, C.~J., and Cook, R.~D. (2016). Web-based Supplementary Materials for ``Enveloping the aster model.'' \emph{in prep}. \cr
\cr
Efron, B. (2014). Estimation and Accuracy After Model Selection. \emph{JASA}, \strong{109:507}, 991-1007.\cr
}
