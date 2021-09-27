\name{fit.boot.Efron}
\alias{fit.boot.Efron}
\title{fit.boot.Efron}
\usage{
  fit.boot.Efron(model, nboot, index, vectors = NULL, dim = NULL,
    data, amat, newdata, modmat.new = NULL, renewdata = NULL, 
    criterion = c("AIC","BIC","LRT"), alpha = 0.05, fit.name = NULL, 
    method = c("eigen","1d"), quiet = FALSE)
}

\description{
  A parametric bootstrap procedure evaluated at an envelope estimator 
  of the submodel mean-value parameter vector \eqn{\tau} that was 
  obtained using reducing subspaces or the 1d algorithm.  
}

\arguments{
  \item{model}{An aster model object.}

  \item{nboot}{The number of bootstrap iterations desired.}

  \item{index}{The indices denoting which components of the 
    canonical parameter vector are parameters of interest.}

  \item{vectors}{The indices denoting which reducing subspace  
    of Fisher information is desired to construct 
    envelope estimators. Must be specified if \code{method = "eigen"}.}

  \item{dim}{The dimension of the envelope space used to construct envelope 
    estimators. Must be specified if \code{method = "1d"}. }  

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

  \item{modmat.new}{A model matrix corresponding to hypothetical 
    individuals in which expected Darwinian fitness is to 
    be estimated.}

  \item{renewdata}{A dataframe in long format corresponding to 
    hypothetical individuals in which expected Darwinian fitness 
    is to be estimated.}  

  \item{criterion}{A model selection criterion of choice.}  

  \item{alpha}{The type 1 error rate desired for the LRT.}

  \item{fit.name}{An expression that appears in the name of the 
    nodes that correspond to Darwinian fitness. This is only
    necessary if \code{renewdata} is not provided.}

  \item{method}{The procedure used to obtain envelope estimators.}    

  \item{quiet}{A logical argument. If FALSE, the function 
    displays how much time it takes to run \code{m} iterations.}

}

\value{


  \item{env.boot.out}{Estimated expected Darwinian fitness from generated 
    data obtained from Steps 3a-3d in the bootstrap procedure using the envelope 
    estimator constructed using reducing subspaces. }
    
  \item{MLE.boot.out}{Estimated expected Darwinian fitness from generated 
    data obtained from Steps 3a-3d in the bootstrap procedure using the MLE. }

  \item{env.1d.boot.out}{Estimated expected Darwinian fitness from generated 
    data obtained from Steps 3a-3d in the bootstrap procedure using the envelope 
    estimator constructed using the 1d algorithm. }

  \item{env.tau.boot}{Estimated mean-value parameter vectors from generated 
    data obtained from Steps 3a-3d in the bootstrap procedure using the envelope 
    estimator constructed using reducing subspaces. }
    
  \item{MLE.tau.boot}{Estimated mean-value parameter vectors from generated 
    data obtained from Steps 3a-3d in the bootstrap procedure using the MLE. }

  \item{env.1d.tau.boot}{Estimated mean-value parameter vectors from generated 
    data obtained from Steps 3a-3d in the bootstrap procedure using the envelope 
    estimator constructed using the 1d algorithm. }

  \item{P.list}{A list of all estimated projections into the envelope 
    space constructed from reducing subspaces for Steps 3a-3d in the bootstrap 
    procedure. }  

  \item{P.1d.list}{A list of all estimated projections into the envelope 
    space constructed using the 1d algorithm for Steps 3a-3d in the bootstrap 
    procedure. }     

  \item{vectors.list}{A list of indices of eigenvectors used to build the 
    projections in P.list. These indices are selected using the user specified 
    model selection criterion as indicated in Steps 3a-3d in the bootstrap 
    procedure. }  

  \item{u.1d.list}{A list of indices of eigenvectors used to build the 
    projections in P.list. These indices are selected using the user specified 
    model selection criterion as indicated in Steps 3a-3d in the bootstrap 
    procedure. }     

}

\details{

  This function implements the first level of the parametric bootstrap 
  procedure given by either Algorithm 1 or Algorithm 2 in Eck (2015) with 
  respect to the mean-value parameterization. This is detailed in Steps 1 
  through 3d in the algorithm below.  This parametric bootstrap  generates 
  resamples from the distribution evaluated at an envelope estimator of \eqn{\tau} 
  adjusting for model selection volatility. 

  The user specifies a model selection criterion which selects vectors that 
  construct envelope estimators using the reducing subspace approach. The user also 
  specifies which method is to be used in order to calculate envelope 
  estimators. When one is using a partial envelope, then this function 
  constructs envelope estimators of \eqn{\upsilon} where we write \eqn{\tau} = \eqn{(\gamma^T,\upsilon^T)^T} 
  and \eqn{\upsilon} corresponds to aster model parameters of interest. 
  In applications, candidate reducing subspaces are indices of eigenvectors of \eqn{\widehat{\Sigma}_{\upsilon,\upsilon}}
  where \eqn{\widehat{\Sigma}_{\upsilon,\upsilon}} is the part of \eqn{\hat{\Sigma}} 
  corresponding to our parameters of interest. These indices are specified 
  by \code{vectors}. When all of the components of \eqn{\tau} are components 
  of interest, then we write \eqn{\widehat{\Sigma}_{\upsilon,\upsilon} = \widehat{\Sigma}}. When data 
  is generated via the parametric bootstrap, it is the indices (not the 
  original reducing subspaces) that are used to construct envelope estimators 
  constructed using the generated data. The algorithm using reducing subspaces 
  is as follows: 

  \enumerate{
  \item[1.] Fit the aster model to the data and obtain 
    \eqn{\hat{\tau} = (\hat{\gamma}^T, \hat{\upsilon}^T)} and \eqn{\hat{\Sigma}} 
    from the aster model fit.
  \item[2.] Compute the envelope estimator of \eqn{\upsilon} in the original sample, given as 
    \eqn{ \hat{\upsilon}_{env} = P_{\hat{G}}\hat{\upsilon}} where \eqn{P_{\hat{G}}} is computed using reducing subspaces 
    and selected via a model selection criterion of choice.
  \item[3.] Perform a parametric bootstrap by generating resamples from the distribution of the 
    aster submodel evaluated at \eqn{\hat{\tau}_{env} = (\hat{\gamma}^T,\hat{\upsilon}_{env}^T)^T}. For iteration 
    \eqn{b=1,...,B} of the procedure: 
    \enumerate{
      \item[(3a)] Compute \eqn{\hat{\tau}^{(b)}} and \eqn{\widehat{\Sigma}_{\upsilon,\upsilon}^{(b)}} from the aster 
        model fit to the resampled data.
      \item[(3b)] Build \eqn{P_{\hat{G}}^{(b)}} using the indices of \eqn{\hat{\Sigma}_{\upsilon,\upsilon}^{(b)}} that 
        are selected using the same model selection criterion as Step 2 to build \eqn{\hat{G}}.
      \item[(3c)] Compute \eqn{\hat{\upsilon}_{env}^{(b)} = P_{\hat{\mathcal{E}}}^{(b)}\hat{\upsilon}^{(b)}} and 
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


  The parametric bootstrap procedure which uses the 1d algorithm to construct 
  envelope estimators is analogous to the above algorithm. To use the 1d 
  algorithm, the user specifies \code{method = "1d"}. A parametric bootstrap 
  generating resamples from the distribution evaluated at the aster model 
  MLE is also conducted by this function. 

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
### see Web-based Supplementary Materials for ``Enveloping the aster model.'' ###
}
