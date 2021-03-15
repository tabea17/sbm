
nbParamSubmodel <- function(submodel, N_Q) {
   dimW <- N_Q
   if (submodel %in% c('Gauss0', 'ExpGamma')){   ####FANNY : rajouter private$submodel
     dimH0 <- 1
     dimH1 <- 2*N_Q
   }
   if (submodel=='Gauss'){
    dimH0 <- 2
    dimH1 <- 2*N_Q
  }
   if (submodel=='Gauss01'){
    dimH0 <- 0
    dimH1 <- 2*N_Q
  }
  if (submodel=='Gauss0Var1'){
    dimH0 <- 0
    dimH1 <- N_Q
  }
  if (submodel=='GaussEqVar'){
   dimH0 <- 2 # = 1 gaussian mean under H0 + 1 common variance for all gaussians in the model
   dimH1 <- N_Q # nb of gaussian means under H1
  }
  if (submodel=='Gauss0EqVar'){
   dimH0 <- 1 # 1 common variance for all gaussians in the model
   dimH1 <- N_Q # nb of gaussian means under H1
  }
  if (submodel=='GaussAffil'){
   dimH0 <- 2
   dimH1 <- 4
  }
  if (submodel=='Gauss2distr'){
   dimH0 <- 2
   dimH1 <- 2
  }
  if (submodel=='Exp'){
   dimH0 <- 1
   dimH1 <- N_Q
  }

  dimParam <- dimW + dimH0 + dimH1
  return(dimParam)
}



wVectToMat <- function(w.vec, Q) {
  M <- matrix(0, Q, Q)
  M[upper.tri(M, diag=TRUE)]= w.vec
  wMat=M + t(M)
  diag(wMat)=diag(M)
  return(wMat)
}

#theta=list(pi=rep(1/2,2), w=c(0.8,0.1,0.9), nu0=c(0,1), nu=matrix(c(0,0,0,2,2,2), ncol=2)
# rnsbmObs <- function(theta, Z, modelFamily='Gauss', directed=FALSE){   # Z matrice n,Q
#   n=nrow(Z)
#   Q=ncol(Z)
#   # adjacency matrix
#   wqlGrand <- Z %*% wVectToMat(theta$w, Q) %*% t(Z)
#   A = rbinom(n^2, size=1, prob=wqlGrand) %>% matrix(n, n)
#   diag(A) <- 0
#   if (!directed) A <- A * lower.tri(A) + t(A * lower.tri(A))
#
#   # noisy observations under the null
#     nu0Grand_1 <-  matrix(theta$nu0[1], n, n)
#     nuqlGrand_1 <- Z %*% wVectToMat(theta$nu[,1], Q) %*% t(Z)
#     nuGrand_1 <- nu0Grand_1 * (A==0) +  nuqlGrand_1 * (A==1)
#
#     if (modelFamily!= 'Poisson'){
#        nu0Grand_2 <-  matrix(theta$nu0[2], n, n)
#        nuqlGrand_2 <- Z %*% wVectToMat(theta$nu[,2], Q) %*% t(Z)
#        nuGrand_2 <- nu0Grand_2 * (A==0) +  nuqlGrand_2 * (A==1)
#     }
#
#   if (modelFamily=='Gauss'){
#     X <- stats::rnorm(n^2, nuGrand_1, nuGrand_2) %>% matrix(n, n)
#   }
#   if (modelFamily=='Gamma'){
#     X <- stats::rgamma(n^2, nuGrand_1, nuGrand_2) %>% matrix(n, n)
#   }
#   if (modelFamily=='Poisson'){
#     X <- stats::rpois(n^2, nuGrand_1) %>% matrix(n, n)
#   }
#
#   diag(X) <- 0
#   if (!directed) X <- X * lower.tri(X) + t(X * lower.tri(X))
#   return(list(dataMatrix=X, latentAdj=A))
# }


rnsbmObs <- function(parametre, Z, modelFamily='Gauss', directed=FALSE){   # Z matrice n,Q; parametre$w : matrice
  n=nrow(Z)
  Q=ncol(Z)
  # adjacency matrix
  wqlGrand <- Z %*% parametre$w %*% t(Z)
  A = rbinom(n^2, size=1, prob=wqlGrand) %>% matrix(n, n)
  diag(A) <- 0
  if (!directed) A <- A * lower.tri(A) + t(A * lower.tri(A))

  # noisy observations under the null
  nu0Grand_1 <-  matrix(parametre$nu0[1], n, n)
  nuqlGrand_1 <- Z %*% wVectToMat(parametre$nu[,1], Q) %*% t(Z)
  nuGrand_1 <- nu0Grand_1 * (A==0) +  nuqlGrand_1 * (A==1)

  if (modelFamily!= 'Poisson'){
    nu0Grand_2 <-  matrix(parametre$nu0[2], n, n)
    nuqlGrand_2 <- Z %*% wVectToMat(parametre$nu[,2], Q) %*% t(Z)
    nuGrand_2 <- nu0Grand_2 * (A==0) +  nuqlGrand_2 * (A==1)
  }

  if (modelFamily=='Gauss'){
    X <- stats::rnorm(n^2, nuGrand_1, nuGrand_2) %>% matrix(n, n)
  }
  if (modelFamily=='Gamma'){
    X <- stats::rgamma(n^2, nuGrand_1, nuGrand_2) %>% matrix(n, n)
  }
  if (modelFamily=='Poisson'){
    X <- stats::rpois(n^2, nuGrand_1) %>% matrix(n, n)
  }

  diag(X) <- 0
  if (!directed) X <- X * lower.tri(X) + t(X * lower.tri(X))
  return(list(dataMatrix=X, latentAdj=A))
}



graphInferenceFromqvalues <- function(qvalues, n, alpha=0.05){

  A <- matrix(0, n, n)
  A[lower.tri(A)] <- (qval_results<alpha)
  A <- A + t(A)
  return(A)
}



# vient de sample.R

#' Sampling of Noisy SBMs
#'
#' This function samples a simple Stochastic Block Models, with various model
#' for the distribution of the edges:  Bernoulli, Poisson, or Gaussian models, and possibly with covariates
#'
#' @param nbNodes number of nodes in the network
#' @param blockProp parameters for block proportions
#' @param connectParam list of parameters for connectivity with a matrix of means 'mean' and an optional matrix of variances 'var', the sizes of which must match \code{blockProp} length
#' @param model character describing the model for the relation between nodes (\code{'bernoulli'}, \code{'poisson'}, \code{'gaussian'}, ...). Default is \code{'bernoulli'}.
#' @param directed logical, directed network or not. Default is \code{FALSE}.
#' @param dimLabels an optional list of labels for each dimension (in row, in column)
#' @param covariates a list of matrices with same dimension as mat describing covariates at the edge level. No covariate per Default.
#' @param covariatesParam optional vector of covariates effect. A zero length numeric vector by default.
#'
#' @return  an object with class \code{\link{SimpleSBM}}
#'
#' @examples
#' ### =======================================
#' Noisy Gaussian model
#' ## Graph parameters
#' nbNodes  <- 10
#' blockProp <- c(.5, .5)      # group proportions
#' connectParam=list(mean=matrix(c(0.8,0.1,0.1,0.9), ncol=2, nrow=2))
#' noiseParam=c(0,1)
#' signalParam=matrix(c(10,10,10,.2,.2,.2), ncol=2)
#' ## Graph Sampling
#' mySampler <- sampleNoisySBM(model = 'Gauss', nbNodes, directed = FALSE, blockProp, connectParam, noiseParam, signalParam)
#' plot(mySampler)
#' plot(mySampler,type='meso')
#' plot(mySampler,type='latentNetwork')
#' hist(mySampler$networkData)
#' @export
sampleNoisySBM <- function(model = 'Gauss',
                           nbNodes,
                           directed = FALSE,
                           blockProp,
                           connectParam,
                           noiseParam,
                           signalParam
                           ) {

  mySampler <- NoisySBM$new(model, nbNodes, directed, blockProp, connectParam, noiseParam, signalParam)
  mySampler$rNetwork(store = TRUE)
  mySampler
}
