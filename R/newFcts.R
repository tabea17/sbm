
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


rnsbmObs <- function(theta, Z, modelFamily='Gauss', directed=FALSE){   # Z matrice n,Q
  n=nrow(Z)

  # adjacency matrix
  wqlGrand <- Z %*% wVectToMat(theta$w) %*% t(Z)
  A = rbinom(n^2, size=1, prob=wqlGrand) %>% matrix(n, n)
  diag(A) <- 0
  if (!directed) A <- A * lower.tri(A) + t(A * lower.tri(A))

  # noisy observations under the null
    nu0Grand_1 <-  matrix(theta$nu0[1], n, n)
    nuqlGrand_1 <- Z %*% wVectToMat(theta$nu[,1]) %*% t(Z)
    nuGrand_1 <- nu0Grand_1 * (A==0) +  nuqlGrand_1 * (A==1)

    if (modelFamily!= 'Poisson'){
       nu0Grand_2 <-  matrix(theta$nu0[2], n, n)
       nuqlGrand_2 <- Z %*% wVectToMat(theta$nu[,2]) %*% t(Z)
       nuGrand_2 <- nu0Grand_2 * (A==0) +  nuqlGrand_2 * (A==1)
    }

  if (modelFamily=='Gauss'){
    X <- stats::rnorm(n^2, nuGrand_1, nuGrand_2)
  }
  if (modelFamily=='Gamma'){
    X <- stats::rgamma(n^2, nuGrand_1, nuGrand_2)
  }
  if (modelFamily=='Poisson'){
    X <- stats::rpois(n^2, nuGrand_1)
  }

  diag(X) <- 0
  if (!directed) X <- X * lower.tri(X) + t(X * lower.tri(X))
  return(list(dataMatrix=X, latentAdj=A))
}



