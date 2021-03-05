
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




