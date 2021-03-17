#' R6 Class definition of a Noisy SBM fit
#'
#' This class is designed to give a representation and adjust an Noisy SBM fitted.
#'
#' @import R6 blockmodels
#' @export
NoisySBM_fit <-
  R6::R6Class(
    classname = "NoisySBM_fit",
    inherit = NoisySBM,
    private = list(
      qvalues_=NULL,
      testLevel_ =NULL,
      submodel=NULL,
      J              = NULL, # approximation of the log-likelihood
      vICL           = NULL, # approximation of the ICL
      noisySBMobject       = NULL, # blockmodels output (used to stored the optimization results when blockmodels is used)
      import_from_noisySBM = function(Q = noisySBM::getBestQ(private$noisySBMobject)$Q) { # a function updating the Class
         fit <- private$noisySBMobject[[Q]]
         private$J     <- fit$convergence$J
         private$vICL  <- fit$sbmParam$ICL
         private$parameters    <- fit$theta   ### verif
         private$parameters$w <- wVectToMat(fit$theta$w, Q)
         private$Z  <- t(fit$sbmParam$clusterProba)
         private$pi <- fit$theta$pi #parameters$pi
        # private$edgeProba <- fit$sbmParam$edgeProba  #N_Q x N matrix
      }
    ),
    public = list(
      #' @description constructor for a NoisySBM fit
      #' @param dataMatrix square (noisy) matrix
      #' @param submodel character (\code{'Gauss'}, \code{'Gauss01'})
      #' @param directed logical, directed network or not. In not, \code{dataMatrix} must be symmetric.
      initialize = function(dataMatrix, submodel, directed) {

        ## SANITY CHECKS (on data)
        stopifnot(is.matrix(dataMatrix))                   # must be a matrix
        stopifnot(all.equal(nrow(dataMatrix),
                            ncol(dataMatrix)))             # matrix must be square
        stopifnot(isSymmetric(dataMatrix))    # symmetry and direction must agree

        ## INITIALIZE THE SBM OBJECT ACCORDING TO THE DATA
       if (submodel %in% c('Gauss','Gauss0','Gauss01','GaussEqVar','Gauss0EqVar','Gauss0Var1','Gauss2distr','GaussAffil'))
            modelFamily <- 'Gauss'
       if (submodel %in% c('Exp','ExpGamma','Gamma'))
            modelFamily <- 'Gamma'

        super$initialize(modelFamily  = modelFamily,
                         nbNodes      = nrow(dataMatrix),
                         directed     = FALSE,
                         blockProp    = vector("numeric", 0),
                         connectParam =  list(mean = matrix(0, 0, 0)), #connectParam,
                         dimLabels    =  c(node="nodeName"), #dimLabels,
                         signalParam = 0,
                         noiseParam = 0,
                         covarList    = NA)
        private$dataMatrix <- dataMatrix
        private$submodel=submodel
        },
      #--------------------------------------------
      #' @description function to perform optimization
      #' @param estimOptions a list of parameters controlling the inference algorithm and model selection. See details.
      #'
      #' @details The list of parameters \code{estimOptions} essentially tunes the optimization process and the variational EM algorithm, with the following parameters
      #'  \itemize{
      #'  \item{"nbCores"}{integer for number of cores used. Default is 2}
      #'  \item{"verbosity"}{integer for verbosity (0, 1). Default is 1}
      #'  \item{"plot"}{boolean, should the ICL by dynamically plotted or not. Default is TRUE}
      #'  \item{"exploreFactor"}{control the exploration of the number of groups}
      #'  \item{"nbBlocksRange"}{minimal and maximal number or blocks explored}
      #'  \item{"fast"}{logical: should approximation be used for Bernoulli model with covariates. Default to \code{TRUE}}
      #' }
      optimize = function(estimOptions = list()){

        currentOptions <- list(
           explorFactor  = 1.5,
           nbBlocksRange = c(1,NULL),
           nbCores       = 2,
           filename = NULL
              )
         currentOptions[names(estimOptions)] <- estimOptions

         sbmSize = list(Qmin= currentOptions$nbBlocksRange[1],
                        Qmax= if(is.na(currentOptions$nbBlocksRange[2])) NULL else currentOptions$nbBlocksRange[2],
                        explor= currentOptions$explorFactor)

        private$noisySBMobject <- fitNSBM(private$dataMatrix, model=private$submodel, sbmSize=sbmSize, filename = currentOptions$filename, nbCores= currentOptions$nbCores)

        private$import_from_noisySBM()   #### ok choisit le bon Q
       },

      graphInference = function(testLevel =0.05, qvalues=NULL){
        private$testLevel_=testLevel
         if (is.null(private$qvalues_)){
             ProcTest = noisySBM::graphInference(private$dataMatrix,nodeClustering = self$memberships, theta= private$parameters, alpha=testLevel, private$modelFamily)
              private$qvalues_=ProcTest$qvalues
              private$Y = ProcTest$A
        }
        else {private$Y = graphInferenceFromqvalues(private$qvalues_, self$nbNodes, private$testLevel_)}
        },

      #' @description method to select a specific model among the ones fitted during the optimization.
      #'  Fields of the current SBM_fit will be updated accordingly.
      #' @param index integer, the index of the model to be selected (row number in storedModels)
      setModel = function(index) {
        stopifnot(!is.null(private$noisySBMobject))
        stopifnot(index %in% seq.int(nrow(self$storedModels)))
        private$import_from_noisySBM(index)
        self$reorder()
      },
      #' @description permute group labels by order of decreasing probability
      reorder = function(){
        o <- order(private$parameters$w %*% private$pi, decreasing = TRUE)
        private$pi <- private$pi[o]
        private$parameters$w <- private$parameters$w[o,o]   #pb si on choisit Q= 1 et w n'est pas une matrice ?
        private$Z <- private$Z[, o, drop = FALSE]
      },
      # #--------------------------------------------
      #' @description show method
      #' @param type character used to specify the type of SBM
      show = function(type = "Fit of a Simple Stochastic Block Model"){
        super$show(type)
        cat("* Additional fields\n")
        cat("  $probMemberships, $loglik, $ICL, $storedModels, \n")
        cat("* Additional methods \n")
        cat("  predict, fitted, $setModel, $reorder \n")
      }
    ),

   active = list(
     #' @field qvalues
      qvalues   = function(value) {  return(private$qvalues_)},   # A VERIF
      #     return(private$qvalues_)
      # else { private$qvalues_ <- value }
      # },

     #' @field testLevel
     testLevel = function(value) {return(private$testLevel_)},    # A VERIF

      #' @field loglik double: approximation of the log-likelihood (variational lower bound) reached
      loglik = function(value) {private$J},
      #' @field ICL double: value of the integrated classification log-likelihood
      ICL    = function(value) {private$vICL},

      #' @field penalty double, value of the penalty term in ICL
      ### FANNY : Revoir : pourquoi nbConnectParam en self ????
      penalty  = function(value) {unname((self$nbConnectParam + self$nbCovariates) * log(self$nbDyads) + (self$nbBlocks-1) * log(self$nbNodes))},
      #' @field entropy double, value of the entropy due to the clustering distribution
      entropy  = function(value) {-sum(.xlogx(private$Z))},

      #' @field storedModels data.frame of all models fitted (and stored) during the optimization
      storedModels = function(value) {
        nbBlocks <- sapply(private$noisySBMobject, function(m) m$sbmParam$Q)
        N_Q <- nbConnectParam <- nbBlocks*(nbBlocks+1)/2         # A changer si dirige
        dimParamSubmodel <- nbParamSubmodel(private$submodel, nbBlocks)
        U <- data.frame(
          indexModel  = nbBlocks,
          # nbParams = nbConnectParam + nbBlocks - 1,
          nbParams = dimParamSubmodel + nbBlocks - 1,
          nbBlocks = nbBlocks,  ### FANNY : pourquoi le réafficher ????????????????????????
          ICL      = sapply(private$noisySBMobject, function(m) m$sbmParam$ICL),
          loglik   = sapply(private$noisySBMobject, function(m) m$convergence$complLogLik)
          )
        U[!is.na(U$nbParams),]
      }
    )
  )
