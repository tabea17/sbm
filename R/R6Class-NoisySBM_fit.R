#' R6 Class definition of a Simple SBM fit
#'
#' This class is designed to give a representation and adjust an SBM fitted with blockmodels.
#'
#' @import R6 blockmodels
#' @export
NoisySBM_fit <-
  R6::R6Class(
    classname = "NoisySBM_fit",
    inherit = NoisySBM,
    private = list(
      J              = NULL, # approximation of the log-likelihood
      vICL           = NULL, # approximation of the ICL
      noisySBMobject       = NULL, # blockmodels output (used to stored the optimization results when blockmodels is used)
      import_from_noisySBM = function(Q = noisySBM::getBestQ(private$noisySBMobject)$Q) { # a function updating the Class
        private$J     <- private$noisySBMobject[[Q]]$convergence$J
        private$vICL  <- private$noisySBMobject[[Q]]$sbmParam$ICL
        parameters    <- private$noisySBMobject[[Q]]$theta
        private$w <- parameters$w
        private$nu0 <- parameters$nu0
        private$nu <- parameters$nu
      #  private$beta  <- parameters$beta ## NULL if no covariates
     #   private$theta <- switch(private$BMobject$model_name,
        #   "bernoulli"                 = list(mean = parameters$pi),
        #   "bernoulli_covariates"      = list(mean = .logistic(parameters$m)),
        #   "bernoulli_covariates_fast" = list(mean = .logistic(parameters$m)),
        #   "poisson"                   = list(mean = parameters$lambda),
        #   "poisson_covariates"        = list(mean = parameters$lambda),
        #   "gaussian"                  = list(mean = parameters$mu, var = parameters$sigma2),
        #   "gaussian_covariates"       = list(mean = parameters$mu, var = parameters$sigma2),
        #   "ZIgaussian"                = list(mean = parameters$mu, var = parameters$sigma2, p0 = parameters$p0),
        # )
        private$Z  <- t(private$noisySBMobject[[Q]]$sbmParam$clusterProba)  #verif n*Q chez eux ?
        private$pi <- parameters$pi
        private$edgeProba <- private$noisySBMobject[[Q]]$sbmParam$edgeProba  #N_Q x N matrix
      }
    ),
    public = list(
      #' @description constructor for a Noisy SBM fit
      #' @param dataMatrix square (noisy) matrix
      #' @param submodel character (\code{'Gauss'}, \code{'Gauss01'})
    #  #' @param directed logical, directed network or not. In not, \code{dataMatrix} must be symmetric.
        initialize = function(dataMatrix, submodel, dimLabels= c(node='nodes')) {  #rajouter directed plus tard

        ## SANITY CHECKS (on data)
        stopifnot(is.matrix(dataMatrix))                   # must be a matrix
        stopifnot(all.equal(nrow(dataMatrix),
                            ncol(dataMatrix)))             # matrix must be square
        stopifnot(isSymmetric(dataMatrix))    # symmetry and direction must agree

        ## INITIALIZE THE SBM OBJECT ACCORDING TO THE DATA
        connectParam <- NA
        signalParam <- NA
        noiseParam <- NA

        # connectParam <- switch(model,
        #   "bernoulli"  = list(mean = matrix(0, 0, 0)),
        #   "poisson"    = list(mean = matrix(0, 0, 0)),
        #   "gaussian"   = list(mean = matrix(0, 0, 0), var = 1),
        #   "ZIgaussian" = list(mean = matrix(0, 0, 0), var = 1, p0 = 0),
        # )
    #    initialize = function(model, nbNodes, directed, blockProp, connectParam, dimLabels=c(node="nodeName"), covarParam=numeric(length(covarList)), noiseParam, signalParam) {
          if (submodel %in% c('Gauss','Gauss0','Gauss01','GaussEqVar','Gauss0EqVar','Gauss0Var1','Gauss2distr','GaussAffil'))
            model <- 'Gauss'
          if (submodel %in% c('Exp','ExpGamma','Gamma'))
            model <- 'Gamma'

        super$initialize(model        = model,
                         directed     = FALSE,
                         nbNodes      = nrow(dataMatrix),
                         blockProp    = vector("numeric", 0),
                         dimLabels    = dimLabels,
                         connectParam = connectParam,
                         signalParam = signalParam,
                         noiseParam = noiseParam)
        private$Y <- dataMatrix
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

        if(private$model == 'ZIgaussian') stop("Inference not yet implemented for ZI gaussian network")

        currentOptions <- list(
          verbosity     = 1,
          plot          = TRUE,
          explorFactor  = 1.5,
          nbBlocksRange = c(1,NULL),
          nbCores       = 2,
          fast          = TRUE
        )
        currentOptions[names(estimOptions)] <- estimOptions

        ## Transform estimOptions to a suited for blockmodels list of options
        blockmodelsOptions <- list(
          verbosity          = currentOptions$verbosity,
          plotting           = if(currentOptions$plot) character(0) else "",
          explore_min        = currentOptions$nbBlocksRange[1],
          explore_max        = currentOptions$nbBlocksRange[2],
          ncores             = currentOptions$nbCores,
          exploration_factor = currentOptions$explorFactor
         )
        fast <- currentOptions$fast

        ## generating arguments for blockmodels call
        args <- list(membership_type =  ifelse(!private$directed_, "SBM_sym", "SBM"), adj = .na2zero(private$Y))
        if (self$nbCovariates > 0) args$covariates <- private$X
        args <- c(args, blockmodelsOptions)

        ## model construction
        model_type <- ifelse(self$nbCovariates > 0, paste0(private$model,"_covariates"), private$model)
        if (model_type == 'bernoulli_covariates' & fast == TRUE) model_type <- 'bernoulli_covariates_fast'
        private$BMobject <- do.call(paste0("BM_", model_type), args)

        ## performing estimation
        private$BMobject$estimate()

        ## Exporting blockmodels output to simpleSBM_fit fields
        private$import_from_BM()

        invisible(private$BMobject)
      },
      #' @description method to select a specific model among the ones fitted during the optimization.
      #'  Fields of the current SBM_fit will be updated accordingly.
      #' @param index integer, the index of the model to be selected (row number in storedModels)
      setModel = function(index) {
        stopifnot(!is.null(private$BMobject))
        stopifnot(index %in% seq.int(nrow(self$storedModels)))
        private$import_from_BM(index)
        self$reorder()
      },
      #' @description permute group labels by order of decreasing probability
      reorder = function(){
        o <- order(private$theta$mean %*% private$pi, decreasing = TRUE)
        private$pi <- private$pi[o]
        private$theta$mean <- private$theta$mean[o,o]
        private$Z <- private$Z[, o, drop = FALSE]
      },
      #--------------------------------------------
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
      #' @field loglik double: approximation of the log-likelihood (variational lower bound) reached
      loglik = function(value) {private$J},
      #' @field ICL double: value of the integrated classification log-likelihood
      ICL    = function(value) {private$vICL},
      #' @field penalty double, value of the penalty term in ICL
      penalty  = function(value) {unname((self$nbConnectParam + self$nbCovariates) * log(self$nbDyads) + (self$nbBlocks-1) * log(self$nbNodes))},
      #' @field entropy double, value of the entropy due to the clustering distribution
      entropy  = function(value) {-sum(.xlogx(private$Z))},
      #' @field storedModels data.frame of all models fitted (and stored) during the optimization
      storedModels = function(value) {
        nbBlocks <- unlist(sapply(private$BMobject$memberships, function(m) ncol(m$Z)))
        nbConnectParam <- unlist(sapply(private$BMobject$model_parameters, function(param) param$n_parameters))
        U <- data.frame(
          indexModel  = nbBlocks,
          nbParams = nbConnectParam + nbBlocks - 1,
          nbBlocks = nbBlocks,
          ICL      = private$BMobject$ICL,
          loglik   = private$BMobject$PL
          )
        U[!is.na(U$nbParams),]
      }
    )
  )
