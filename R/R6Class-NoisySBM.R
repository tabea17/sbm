#' R6 class for Simple SBM
#'
#'
#' @export
NoisySBM <-
  R6::R6Class(
    classname = "NoisySBM",
    inherit = SBM,

    private=list(
      modelFamily=NULL,
      parameters = NULL,
      dataMatrix = NULL #,
     ),

    public = list(
      #' @description constructor for noisySBM
      #' @param modelFamily character describing the type of model  (gauss ou exp ou poisson)
      #' @param nbNodes number of nodes in the network
      #' @param directed logical, directed network or not. Code uniquement en non dirige pour le moment
      #' @param blockProp parameters for block proportions (vector of list of vectors)
      #' @param connectParam list of parameters for connectivity with a matrix of means 'mean' and an optional scalar for the variance 'var'. The size of mu must match \code{blockProp} length
      #' @param covarParam optional vector of covariates effect
      #' @param covarList optional list of covariates data
      #' @param noiseParam list of parameters for connectivity ...
      #' @param signalParam list of parameters for connectivity ...
      #' @param dimLabels optional label for the node
      initialize = function(modelFamily, nbNodes, directed=FALSE, blockProp, connectParam, noiseParam, signalParam, dimLabels=c("nodeName"), covarParam=numeric(length(covarList)), covarList=list()) {

        ## SANITY CHECKS (on parameters)
       stopifnot(length(dimLabels) == 1)
       stopifnot(is.atomic(blockProp), all(blockProp > 0), all(blockProp < 1)) # positive proportions
       stopifnot(all.equal(length(blockProp), ncol(connectParam$mean)),        # dimensions match between vector of
       all.equal(length(blockProp), nrow(connectParam$mean)))        # block proportion and connectParam$mean
       # if (!directed) stopifnot(isSymmetric(connectParam$mean)) # connectivity and direction must agree

        private$modelFamily  = modelFamily
        private$parameters = list(pi = blockProp, w=connectParam$mean, nu0=noiseParam, nu=signalParam)
        super$initialize(model='bernoulli', directed, nbNodes, dimLabels, blockProp, connectParam, covarParam, covarList)
         },

     #' #' @description a method to sample new block memberships for the current SBM
     #' #' @param store should the sampled blocks be stored (and overwrite the existing data)? Default to FALSE
     #' #' @return the sampled blocks
      rMemberships = function(store = FALSE) {
        Z <- t(rmultinom(private$dim, size = 1, prob = private$pi))
        if (store) private$Z <- Z
        Z
      },
     #'  #' @description a method to sample a network data (edges) for the current SBM
     #'  #' @param store should the sampled edges be stored (and overwrite the existing data)? Default to FALSE
     #'  #' @return the sampled network
      rEdges = function(store = FALSE) {
        YandX <- rnsbmObs(private$parameters, private$Z, private$modelFamily, private$directed_)
        if (store) private$Y <- YandX$latentAdj
        if (store) private$dataMatrix <- YandX$dataMatrix
        YandX   # attention liste chez nous
      },

      #'  #--------------------------------------------
     #'  #' @description prediction under the currently parameters
     #'  #' @param covarList a list of covariates. By default, we use the covariates with which the model was estimated
     #'  #' @param theta_p0 a threshold...
     #'  #' @return a matrix of expected values for each dyad
      predict = function(covarList = NULL, theta_p0 = NULL){          # donne grande matrice des w_ql. On a laisse covarList et theta_p0 pour que ca fonctionne avec predict.SBM mais aucun sens pour nous.
        mu <- private$Z %*% private$parameters$w %*% t(private$Z)
        mu
      },

      #' @description show method
      #' @param type character used to specify the type of SBM
      show = function(type = "Noisy Stochastic Block Model") {super$show(type)},

      #' @description basic matrix plot method for SimpleSBM object or mesoscopic plot
      #' @param type character for the type of plot: either 'data' (true connection), 'expected' (fitted connection) or 'meso' (mesoscopic view). Default to 'data'.
      #' @param ordered logical: should the rows and columns be reordered according to the clustering? Default to \code{TRUE}.
      #' @param plotOptions list with the parameters for the plot. See help of the corresponding S3 method for details.
      #' @return a ggplot2 object for the \code{'data'} and \code{'expected'}, a list with the igraph object \code{g}, the \code{layout} and the \code{plotOptions} for the \code{'meso'}
      #' @import ggplot2
      plot = function(type = c('data','expected','meso', 'latentNetwork'), ordered = TRUE, plotOptions = list()) {

        if(is.null(self$memberships)) {ordered <- FALSE; type <- 'data'}
        if (ordered & !is.null(self$memberships))
          clustering <- list(row = self$memberships)
        else
          clustering <- NULL

        switch(match.arg(type),
          "meso" =
            plotMeso(
              thetaMean  = private$parameters$w, #private$theta$mean,
              pi         = private$pi,
              model      = private$model,
              directed   = private$directed_,
              bipartite  = FALSE,
              nbNodes    = self$nbNodes,
              nodeLabels = as.list(private$dimlab),
              plotOptions),
          "data" =
            plotMatrix(self$networkObs,
                       private$dimlab,
                       clustering, plotOptions),

          "latentNetwork" =
            plotMatrix(self$networkData,
                       private$dimlab,
                       clustering, plotOptions),
          "expected" =
            plotMatrix(self$expectation,
                       private$dimlab,
                       clustering, plotOptions)
        )
      }
    ),

    active = list(
    ### field with write access
      #' @field dimLabels a single character giving the label of the nodes
       dimLabels    = function(value) {
         if (missing(value))
           return(private$dimlab)
         else {
           stofifnot(is.atomic(value), is.character(value), length(value) == 1)
           if (is.null(names(value))){names(value)  = c('node')}
           private$dimlab <- value
         }
       },
      #' @field blockProp vector of block proportions (aka prior probabilities of each block)
      blockProp   = function(value) {
        if (missing(value))
          return(private$pi)
        else {
     #     stopifnot(is.numeric(value), is.atomic(value),
     #               all(value > 0), all(value < 1)) # positive proportions
          private$pi <- value
          private$parameters$pi <- value
        }
      },
      #' @field connectParam parameters associated to the connectivity of the SBM, e.g. matrix of inter/inter block probabilities when model is Bernoulli
      connectParam   = function(value) {
        if (missing(value))
        return(list(mean=private$parameters$w))
        else {
           stopifnot(is.list(value))
           # if (!self$directed) stopifnot(isSymmetric(value$mean)) # connectivity and direction must agree
          private$parameters$w <- value
        }
      },
      #' @field signalParam parameters
      signalParam   = function(value) {
        if (missing(value))
          return(private$parameters$nu)
        else { private$parameters$nu <- value }
      },
      #' @field noiseParam parameters
      noiseParam   = function(value) {
        if (missing(value))
          return(private$parameters$nu0)
        else { private$parameters$nu0 <- value }
      },
      #' @field probMemberships  matrix of estimated probabilities for block memberships for all nodes
      probMemberships = function(value) {
        if (missing(value))
          return(private$Z)
        else {
          stopifnot(nrow(value) == private$dim)
          private$Z <- value
        }
      },


     ### field with access only
      #' @field nbBlocks number of blocks
      nbBlocks    = function(value) {length(private$pi)},
      #' @field nbDyads number of dyads (potential edges in the network)
      nbDyads     = function(value) {ifelse(private$directed_, self$nbNodes*(private$dim - 1), private$dim*(private$dim - 1)/2)},
      #' @field nbConnectParam number of parameter used for the connectivity
      nbConnectParam = function(value) {ifelse(private$directed_, self$nbBlocks^2, self$nbBlocks*(self$nbBlocks + 1)/2)},
      #' @field memberships vector of clustering
      memberships = function(value) {if (!is.null(private$Z)) as_clustering(private$Z)},
      #' @field indMemberships matrix for clustering memberships
      indMemberships = function(value) {as_indicator(as_clustering(private$Z))},
      networkObs = function(value) {return(private$dataMatrix)} #,   ##FANNY : j'ai rajoute
    )
  )


