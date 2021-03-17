source("R/R6Class-NoisySBM_fit.R")
source("R/R6Class-NoisySBM.R")
source("R/R6Class-SBM.R")
source("R/utils_plot.R")
source("R/utils.R")
source("R/plotMyMatrix.R")
source("R/newFcts.R")

library(ggplot2)
library(purrr)
#library(noisySBM)
#simu=rnsbm(n=10,theta=list(pi=rep(1/2,2), w=c(0.8,0.1,0.9), nu0=c(0,1), nu=matrix(c(10,10,10,.2,.2,.2), ncol=2)))
#dataMatrix=simu$dataMatrix

 nbNodes  <- 20
 blockProp <- c(.5, .5)      # group proportions
 connectParam=list(mean=matrix(c(0.8,0.1,0.1,0.9), ncol=2, nrow=2))
 noiseParam=c(0,1)
 signalParam=matrix(c(10,10,10,.2,.2,.2), ncol=2)
 ## Graph Sampling

 nbNodes  <- 10
 blockProp <- c(.5, .5)      # group proportions
 connectParam=list(mean=matrix(c(0.5,0.9,0.9,0.3), ncol=2, nrow=2))
 noiseParam=c(0,.1)
 signalParam=matrix(c(-2,10,-2,1,1,1), ncol=2)
 signalParam=matrix(c(0,0,1,1,1,1), ncol=2)
 ## Graph Sampling


 mySampler <- sampleNoisySBM(model = 'Gauss', nbNodes, directed = FALSE, blockProp, connectParam, noiseParam, signalParam)
 plot(mySampler)
 plot(mySampler,type='meso')
#plot(mySampler,type='latentNetwork')     # ne fonctionne pas : changer plot. rajouter fonction dans utils_plot.R ?
 hist(mySampler$networkObs)   # X
 hist(mySampler$networkData)  # latentNetwork
 plotMyMatrix(mySampler$networkObs)
 plotMyMatrix(mySampler$networkData)

 dataMatrix= mySampler$networkObs
nous= fitNSBM(dataMatrix,'Gauss01', sbmSize=list(Qmax=2))


 library(noisySBM)
 tata = NoisySBM_fit$new(dataMatrix, submodel="Gauss01")
 #tata$optimize()   #param retire de active
 tata$optimize(list(nbBlocksRange=c(1,2)))
 tata$blockProp
 tata$storedModels
tata$connectParam

tata$graphInference()
tata$testLevel
tata$qvalues
tata$networkData

tata$graphInference(0.5)[1:10,1:10]
tata$testLevel
tata$qvalues

tata$graphInference(0.8, rep(0, 100*99/2))[1:10,1:10]
tata$testLevel
tata$qvalues




######## VIEUX
toto = NoisySBM$new(model="Gauss", nbNodes=100, directed=F, blockProp = rep(1/2,2), connectParam=list(mean=matrix(c(0.8,0.1,0.1,0.9), ncol=2, nrow=2)), noiseParam=c(0,1), signalParam=matrix(c(10,10,10,.2,.2,.2), ncol=2))
toto$connectParam
toto$rMemberships(TRUE)
toto$rEdges()
toto$rNetwork()
toto$rEdges(TRUE)
toto$networkData   #latentAdj : même nom que eux mais foireux chez nous
toto$networkObs    #dataMatrix
toto$plot(type="data")
toto$plot(type="data", ordered=FALSE)
toto$plot("latentNetwork")
toto$plot("expected")
toto$plot("meso")
#parametre = list(pi = toto$blockProp, w=toto$connectParam, nu0=toto$noiseParam, nu=toto$signalParam)
#Z <- t(rmultinom(n=10, size = 1, prob = toto$blockProp))
#rnsbmObs(parametre, Z, modelFamily='Gauss', directed=FALSE)

dataMatrix= toto$networkObs

library(noisySBM)
#simu=rnsbm(n=10,theta=list(pi=rep(1/2,2), w=c(0.8,0.1,0.9), nu0=c(0,1), nu=matrix(c(10,10,10,.2,.2,.2), ncol=2)))
#dataMatrix=simu$dataMatrix

tata = NoisySBM_fit$new(dataMatrix, submodel="Gauss01")
tata$optimize()   #param retire de active
tata$optimize(list(nbBlocksRange=c(1,2)))
tata$plot("data")
tata$plot("latentNetwork")    # ne devrait pas fonctionner tant qu'on n'a pas codé tests multiples
tata$plot("meso")
tata$plot("expected")
tata$blockProp
tata$storedModels
tata$ICL
tata$reorder()
print(tata$reorder())

tata$setModel(2)
tata$ICL   # attention, ca prend le modele 2 !!!
tata$plot("meso")
tata$plot("expected")
tata$blockProp
tata$param
print(tata$reorder())

tata$graphInference()
tata$graphInference(0.5)

#coef(tata)        #REVOIR

# plot(tata$networkData, type="data")
# toto$plot()
# toto$plot(type="data")
# plot(toto)
# tata$plot()
# plotMyMatrix(simu$latentAdj)
# plotMatrix(simu$latentAdj)
# plotMatrix(simu$latentAdj, dimLabels = c(rep("a",2)))
# plotMatrix(simu$dataMatrix, dimLabels = c(rep("a",2)))


library(purrr)
titi = SBM$new(
  model = "bernoulli",
  dimension = c(10,10),
  dimLabels = rep("nodes",2),
  blockProp = rep(1/2,2),
  connectParam = list(mean = matrix(c(0.8,0.1,0.8,0.1), ncol=2)),
  # covarParam = numeric(length(covarList)),
  covarList = list()
)


source("R/R6Class-BipartiteSBM.R")
titi = BipartiteSBM$new(
   model = "bernoulli",
   nbNodes =c(10,10),
   blockProp = list(rep(1/2,2),rep(1/2,2)),
   connectParam = list(mean = matrix(c(0.8,0.1,0.8,0.1), ncol=2)),
   covarParam = numeric(length(covarList)),
   covarList = list()
)



source("R/R6Class-MultipartiteSBM.R")
titi=MultipartiteSBM$new(
  model = "bernoulli",
  architecture = matrix(0, 0, 2),
  directed = FALSE,
  dimension = c(10,10),
  dimLabels = character(0),
  blockProp = list(),
  connectParam = list()
)
