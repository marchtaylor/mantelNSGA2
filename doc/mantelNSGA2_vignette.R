## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  echo = TRUE, 
  collapse = TRUE,
  comment = "#>",
  fig.width = 6, fig.height = 5
)

## ----setup, warning=FALSE, message=FALSE--------------------------------------
library(mantelNSGA2)
library(vegan)
data("varechem")
data("varespec")

## -----------------------------------------------------------------------------
# chance of having few or many "genes" turned on is low, so one ends up
# mainly exploring the middle
nvar <- 20
tmp <- replicate(1000, expr = sum(sample(x = c(0,1), nvar, replace = TRUE)))
hist(tmp, breaks = seq(0, nvar))

# more extreme when number of genes is high
nvar <- 100
tmp <- replicate(1000, expr = sum(sample(x = c(0,1), nvar, replace = TRUE)))
hist(tmp, breaks = seq(0, nvar))


## ----envbio_nsga2, message=FALSE----------------------------------------------
# biological community variables that best correlate with environmental data
fix.mat = varechem
var.mat = wisconsin(varespec)
fix.dist.method = "euclidean"
var.dist.method = "bray"
scale.fix = TRUE
scale.var = FALSE
p.feat <- 0.1

# mantelNSGA2
set.seed(1111)
fitGA <- mantelNSGA2(
  fix.mat = fix.mat, var.mat = var.mat,
  fix.dist.method = fix.dist.method, var.dist.method = var.dist.method,
  scale.fix = scale.fix, scale.var = scale.var, p.feat = p.feat,
  pop.size = 50, max.iter = 100, stop.criterion = 50,
  mutation.rate = 0.2, crossover.rate = 0.8, verbose = FALSE
)
# plot(fitGA$generation.fitness)
# fitGA$pareto.solution
plot(fitGA, parFrontT = "o")


## ----envbio_nsga2_vs_bvstep, warning=FALSE------------------------------------
# Compare to bvStep
set.seed(1111)
fitBV <- bvStep(
  fix.mat = fix.mat, var.mat = var.mat,
  fix.dist.method = fix.dist.method, var.dist.method = var.dist.method,
  scale.fix = scale.fix, scale.var = scale.var,
  num.restarts = 50, 
  random.selection = TRUE,
  prop.selected.var = 0.4, 
  verbose = FALSE
)
# fitBV$order.by.best # mantelGA not looking in simpler solutions enough
# fitBV$order.by.i.comb
plot(fitGA, parFrontT = "o")
points(-fitBV$order.by.i.comb$rho, fitBV$order.by.i.comb$n.var/ncol(var.mat), pch = 1, col = 4, t = "o")

## ----biobio_nsga2, message=FALSE----------------------------------------------
# biological community variables that best correlate with the
# overall biological community
fix.mat = wisconsin(varespec)
var.mat = wisconsin(varespec)
fix.dist.method = "bray"
var.dist.method = "bray"
scale.fix = FALSE
scale.var = FALSE
p.feat <- 0.1

# mantelNSGA2
set.seed(1111)
fitGA2 <- mantelNSGA2(
  fix.mat = fix.mat, var.mat = var.mat,
  fix.dist.method = fix.dist.method, var.dist.method = var.dist.method,
  scale.fix = scale.fix, scale.var = scale.var, p.feat = p.feat,
  pop.size = 50, max.iter = 100, stop.criterion = 50,
  mutation.rate = 0.2, crossover.rate = 0.8, verbose = FALSE
)
# plot(fitGA2$generation.fitness)
# fitGA2$pareto.solution
plot(fitGA2, parFrontT = "o")


## ----biobio_nsga2_vs_bvstep, warning=FALSE------------------------------------
# Compare to bvStep
set.seed(1111)
fitBV2 <- bvStep(
  fix.mat = fix.mat, var.mat = var.mat,
  fix.dist.method = fix.dist.method, var.dist.method = var.dist.method,
  scale.fix = scale.fix, scale.var = scale.var,
  num.restarts = 50, 
  random.selection = TRUE,
  prop.selected.var = 0.3, 
  var.always.include = c(15,23,26), 
  verbose = FALSE
)
# fitBV2$order.by.best # mantelGA not looking in simpler solutions enough
# fitBV2$order.by.i.comb
plot(fitGA2, parFrontT = "o")
points(-fitBV2$order.by.i.comb$rho, fitBV2$order.by.i.comb$n.var/ncol(var.mat), pch = 1, col = 4, t = "o")

