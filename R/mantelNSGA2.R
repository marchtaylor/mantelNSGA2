#' Mantel test permutation search with NSGA2 (Non-dominated
#'   sorting genetic algorithm)
#'
#' @param fix.mat The "fixed" matrix of community or environmental sample by variable values
#' @param var.mat A "variable" matrix of community or environmental sample by variable values
#' @param fix.dist.method The method of calculating dissimilarity indices between samples in the fixed
#' matrix (Uses the \code{\link[vegan]{vegdist}} function from the vegan package to calculate distance matrices. See
#' the documentation for available methods.). Defaults to Bray-Curtis dissimilarity \code{"bray"}.
#' @param var.dist.method The method of calculating dissimilarity indices between samples in the variable
#' matrix. Defaults to Euclidean dissimilarity \code{"euclidean"}.
#' @param scale.fix Logical. Should fixed matrix be centered and scaled (Defaults to \code{FALSE},
#' recommended for biologic data).
#' @param scale.var Logical. Should fixed matrix be centered and scaled (Defaults to \code{TRUE},
#' recommended for environmental data to correct for differing units between variables).
#' @param pop.size Parent population size for the genetic algorithm.
#' @param offspring.size Offspring population size for the genetic algorithm.
#' @param max.iter Maximum number of iterations the algorithm should run.
#' @param crossover.rate Fraction of crossover between parent solutions.
#' @param mutation.rate Probability of a mutation (bitflip) within the feature vector.
#' @param allowZeroGenes Logical. If individuals containing only zeros
#'   (no features) should be included in the search (\code{TRUE}) or not (\code{FALSE}).
#' @param remove.overlap Logical. If overlapping solutions should be removed (see Nojima et al. 2005)?
#' @param p.feat Numeric. Value between 0-1 (but non-zero), indicating the maximum proportion
#'   of explanatory variables to turn on (i.e. genes) in individuals of the
#'   initial population. A lower value will assist in the search of simpler
#'   models in the earlier generations of the search.
#' @param stop.criterion Numeric. Maximum number of iterations without improvement
#'   of the solution.
#' @param memoisation Logical. If memoisation should be
#'   used (\code{TRUE}) or not (\code{FALSE}). If memoisation is used,
#'   function calls with the same input paramet
#' @param ref.point Vector of length 2, defining the Reference point for the
#'   Hypervolume Indicator.
#' @param verbose logical. Should genetic algorithm progress be printed.
#'
#' @return object of class mantelNSGA2
#' @export
#'
#' @importFrom vegan vegdist
#'
#' @examples
#'
#' library(vegan)
#' data("varechem")
#' data("varespec")
#'
#'
#' ### envbio
#' # biological community variables that best correlate with environmental data
#' fix.mat = varechem
#' var.mat = wisconsin(varespec)
#' fix.dist.method = "euclidean"
#' var.dist.method = "bray"
#' scale.fix = TRUE
#' scale.var = FALSE
#' p.feat <- 0.1
#' 2^ncol(var.mat)-1 # total combinations tested by bioEnv()
#'
#' # mantelNSGA2
#' set.seed(1111)
#' fitGA <- mantelNSGA2(
#'   fix.mat = fix.mat, var.mat = var.mat,
#'   fix.dist.method = fix.dist.method, var.dist.method = var.dist.method,
#'   scale.fix = scale.fix, scale.var = scale.var, p.feat = p.feat,
#'   pop.size = 50, max.iter = 100, stop.criterion = 50,
#'   mutation.rate = 0.2, crossover.rate = 0.8
#' )
#' plot(fitGA$generation.fitness)
#' fitGA$pareto.solution
#' plot(fitGA, parFrontT = "o")
#'
#'
#' # Compare to bvStep
#' set.seed(1111)
#' fitBV <- bvStep(
#'   fix.mat = fix.mat, var.mat = var.mat,
#'   fix.dist.method = fix.dist.method, var.dist.method = var.dist.method,
#'   scale.fix = scale.fix, scale.var = scale.var,
#'   num.restarts = 50,
#'   random.selection = TRUE,
#'   prop.selected.var = 0.4
#' )
#' fitBV$order.by.best # mantelGA not looking in simpler solutions enough
#' fitBV$order.by.i.comb
#' points(-fitBV$order.by.i.comb$rho, fitBV$order.by.i.comb$n.var/ncol(var.mat),
#'   pch = 2, col = 4, t = "o", lty = 2)
#'
#'
#' ### biobio
#' # biological community variables that best correlate with the
#' # overall biological community
#' fix.mat = wisconsin(varespec)
#' var.mat = wisconsin(varespec)
#' fix.dist.method = "bray"
#' var.dist.method = "bray"
#' scale.fix = FALSE
#' scale.var = FALSE
#' p.feat <- 0.1
#'
#' # mantelNSGA2
#' set.seed(1111)
#' fitGA <- mantelNSGA2(
#'   fix.mat = fix.mat, var.mat = var.mat,
#'   fix.dist.method = fix.dist.method, var.dist.method = var.dist.method,
#'   scale.fix = scale.fix, scale.var = scale.var, p.feat = p.feat,
#'   pop.size = 50, max.iter = 100, stop.criterion = 50,
#'   mutation.rate = 0.2, crossover.rate = 0.8
#' )
#' plot(fitGA$generation.fitness)
#' fitGA$pareto.solution
#' plot(fitGA, parFrontT = "o")
#'
#'
#' # Compare to bvStep
#' set.seed(1111)
#' fitBV <- bvStep(
#'   fix.mat = fix.mat, var.mat = var.mat,
#'   fix.dist.method = fix.dist.method, var.dist.method = var.dist.method,
#'   scale.fix = scale.fix, scale.var = scale.var,
#'   num.restarts = 50,
#'   random.selection = TRUE,
#'   prop.selected.var = 0.3, var.always.include = c(15,23,26)
#' )
#' fitBV$order.by.best # mantelGA not looking in simpler solutions enough
#' fitBV$order.by.i.comb
#' points(-fitBV$order.by.i.comb$rho, fitBV$order.by.i.comb$n.var/ncol(var.mat),
#'   pch = 2, col = 4, t = "o", lty = 2)
#'
#'
#'
#'
#'
#'
#'
mantelNSGA2 <- function(
  fix.mat,
  var.mat,
  fix.dist.method = "bray",
  var.dist.method = "euclidean",
  scale.fix = FALSE,
  scale.var = TRUE,
  pop.size = 100,
  offspring.size = NULL,
  crossover.rate = 0.8,
  mutation.rate = 0.1,
  max.iter = 150,
  allowZeroGenes = FALSE,
  remove.overlap = TRUE,
  p.feat = 1,
  stop.criterion = NULL,
  memoisation = TRUE,
  ref.point = c(0,1),
  verbose = TRUE
){

  if (is.null(offspring.size)) {
    offspring.size = pop.size
  }

  if(scale.fix){fix.mat <- scale(fix.mat)}else{fix.mat <- fix.mat}
	if(scale.var){var.mat <- scale(var.mat)}else{var.mat <- var.mat}

  fix.dist <- vegan::vegdist(as.matrix(fix.mat), method=fix.dist.method)

  PARS <- rep(0, ncol(var.mat))
  n.feat <- length(PARS)

  fn <- function(feat, var.mat, fix.dist, var.dist.method){
    vars.incl <- seq(ncol(var.mat))[as.logical(feat)]
    if(length(vars.incl) != 0){
      var.dist <- suppressWarnings(
        vegan::vegdist(as.matrix(var.mat[,vars.incl]),
          method = var.dist.method))
      temp <- suppressWarnings(cor.test(fix.dist, var.dist, method="spearman"))
      score <- temp$estimate
    } else {
      score <- 0
    }
    ret <- c(metric = -as.vector(score), frac.feat = length(vars.incl)/length(feat))
    return(ret)
  }


  if (memoisation) {
    fn = memoise::memoise(f = fn)
  }


  ctrl = ecr::initECRControl(fitness.fun = fn, n.objectives = 2, minimize = c(TRUE, TRUE))
  ctrl = ecr::registerECROperator(ctrl, "mutate", ecr::mutBitflip)
  ctrl = ecr::registerECROperator(ctrl, "recombine", ecr::recUnifCrossover)
  ctrl = ecr::registerECROperator(ctrl, "selectForMating", ecr::selSimple)
  ctrl = ecr::registerECROperator(ctrl, "selectForSurvival", ecr::selNondom)

  population = ecr::genBin(pop.size, n.feat)

  if(p.feat < 1){
    max.feat <- round(n.feat*p.feat)
    if(any(sapply(population, sum) > max.feat)){
      indx = which(sapply(population, sum) > max.feat)
      for(i in indx){
        feat.i = which(population[[i]]==1)
        off <- sample(feat.i, sum(population[[i]])-max.feat)
        population[[i]][off] <- 0
      }
    }
  }

  if (allowZeroGenes == F) {
    if (any(sapply(population, sum) == 0)) {
      indx = which(sapply(population, sum) == 0)
      for(i in indx){
        on <- sample(n.feat, 1)
        population[[i]][on] <- 1
      }
    }
  }


  fitness = ecr::evaluateFitness(inds = population, control = ctrl, var.mat = var.mat, fix.dist = fix.dist,
    var.dist.method = var.dist.method)

  logger = ecr::initLogger(control = ctrl, log.stats = list(fitness = list(HV = list(fun = ecr::computeHV,
    pars = list(ref.point = ref.point)))),
    log.pop = TRUE, init.size = max.iter + 1L)
  ecr::updateLogger(log = logger, population = population,
    fitness = fitness, n.evals = offspring.size)
  store.offspring = rep(list(NA), max.iter)
  for (i in seq_len(max.iter)) {
    offspring = ecr::recombinate(control = ctrl, inds = population,
      fitness = fitness, lambda = offspring.size, p.recomb = crossover.rate)
    offspring = ecr::mutate(control = ctrl, inds = offspring,
      p.mut = mutation.rate)

    if (!allowZeroGenes) {
      if (any(sapply(offspring, sum) == 0)) {
        indx = which(sapply(offspring, sum) == 0)
        for (zeroGenes in indx) {
          select.offspring = sample(1:length(offspring),
            1, replace = T)
          replace.offspring = offspring[[select.offspring]]
          rnd.gene = sample(1:n.feat, 1, replace = T)
          replace.offspring[rnd.gene] = 1
          offspring[[zeroGenes]] = replace.offspring
        }
      }
    }

    if (remove.overlap) {
      pop.tmp = do.call(rbind, population)
      offspr.tmp = do.call(rbind, offspring)
      offspr.tmp = unique(offspr.tmp)
      dubs = rep(NA, nrow(offspr.tmp))
      for (mm in 1:nrow(offspr.tmp)) {
        dubs[mm] = any(apply(pop.tmp, 1, function(x) identical(x,
          offspr.tmp[mm, ])))
      }
      index.rm = which(dubs == T)
      if (!pracma::isempty(index.rm)) {
        offspr.tmp = offspr.tmp[-index.rm, ]
      }
      offspring = split(offspr.tmp, seq(nrow(offspr.tmp)))
    }

    if(verbose){
      cat(paste("NSGA-II | iter =", i, "| pop.size =", length(population),
        "| offspr.size =", length(offspring)))
      cat("\n")
      utils::flush.console()
    }

    fitness.o = ecr::evaluateFitness(inds = offspring, control = ctrl,
      var.mat = var.mat, fix.dist = fix.dist, var.dist.method = var.dist.method)

    sel = ecr::replaceMuPlusLambda(control = ctrl, population = population,
      offspring = offspring, fitness = fitness, fitness.offspring = fitness.o)
    population = sel$population #
    fitness = sel$fitness # update fitness for remaining population
    ecr::updateLogger(logger, population, fitness = fitness,
      n.evals = length(offspring))
    store.offspring[[i]] = list(offspring = offspring, fitness.o = fitness.o)
    if (!is.null(stop.criterion)) {
      if (i == 1) {
        best.fitness = ecr::getStatistics(logger)$fitness.HV[1]
        n = 1
      } else {
        new.fitness = ecr::getStatistics(logger)$fitness.HV[i]
        if (best.fitness == new.fitness) {
          n = n + 1
        }
        else {
          n = 1
        }
        best.fitness = new.fitness
        if (n >= stop.criterion) {
          (break)()
        }
      }
    }
  }

  stats.fitness = ecr::getStatistics(logger)
  stats.populations = ecr::getPopulations(logger)
  Pareto.set = ecr::initParetoArchive(control = ctrl)
  ecr::updateParetoArchive(Pareto.set, inds = population,
    fitness = fitness)
  pareto.front = t(ecr::getFront(Pareto.set))

  fit.individuals = do.call(rbind, ecr::getIndividuals(Pareto.set))
  indx.uniq = which(duplicated(pareto.front) == F)
  pareto.front.uniq = data.frame(ID = indx.uniq, pareto.front[indx.uniq,])
  pareto.individuals = fit.individuals[indx.uniq, ]
  if (is.null(dim(pareto.individuals))) {
    pareto.front.uniq$nr.of.params = sum(pareto.individuals)
    pareto.vars = list(seq(ncol(var.mat))[which(pareto.individuals ==
      1)])
  }else{
    pareto.front.uniq$nr.of.params = apply(pareto.individuals,
      1, sum)
    pareto.vars = list()
    for (j in 1:nrow(pareto.individuals)) {
      pareto.vars[[j]] = seq(ncol(var.mat))[which(pareto.individuals[j,
        ] == 1)]
    }
  }
  pareto.front.uniq$dist.to.CoordOrigin = apply(pareto.front.uniq[,
    c(2, 3)], 1, function(x) sqrt(sum((x - c(-1, 0))^2)))

  sol <- pareto.front.uniq
  sol$vars <- unlist(lapply(pareto.vars, function(x){paste0(x, collapse = ",")}))
  sol <- sol[order(sol$nr.of.params),]


  ret <- list(generation.fitness = stats.fitness, generation.populations = stats.populations,
    generation.offspring = store.offspring, pareto.individuals = pareto.individuals,
    pareto.vars = pareto.vars, pareto.solution = sol)

  if(memoisation) memoise::forget(fn)

  class(ret) <- "mantelNSGA2"

  return(ret)

}
