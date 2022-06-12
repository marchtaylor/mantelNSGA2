#' plot results from mantelNSGA2
#'
#' @param obj bla
#' @param xlim  bla
#' @param ylim  bla
#' @param genCol  bla
#' @param genPch  bla
#' @param genCex bla
#' @param parFrontCol  bla
#' @param parFrontBg  bla
#' @param parFrontPch  bla
#' @param parFrontCex  bla
#' @param parFrontT  bla
#'
#' @return a plot
#' @method plot mantelNSGA2
#'
#' @export
#'
#' @importFrom graphics points
#' @importFrom grDevices hcl.colors
#'
#' @examples
#'
#' library(vegan)
#' data("varechem")
#' data("varespec")
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
plot.mantelNSGA2 <- function(obj, xlim = c(-1,0), ylim = c(0,1),
  genCol = hcl.colors(12, "viridis", alpha = 0.3), genPch = ".", genCex = 3,
  parFrontCol = 1, parFrontBg = 2, parFrontPch = 21, parFrontCex = 1, parFrontT = "p"
    ){

  plot(obj$pareto.solution[,c("metric", "frac.feat")], xlim = xlim, ylim = ylim, t = "n")
  ngen <- length(obj$generation.populations)
  COL <- val2col(z = seq(ngen), col = genCol)
  for(i in seq(ngen)){
    points(t(obj$generation.populations[[i]]$fitness), pch = genPch, cex = genCex, col = COL[i])
  }
  ord <- obj$pareto.solution[,c("metric", "frac.feat")]
  ord <- ord[order(ord$frac.feat),]
  points(ord[,c("metric", "frac.feat")], pch = parFrontPch, col = parFrontCol, bg = parFrontBg, cex = parFrontCex,
    t = parFrontT)

}


