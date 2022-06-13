## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  cache = TRUE, # TRUE for speed. Remember to empty cache after changes to cached R code sections
  cache.path = "cache/", fig.path = "tex/",
  echo = TRUE, 
  collapse = TRUE,
  comment = "#>",
  fig.width = 6, fig.height = 5, dpi=300, 
  out.width="600px", out.height="500px"
)

## -----------------------------------------------------------------------------
trapzfun <- function(x, y){
  idx = 2:length(x)
  return(as.double((x[idx] - x[idx - 1]) %*% (y[idx] + y[idx - 
    1]))/2)
} 
nvar <- c(8, 10, 12, 14, 16)
RES <- vector("list", length(nvar))
names(RES) <- nvar
for(j in seq(nvar)){
  res <- seq(nvar[j])*NaN
  for(i in seq(res)){
    res[i] <- ncol(combn(seq(length(res)), i))
  }
  RES[[j]] <- res/trapzfun(x = seq(res)/length(res), y = res)
}

ylim = c(0, max(do.call("c", RES)))
for(j in seq(RES)){
  if(j == 1){
    plot(seq(RES[[j]])/length(RES[[j]]), RES[[j]], ylim = ylim, xlim = c(0,1),
      xlab = "Parameters used [%]", ylab = "Combinations [density]", t = "n")
  }
  lines(seq(RES[[j]])/length(RES[[j]]), RES[[j]], col = j)
}
legend("topright", legend = rev(nvar), col = rev(seq(nvar)), lty = 1, title = "Variables [n]", bty = "n")

