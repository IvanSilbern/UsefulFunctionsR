# from https://stats.stackexchange.com/questions/259208/how-to-perform-isometric-log-ratio-transformation/270203#270203

ilr <- function(x, p=0) {
  y <- log(x)
  if (p != 0) y <- (exp(p * y) - 1) / p       # Box-Cox transformation
  y <- y - rowMeans(y, na.rm=TRUE)            # Recentered values
  k <- dim(y)[2]
  H <- contr.helmert(k)                       # Dimensions k by k-1
  H <- t(H) / sqrt((2:k)*(2:k-1))             # Dimensions k-1 by k
  if(!is.null(colnames(x)))                   # (Helps with interpreting output)
    colnames(y) <- paste0(colnames(x)[-1], ".ILR")
  return(y %*% t(H))                          # Rotated/reflected values
}
