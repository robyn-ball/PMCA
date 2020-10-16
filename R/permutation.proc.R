#' Permutation procedure
#'
#' @param X a (p x n) matrix
#' @param Y a (q by n) matrix
#' @param method a character string, "overall" or "each" denoting the FPR method
#' @param B an integer denoting the number of permutations
#'
#' @return a 4 dimensional array of distance scores for each of the B permutations
#' @export
#'
#' @examples see README
permutation.proc <- function(X,Y,method="overall",B=1000) {
  p <- nrow(X); q <- nrow(Y); N <- ncol(X);
  J <- min(p,q)-1
  scores <- array(dim=c(q,J,p,B))
  shuffle <- function(x) {
    x[sample(length(x))]
  }
  if (method=="overall") {
    for (b in 1:B) {
      Xstar <- apply(X,2,shuffle)
      mca <- get.mca(Xstar,Y)
      scores[,,,b] <- get.scores(mca$Zx,mca$Zy)
    }
  } else { #method=FPR for each rowterm of X
    for (b in 1:B) {
      Ystar <- t(apply(Y,1,shuffle))
      mca <- get.mca(X,Ystar)
      scores[,,,b] <- get.scores(mca$Zx,mca$Zy)
    }
  }
  scores
}
