#' Calculates the distance scores
#'
#' @param Zx a (p x n+1) matrix the represents the amplitudes of the co-varying patterns in X
#' @param Zy a (q x n+1) matrix that represents the amplitudes of the co-varying patterns in Y
#'
#' @return a (q x n x p) array of distance scores of the elements of Zx and Zy.
#' @export
#'
#' @examples see README
get.scores <- function(Zx,Zy) {
  p <- nrow(Zx)
  q <- nrow(Zy)
  n <- ncol(Zx)-1
  score <- array(dim=c(q,n,p))
  for (i in 1:p) {
    for (j in 1:n) {
      zvec <- rep(Zx[i,j],q)
      score[,j,i] <- abs(zvec-Zy[,j])
    }
  }
  score
}
