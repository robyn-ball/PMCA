#' Produce the lowest pairwise FPR associations and anti-associations
#'
#' @param g a list that maps rowterms of Y onto rowterms of X for each component
#' @param fpr a matrix (p x p-1) of FPRs where the rows represent the rowterms of X and
#'  the columns represent the p-1 components of the covariance matrix related to X
#' @param X a matrix (p x n) where the columns represent the coupled data (i.e., sample)
#'
#' @return a dataframe of FPRs for each pairwise relationship contained in g
#'
#' @export
#'
#' @examples see README
get.mapped.fpr <- function(g, fpr, X) {
  # number of components (p-1)
  J <- length(g[[1]])
  mapped <- NULL
  for (j in J:2) {
    # for the jth component
    for (i in 1:nrow(X)) {
      # for each rowterm i in X, cell contains the rowterms of Y that mapped to the i-th rowterm of X with a FPR <= fpr[i,j]
      cell <- g[[i]][[j]]
      if (length(cell)>0) {
        for (c in cell) {
          # for each of these rowterms of Y, add it as a column of mapped if it does not exist
          # enter the FPR for this relationship
          if (!is.element(c,colnames(mapped))) {
            mapped <- cbind(mapped,rep(NA,nrow(X)))
            colnames(mapped)[ncol(mapped)] <- c
            mapped[i,ncol(mapped)] <- fpr[i,j]
          } else {
            # else, the column for this rowterm of Y exists. Enter the FPR for this relationship
            is.cell <- which(colnames(mapped)==c)
            if (is.na(mapped[i,is.cell])) {
              mapped[i,is.cell] <- fpr[i,j]
            }
          }
        }
      }
    }
  }
  mapped
}
