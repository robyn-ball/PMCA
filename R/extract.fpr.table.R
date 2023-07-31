#' Extract dataframes that contain the pairwise FPR associations and anti-associations
#'
#' @param X a matrix (p x n) where the columns represent the coupled data (i.e., sample)
#' @param Y a matrix (q x n) where the columns represent the coupled data (i.e., sample)
#' @param it.result a list of the FPR (array of FPRs), wopt (optimal window width), and tau (end scaling factor)
#'
#' @return a list of two dataframes, assoc.fpr.table FPR for each PMCA association and
#'  antiassoc.fpr.table FPR for each PMCA anti-association.
#'  row names of each dataframe are  rownames(X) and column names are the rowterms of Y that mapped onto each rowterm of X
#' @export
#'
#' @examples see README
extract.fpr.table <- function(X, Y, it.result) {
  options(stringsAsFactors = FALSE)
  source("R/get.mca.R")
  source("R/match.patterns.R")
  source("R/get.mapped.fpr.R")
  #
  # calculate the MCA of X, Y
  mca.real <- get.mca(X,Y)
  # extract the associated scores, distances between each rowterm of X and each rowterm of Y
  scores.assoc <- get.scores(mca.real$Zx, mca.real$Zy)
  # extract the anti-associated scores
  scores.antiassoc <- get.scores(-mca.real$Zx, mca.real$Zy)
  # list of associations for each rowterm of X and their FPRs
  assoc.list <- match.patterns(scores.assoc, w = it.result$wopt, rowterms = rownames(Y))
  # list of anti-associations for each rowterm of X and their FPRs
  antiassoc.list <- match.patterns(scores.antiassoc, w = it.result$wopt, rowterms = rownames(Y))
  #
  # FPR table from the permutations
  fpr <- it.result$FPR
  rownames(fpr) <- rownames(X)
  colnames(fpr) <- paste0("component=", 1:ncol(fpr))
  #
  # dataframe of associated FPRs with rowterms of X as rows, columns as rowterms of Y that mapped onto X
  assoc.fpr.table <- get.mapped.fpr(assoc.list, fpr, X)
  # dataframe of anti-associated FPRs with rowterms of X as rows, columns as rowterms of Y that mapped onto X
  antiassoc.fpr.table <- get.mapped.fpr(antiassoc.list, fpr, X)
  rownames(assoc.fpr.table) <- rownames(antiassoc.fpr.table) <- rownames(X)
  list(assoc.fpr.table=assoc.fpr.table, antiassoc.fpr.table=antiassoc.fpr.table)
}

