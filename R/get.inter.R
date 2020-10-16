#' Counts of mappings (and their intersections) for a given list, g, and component, J.
#'
#' @param g a nested list
#' @param J component
#'
#' @return an array of dimensions length(g) x length(g)
#' @export
#'
#' @examples see README
get.inter <- function(g,J=by) {
  inter <- array(dim=c(length(g),length(g)))
  for (i in 1:length(g)) {
    for (j in 1:length(g)) {
      inter[i,j] <- length(intersect(g[[i]][[J]],g[[j]][[J]]))
    }
  }
  inter
}
