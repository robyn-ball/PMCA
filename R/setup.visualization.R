#' Simple wrapper for setup functions
#'
#' @export
#'
#' @examples see README

setup.visualization <- function(){
  
  library(reticulate) # We use reticulate to run Python functions in R
  library(plotly) # Plotly is the visualization tool
  
  # Set up the Python virtual environment
  reticulate::use_virtualenv("iggi")
  
  # Load the Python module with reticulate
  kp <- import("kpartite")
  
}