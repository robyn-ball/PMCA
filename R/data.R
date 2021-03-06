#' @title Example data for PMCA
#'
#' @description
#'
#' @format Two coupled matrices, Xexample (6 x 28) and Yexample (1000 x 28)
#' \describe{
#'   \item{Xexample}{A matrix of substage proportions for each mouse (28 columns)
#'   and substage (6 rows). The rownames of Xexample are the 6 substage names.}
#'   \item{Yexample}{A matrix of gene expression values for each mouse (28 columns)
#'   and gene (1000 rows). The rownames of Yexample are the 1000 gene names.}
#' }
#' @source Ball, R.L., Fujiwara, Y., Sun, F. et al.
#' Regulatory complexity revealed by integrated cytological and RNA-seq analyses
#' of meiotic substages in mouse spermatocytes. BMC Genomics 17, 628 (2016).
#' <https://doi.org/10.1186/s12864-016-2865-1>
"Xexample"
