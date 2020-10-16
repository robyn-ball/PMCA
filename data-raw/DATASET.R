## code to prepare `DATASET` dataset goes here

load("X.RData")
load("Y.RData")
Xexample <- X
Yexample <- Y
usethis::use_data(Xexample)
usethis::use_data(Yexample)
