
set.seed(12345)

#===============================================================
# load all the necessary libraries
#===============================================================
needed.packages <- c("evd", "Matrix", "fdrtool", "shape", "corpcor", "RColorBrewer",
"doParallel", "foreach", "caTools", "stringr", "abind", "propagate", "here",
"testthat", "regress", "tidyverse", "qtl2", "tidyr", "qtl2convert", "devtools",
"yaml","jsonlite", "data.table", "Rcpp", "RcppEigen", "RSQLite", "qtl", "shape")
for(i in 1:length(needed.packages)){library(needed.packages[i], character.only = TRUE)}


#===============================================================
# source code
#===============================================================
cape.fun <- list.files(file.path("..", "capempp"), full.names = TRUE)
for(i in 1:length(cape.fun)){source(cape.fun[i])}

#===============================================================
# read in data and convert to new cape format
#===============================================================
old.data.obj <- read.population("CAPE_data.csv", id.col = 1)
cape.obj <- cape2mpp(old.data.obj)

data.obj <- cape.obj$data.obj
geno.obj <- cape.obj$geno.obj

#===============================================================
# run.cape
#===============================================================
final.obj <- run.cape(data.obj, geno.obj, path = "results")
