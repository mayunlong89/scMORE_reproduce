

#网址
#https://blog.csdn.net/nixiang_888/article/details/122744018
#https://bioconductor.org/packages/devel/bioc/vignettes/BiocStyle/inst/doc/AuthoringRmdVignettes.html

setwd("/Users/mayunlong/Documents/GitHub/")
#install.packages("devtools")
library(devtools)

#install.packages("roxygen2")
library(roxygen2)
library("usethis")
library(testthat)

usethis::create_package(path = "scMORE")
roxygen2::roxygenize("/Users/mayunlong/Documents/GitHub/scMORE")


usethis::use_package("fitdistrplus")
usethis::use_package("Pando")
usethis::use_package("Signac")
usethis::use_package("ArchR")
usethis::use_package("Seurat")
usethis::use_package("BSgenome.Hsapiens.UCSC.hg38")
usethis::use_package("tidyr")
usethis::use_package("dotgen")
usethis::use_package("GenomicRanges")
usethis::use_package("AnnotationDbi")
usethis::use_package("org.Hs.eg.db")
usethis::use_package("dplyr")
usethis::use_package("COSG")
usethis::use_package("IRanges")


devtools::load_all()
usethis::use_mit_license()
devtools::check() # 不需要任何参数

usethis::use_readme_rmd()
build_readme()

library(Pando)
usethis::use_data(phastConsElements20Mammals.UCSC.hg38, overwrite = TRUE)
usethis::use_data(motifs, overwrite = TRUE)
usethis::use_data(motif2tf, overwrite = TRUE)


#update NAMESPACE
devtools::document()


devtools::build() # 创建R包：*.tar.gz

# 需要重新启动R编辑器
devtools::install() # 安装R包

