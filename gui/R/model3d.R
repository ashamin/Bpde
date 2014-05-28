# call via: Rscript model3d.R model.csv (does not work for now)

library("MASS")
library("lattice")
library("plyr")
library("emdbook")
library("rgl")
library("fields")

# args <- commandArgs(trailingOnly = TRUE)
csv = read.csv2(args[1])
x = csv$fx[csv$fx>=0]
y = csv$fy[csv$fy>=0]
H = matrix(0, length(x), length(y))
for (j in 1:length(y)) {
  for (i in 1:length(x)) {
    H[i,j] = csv$fH[(j-1)*length(x)+(i-1) + 1]
  }
}

open3d()
bg3d("white")
material3d(col="black")
persp3d(x, y, H, aspect=c(1, 1, 0.5), col = "lightblue",
        xlab = "X", ylab = "Y")
