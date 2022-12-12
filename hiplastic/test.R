

library(dplyr)
library(ggplot2)
library(reshape2)
library(tidyr)
library(data.table)
library(Matrix)
library(RColorBrewer)
library(dcon)


load('/rafalab/lzou/dcon/hiplastic/shiny_data.RData')
b1 <- summary(b1_chr18) |> mutate(i=i*1e5, j=j*1e5)
l1 <- summary(l1_chr18) |> mutate(i=i*1e5, j=j*1e5)
tot <- summary(total_chr18) |> mutate(i=i*1e5, j=j*1e5)


r1 <- c(35.5e6,39.3e6)
r2 <- c(61e6,66e6)

r1 <- r1/1e5
r2 <- r2/1e5
loop_coords <- c(r1[1], r1[2], r2[1], r2[2])

little_L1 <- as.matrix(l1_chr18)[loop_coords[1]:loop_coords[2],
                                 loop_coords[3]:loop_coords[4]]

d <- as.matrix(total_chr18)[loop_coords[1]:loop_coords[2],
                                    loop_coords[1]:loop_coords[2]]
D_r1 <- d+t(d)-diag(diag(d))

d <- as.matrix(total_chr18)[loop_coords[3]:loop_coords[4],
                            loop_coords[3]:loop_coords[4]]
D_r2 <- d+t(d)-diag(diag(d))

offset_y <- rowSums(D_r1)
offset_x <- rowSums(D_r2)


y <- rowSums(little_L1)
D <- D_r1
D <- dcon:::normalize_hic(D, gamma=0.5)
B <- dcon:::construct_basis(1:length(y), df = length(y))
fit <- dcon:::fit_decon(y, D, df = length(y), offset=offset_y)
fit2 <- dcon:::fit_decon(y, D, df = length(y))
a <- fit$a
a2 <- fit2$a


dcon:::loglik(a, y, D, B)

####### loglik function:

if (is.null(offset)) {
  eBa <- exp(B%*%a)
} else {
  eBa <- offset*exp(B%*%a)
}
DeBa <- D%*%eBa
first_term <- sum(y*log(DeBa))
second_term <- sum(DeBa)
third_term <- sum(lfactorial(y))
return(-first_term+second_term+third_term)






