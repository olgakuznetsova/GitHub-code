# First load the '.Rdata' files

#install.packages("MASS");
#install.packages("geometry");
#install.packages("scatterplot3d");
#install.packages("plot3D");
#install.packages("rgl");
#install.packages("manipulateWidget");
#install.packages("matlib");
#install.packages("LogConcDEAD");
library(MASS)
library(geometry)
library(scatterplot3d)
library(plot3D)
library(rgl)
library(matlib)
library(LogConcDEAD)

# import all functions needed
source("<add directory of the file 'functions.R'>/functions.R")

###################################
# Example for the graph {(1,2)}
# comparing with LOGCONCDEAD on the data sample
###################################

# data sample
npoints <- 50;
d <- 2;
x <- L_ind_our[[1]];
w <- 1/npoints*rep(1, npoints);
cliques_ind <- list(c(1),c(2)); # for the disconnected graph on two vertices
cliques_complete <- list(c(1,2)); # for the complete graph on two vertices

# 1) load mle optimization with our method and independence model
mle_ind_our <- L_ind_our[[2]];

# check that the integral is one, i.e., the optimal solution is a density
sigmafunc(cliques_ind, x, mle_ind_our$optimization$par, w, 0.001) + mle_ind_our$optimization$par %*% rep(w, length(cliques_ind));

# load tent poles and heights at tent poles
tentpoles_ind_our <- L_ind_our[[3]];
y_ind_our <- L_ind_our[[4]];

# plot the tent function
g_ind_our <- L_ind_our[[5]];
M_ind_our = mesh(g_ind_our$x, g_ind_our$y);
scatterplot3d(c(M_ind_our$x), c(M_ind_our$y), c(g_ind_our$z), highlight.3d=TRUE, xlab="", ylab="", zlab="");

# 2) load mle optimization with our method and complete graph
mle_complete_our <- L_complete_our[[2]];

# check that the integral is one, i.e., the optimal solution is a density
sigmafunc(cliques_complete, x, mle_complete_our$optimization$par, w, 0.001) + mle_complete_our$optimization$par %*% rep(w, length(cliques_complete));

# compute the tent poles and heights at tent poles
tentpoles_complete_our <- L_complete_our[[3]];
y_complete_our <- L_complete_our[[4]];

# plot the tent function
g_complete_our <- L_complete_our[[5]];
M_complete_our = mesh(g_complete_our$x, g_complete_our$y);
scatterplot3d(c(M_complete_our$x), c(M_complete_our$y), c(g_complete_our$z), highlight.3d=TRUE, xlab="", ylab="", zlab="");


# 3) LogConcDEAD for the independence model
y_cliques_ind <- L_ind_lcd[[2]];

# optimal value
sigmafunc(cliques_ind, x, y_cliques_ind, w, 0.001);
# value of the integral
sigmafunc(cliques_ind, x, y_cliques_ind, w, 0.001) + y_cliques_ind %*% rep(w, length(cliques_ind));

# tentpoles and heights of tentpoles
tentpoles_ind_lcd <- L_ind_lcd[[3]];
y_ind_lcd <- L_ind_lcd[[4]];

# plot the result
g_ind_lcd <- L_ind_lcd[[5]];
M_ind_lcd = mesh(g_ind_lcd$x, g_ind_lcd$y);
scatterplot3d(c(M_ind_lcd$x), c(M_ind_lcd$y), c(g_ind_lcd$z), highlight.3d=TRUE, xlab="", ylab="", zlab="");


# 4) LogConcDEAD with complete graph
y_cliques_complete <- L_complete_lcd[[2]];

# optimal value
sigmafunc(cliques_complete, x, y_cliques_complete, w, 0.001);
# value of the integral
sigmafunc(cliques_complete, x, y_cliques_complete, w, 0.001) + y_cliques_complete %*% rep(w, length(cliques_complete));

# tentpoles and heights of tentpoles
tentpoles_complete_lcd <- L_complete_lcd[[3]];
y_complete_lcd <- L_complete_lcd[[4]];

# plot the result
g_complete_lcd <- L_complete_lcd[[5]];
M_complete_lcd = mesh(g_complete_lcd$x, g_complete_lcd$y);
scatterplot3d(c(M_complete_lcd$x), c(M_complete_lcd$y), c(g_complete_lcd$z), highlight.3d=TRUE,xlab="", ylab="", zlab="");

q()
