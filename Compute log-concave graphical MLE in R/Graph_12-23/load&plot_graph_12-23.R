#install.packages("MASS");
#install.packages("geometry");
#install.packages("scatterplot3d");
#install.packages("plot3D");
#install.packages("rgl");
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
# Example for the graph {(1,2),(2,3)}
###################################

# data sample
# First, remember to load the dataset LS = list(x, mle, tentpoles, y)

npoints <- 60;
d <- 3;
x <- LS[[1]];
w <- 1/npoints*rep(1, npoints);
cliques <- list(c(1,2),c(2,3));

# load optimization
mle <- LS[[2]];

# check that the integral is one, i.e., the optimal solution is a density
sigmafunc(cliques,x,mle$optimization$par,w,0.001) + mle$optimization$par %*% rep(w,length(cliques));

# load tent poles and heights at tent poles
tentpoles <- LS[[3]];
y <- LS[[4]];

# compute the support S_{G,X}
SGX <- compute_SGX_vertices(cliques,x);

# restriction to the clique 12
x12 = x[,c(1,2)];
plot(x12, xlab="", ylab="");
y12 = mle$optimization$par[1:npoints];

# plot the tent function on the clique 12
g12 <- interp_tent_mle(x12,list("par"=y12),200);
M12 = mesh(g12$x, g12$y);
scatterplot3d(c(M12$x), c(M12$y), c(g12$z), highlight.3d=TRUE, xlab="", ylab="", zlab="");

# restriction to the clique 23
x23 = x[,c(2,3)];
plot(x23, xlab="", ylab="");
y23 = mle$optimization$par[(npoints+1):(2*npoints)];

# plot the tent function on the clique 23
g23 <- interp_tent_mle(x23,list("par"=y23),200);
M23 = mesh(g23$x, g23$y);
scatterplot3d(c(M23$x), c(M23$y), c(g23$z), highlight.3d=TRUE, xlab="", ylab="", zlab="");

q()