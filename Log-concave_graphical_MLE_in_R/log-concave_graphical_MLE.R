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
source("<add directory of the file 'functions.R'>")

# data sample
npoints <- 10;
d <- 2;
x <- matrix(as.numeric(rnorm(npoints*d, mean=0, sd=1)), npoints, d);
# weights
w <- 1/npoints*rep(1, npoints);
# maximal cliques of the graph G
cliques <- list(c(1),c(2));

# compute the log-concave graphical MLE
mle <- graphical_cgm(cliques, x, w, "BFGS", list(maxit=20000, trace=TRUE, REPORT=1));

# check that the integral is one, i.e., the optimal solution is a density
sigmafunc(cliques, x, mle$optimization$par, w, 0.001) + mle$optimization$par %*% rep(w, length(cliques));

# compute the tent poles and the heights at tent poles
xy <- compute_xy(cliques, x, mle$optimization$par);
tentpoles <- xy$tentpoles;
y <- xy$y;

# plot the tent function
g <- interp_tent_mle(tentpoles, list("par"=y), 200);
M = mesh(g$x, g$y);
scatterplot3d(c(M$x), c(M$y), c(g$z), highlight.3d=TRUE);

# compute the distance from the true density (if the data was sampled from a standard normal distribution)
dist <- squared_L2(standard2dGaussian, tentFunction(tentpoles, y), -10, 10, 50);

# save the results
dir_name="<add directory here>";
# use another dir_name if necessary

L = list(x, mle, tentpoles, y, g, dist);
filename = paste(dir_name,"npoints",npoints,"graph-G-","date",toString(Sys.Date()), ".Rdata", sep="", collapse=NULL);
save(L, file = filename);
q()