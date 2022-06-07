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
npoints <- 10;
d <- 2;
x <- matrix(as.numeric(rnorm(npoints*d, mean=0, sd=1)), npoints, d);
w <- 1/npoints*rep(1, npoints);
cliques_ind <- list(c(1),c(2)); # for the disconnected graph on two vertices
cliques_complete <- list(c(1,2)); # for the complete graph on two vertices

# 1) optimization with our method and independence model
mle_ind_our <- graphical_cgm(cliques_ind, x, w, "BFGS", list(maxit=20000, trace=TRUE, REPORT=1));

# check that the integral is one, i.e., the optimal solution is a density
sigmafunc(cliques_ind, x, mle_ind_our$optimization$par, w, 0.001) + mle_ind_our$optimization$par %*% rep(w, length(cliques_ind));

# compute the tent poles and heights at tent poles
xy_ind_our <- compute_xy(cliques_ind, x, mle_ind_our$optimization$par);
tentpoles_ind_our <- xy_ind_our$tentpoles;
y_ind_our <- xy_ind_our$y;

# plot the tent function
g_ind_our <- interp_tent_mle(tentpoles_ind_our, list("par"=y_ind_our), 200);
M_ind_our = mesh(g_ind_our$x, g_ind_our$y);
scatterplot3d(c(M_ind_our$x), c(M_ind_our$y), c(g_ind_our$z), highlight.3d=TRUE);

# 2) optimization with our method and complete graph
mle_complete_our <- graphical_cgm(cliques_complete, x, w, "BFGS", list(maxit=20000, trace=TRUE, REPORT=1));

# check that the integral is one, i.e., the optimal solution is a density
sigmafunc(cliques_complete, x, mle_complete_our$optimization$par, w, 0.001) + mle_complete_our$optimization$par %*% rep(w, length(cliques_complete));

# compute the tent poles and heights at tent poles
xy_complete_our <- compute_xy(cliques_complete, x, mle_complete_our$optimization$par);
tentpoles_complete_our <- xy_complete_our$tentpoles;
y_complete_our <- xy_complete_our$y;

# plot the tent function
g_complete_our <- interp_tent_mle(tentpoles_complete_our, list("par"=y_complete_our), 200);
M_complete_our = mesh(g_complete_our$x, g_complete_our$y);
scatterplot3d(c(M_complete_our$x), c(M_complete_our$y), c(g_complete_our$z), highlight.3d=TRUE);

# 3) LogConcDEAD for the independence model
y_cliques_ind <- c();
for (i in 1:length(cliques_ind)) {
  xc <- x[,cliques_ind[[i]]]
  y_cliques_ind <- append(y_cliques_ind, mlelcd(xc)$logMLE)
};

# optimal value
sigmafunc(cliques_ind, x, y_cliques_ind, w, 0.001);
# value of the integral
sigmafunc(cliques_ind, x, y_cliques_ind, w, 0.001) + y_cliques_ind %*% rep(w, length(cliques_ind));

# tentpoles and heights of tentpoles
xy_ind_lcd <- compute_xy(cliques_ind, x, y_cliques_ind);
tentpoles_ind_lcd <- xy_ind_lcd$tentpoles;
y_ind_lcd <- xy_ind_lcd$y;

# plot the result
g_ind_lcd <- interp_tent_mle(tentpoles_ind_lcd, list("par"=y_ind_lcd), 200);
M_ind_lcd = mesh(g_ind_lcd$x, g_ind_lcd$y);
scatterplot3d(c(M_ind_lcd$x), c(M_ind_lcd$y), c(g_ind_lcd$z), highlight.3d=TRUE);

# 4) LogConcDEAD with complete graph
y_cliques_complete <- c();
for (i in 1:length(cliques_complete)) {
  xc <- x[,cliques_complete[[i]]]
  y_cliques_complete <- append(y_cliques_complete, mlelcd(xc)$logMLE)
};

# optimal value
sigmafunc(cliques_complete, x, y_cliques_complete, w, 0.001);
# value of the integral
sigmafunc(cliques_complete, x, y_cliques_complete, w, 0.001) + y_cliques_complete %*% rep(w, length(cliques_complete));

# tentpoles and heights of tentpoles
xy_complete_lcd <- compute_xy(cliques_complete, x, y_cliques_complete);
tentpoles_complete_lcd <- xy_complete_lcd$tentpoles;
y_complete_lcd <- xy_complete_lcd$y;

# plot the result
g_complete_lcd <- interp_tent_mle(tentpoles_complete_lcd,list("par"=y_complete_lcd), 200);
M_complete_lcd = mesh(g_complete_lcd$x, g_complete_lcd$y);
scatterplot3d(c(M_complete_lcd$x), c(M_complete_lcd$y), c(g_complete_lcd$z), highlight.3d=TRUE);

# compute distances from true density
# 1) "(tentpoles_ind_our, y_ind_our)"
dist_ind_our <- squared_L2(standard2dGaussian, tentFunction(tentpoles_ind_our, y_ind_our), -10, 10, 50);
# 2) "(tentpoles_complete_our, y_complete_our)"
dist_complete_our <- squared_L2(standard2dGaussian, tentFunction(tentpoles_complete_our, y_complete_our), -10, 10, 50);
# 3) "(x, y_cliques_ind)"
dist_ind_lcd <- squared_L2(standard2dGaussian, tentFunction(tentpoles_ind_lcd, y_ind_lcd), -10, 10, 50);
# 4) "(x, y_cliques_complete)"
dist_complete_lcd <- squared_L2(standard2dGaussian, tentFunction(tentpoles_complete_lcd, y_complete_lcd), -10, 10, 50);

# save the results
dir_name="<add directory here>";
# use another dir_name if necessary

L_ind_our = list(x, mle_ind_our, tentpoles_ind_our, y_ind_our, g_ind_our, dist_ind_our);
filename_ind_our = paste(dir_name,"npoints",npoints,"graph-12-ind-our-","date",toString(Sys.Date()), ".Rdata", sep="", collapse=NULL);
save(L_ind_our, file = filename_ind_our);

L_complete_our = list(x, mle_complete_our, tentpoles_complete_our, y_complete_our, g_complete_our, dist_complete_our);
filename_complete_our = paste(dir_name,"npoints",npoints,"graph-12-complete-our-","date",toString(Sys.Date()), ".Rdata", sep="", collapse=NULL);
save(L_complete_our, file = filename_complete_our);

L_ind_lcd = list(x, y_cliques_ind, tentpoles_ind_lcd, y_ind_lcd, g_ind_lcd, dist_ind_lcd);
filename_ind_lcd = paste(dir_name,"npoints",npoints,"graph-12-ind-lcd-","date",toString(Sys.Date()), ".Rdata", sep="", collapse=NULL);
save(L_ind_lcd, file = filename_ind_lcd);

L_complete_lcd = list(x, y_cliques_complete, tentpoles_complete_lcd, y_complete_lcd, g_complete_lcd, dist_complete_lcd);
filename_complete_lcd = paste(dir_name,"npoints",npoints,"graph-12-complete-lcd-","date",toString(Sys.Date()), ".Rdata", sep="", collapse=NULL);
save(L_complete_lcd, file = filename_complete_lcd);
q()