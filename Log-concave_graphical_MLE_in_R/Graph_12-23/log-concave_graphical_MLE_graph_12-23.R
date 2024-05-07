#install.packages("MASS");
#install.packages("geometry");
#install.packages("scatterplot3d");
#install.packages("plot3D");
#install.packages("rgl");
#install.packages("manipulateWidget");
#install.packages("matlib");
#install.packages("LogConcDEAD");
library(MASS);
library(geometry);
library(rgl);
library(LogConcDEAD);

# import all functions needed
source("<add directory of the file 'functions.R'>")

###################################
# Example for the graph {(1,2),(2,3)}
###################################

# data sample
npoints <- 10;
d <- 3;
x <- matrix(as.numeric(rnorm(npoints*d, mean=0, sd=5)), npoints, d);
w <- 1/npoints*rep(1, npoints);
cliques <- list(c(1,2),c(2,3));

# optimization
mle <- graphical_cgm(cliques, x, w, "BFGS", list(maxit=20000,trace=TRUE,REPORT=1));

# check that the integral is one, i.e., the optimal solution is a density
sigmafunc(cliques,x,mle$optimization$par,w,0.001) + mle$optimization$par %*% rep(w,length(cliques));

# compute the tent poles and heights at tent poles
xy <- compute_xy(cliques,x,mle$optimization$par);
tentpoles <- xy$tentpoles;
y<-xy$y;

# restriction to the clique 12
x12=x[,c(1,2)];
y12=mle$optimization$par[1:npoints];

# restriction to the clique 23
x23=x[,c(2,3)];
y23=mle$optimization$par[(npoints+1):(2*npoints)];

# save the results
dir_name="<add directory here>";
# use another dir_name if necessary

LS = list(x, mle, tentpoles, y);
filename=paste(dir_name,"npoints",npoints,"graph-12-23-","date",toString(Sys.Date()), ".Rdata", sep="", collapse=NULL);
save(LS, file=filename);
q()
