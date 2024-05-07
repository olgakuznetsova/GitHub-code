#Libraries
library(MASS);
library(geometry);
library(grr);
#library(LogConcDEAD); #uncomment, if you want to initialise graphical_cgm2 using an estimate from LogConcDEAD

#Parameters CAN BE UPDATED
## Experiment parameters
path<-"" ; #directory where the data file is located and where the results will be saved
datafile<-"Wisconsin_breast_cancer_data.csv";
data<-read.csv(file.path(path,datafile));
intvaltol<-0.01 #how much the integral of the optimal function can deviate from density. 
                #In experiments on 10 points, we usually (but not always) get an optimal solution with integral 1+/- 0.005
                #If the integral is outside of this tolerance region, the optimisation is automatically restarted again 
                #with a initialisation

sourcefile<-"graphical-logconcave.R"
source(file.path(path,sourcefile))
### data sample
type<-"M" #"B" or "M"
toScale<-FALSE
npoints <- 10; #357 benign and 212 malignant instances
d <- 3;
w <- 1/npoints*rep(1, npoints);
graph_name=c("12-23") #for saving output file
cliques <- list(c(1,2),c(2,3));

### choose from the dataset which columns to consider
# The following are the possible parameters:
# "M_Radius","M_Texture","M_Perimeter","M_Area","M_Smoothness","M_Compactness","M_Concavity","M_Concave_Points","M_Symmetry","M_Fractal_Dimension",
# "SE_Radius","SE_Texture","SE_Perimeter","SE_Area","SE_Smoothness","SE_Compactness","SE_Concavity","SE_Concave_Points","SE_Symmetry","SE_Fractal_Dimension",
# "W_Radius","W_Texture","W_Perimeter","W_Area","W_Smoothness","W_Compactness","W_Concavity","W_Concave_Points","W_Symmetry","W_Fractal_Dimension"
vars <- c("M_Radius","M_Texture","M_Perimeter")
vars_name<-"RTP" #for saving output file


## Optimisation parameters
cl=list(maxit=800,trace=TRUE,REPORT=1, reltol=1e-6) #list of control parameters for the optimisation problem (see R optim))
                                                    #reltol controls convergence. Note that the default reltol in R is usually smaller ('1e-8'),
                                                    #however, in the experiments, setting reltol to '1e-8' did not improve the convergence of 
                                                    #integral to 1 (in some cases made even worse) or the value of the optimal function, 
                                                    #but considerably slowed down the optimisation algorithm. This is because in the final iterations,
                                                    #heights are very similar.

# DO NOT CHANGE BELOW
#Experiment
## Data prep
if (toScale){data<-scale(data[,vars])}

if (type =="B"){
  data<-data[2:(npoints+1),vars]
} else if (type =="M"){
  data<-data[359:(358+npoints),vars]
} else {
  print("Use B or M");
  q()
}
x <- as.matrix(data)

## optimization (with optimized algorithm graphical_cgm2)
## repeat the optimization if the optimal solution is too far from being a density
intval<-0
while(abs(intval-1)>intvaltol){
  tStart<-Sys.time();
  mle <- graphical_cgm2(cliques, x, w);
  tEnd<-Sys.time();
  print(tEnd-tStart);
  intval<-sigmafunc(cliques,x,mle$optimization$par,w,0.001) + mle$optimization$par %*% rep(w,length(cliques));
  print(intval)
}

# compute the tent poles and heights at tent poles
xy <- compute_xy(cliques,x,mle$optimization$par);
tentpoles <- xy$tentpoles;
y<-xy$y;

# save the result
LS = list(x, mle, tentpoles, y,intval);
filename=file.path(path,paste(type,npoints,graph_name,vars_name,paste(toString(Sys.Date()),".Rdata", sep=""),sep = "_", collapse=NULL));
save(LS, file=filename)
#q()