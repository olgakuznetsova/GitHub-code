# 1-dim versions of geometry package functions
# An inequality x+a>=0 is in the matrix format  (-x -a)

# Compute halfspace intersection about a point
# p is an M-by-(N+1) matrix. Each row of p represents a halfspace by a N-dimensional normal
# to a hyperplane and the offset of the hyperplane.
# fp is a “feasible” point that is within the space contained within all the halfspaces.
halfspacenNew <- function(p,fp,options="Tv") {
  n <- nrow(p)
  d <- ncol(p)
  # 1-dimensional case is not considered by the geometry package
  if (d==2) {
    minv <- c()
    maxv <- c()
    for (i in 1:n) {
      if(p[i,1] > 0) {
        maxv <- append(maxv,-p[i,2]/p[i,1])
      } else if (p[i,1] < 0) {
        minv <- append(minv,-p[i,2]/p[i,1])
      }
    }
    v <- matrix(c(min(maxv),max(minv)),nrow=2)
    return(v)
  } else {
    return(halfspacen(p,fp,options))
  }
};

# Returns information about the smallest convex complex of a set of input points in N-dimensional space 
# (the convex hull of the points). By default, indices to points forming the facets of the hull are returned; 
# optionally normals to the facets and the generalised surface area and volume can be returned.
# p is an M-by-N matrix. The rows of p represent M points in N-dimensional space.
convhullnNew <- function(p, options = "Tv", output.options = NULL, return.non.triangulated.facets = FALSE) {
  d <- ncol(p)
  n <- nrow(p)
  if (d==1) {
    maxValue <- max(p)
    minValue <- min(p)
    maxVertices <- c()
    minVertices <- c()
    for (i in 1:n) {
      if(p[i,1]==maxValue) {
        maxVertices <- append(maxVertices,i)
      }
      if(p[i,1]==minValue) {
        minVertices <- append(minVertices,i)
      }
    }
    chull <- matrix(c(maxVertices,minVertices),nrow=2,byrow=TRUE)
    if(is.null(output.options)) {
      return(chull)
    } else {
      maxNormal <- c(1,-maxValue)
      minNormal <- c(-1,minValue)
      normals <- matrix(c(maxNormal,minNormal),nrow=2,byrow=TRUE)
      output <- list("hull" = chull, "normals"  = normals, "p" = p)  
      return(output)
    }
  } else {
    return(convhulln(p,options,output.options,return.non.triangulated.facets))
  }
};

# Given the Cartesian coordinates of one or more points, 
# compute the barycentric coordinates of these points with respect to a simplex.
# X is the reference simplex in N dimensions represented by a (N+1)-by-N matrix
# P is an M-by-N matrix in which each row is the Cartesian coordinates of a point.
# Output: M-by-(N+1) matrix in which each row is the barycentric coordinates of 
# corresponding row of P. If the simplex is degenerate a warning is issued and the function returns NULL.
cart2baryNew <- function(X,P) {
  d <- ncol(P)
  if(d==1) {
    beta1 <- (P-X[2,1])/(X[1,1]-X[2,1])
    beta2 <- 1-beta1
    return(cbind(beta1,beta2))
  } else {
    return(cart2bary(X,P))
  }
};

# Computes (an approximation of) the integral of 
# exp(tent function) above a simplex where the tent function is constant
integral_appr <- function(y) {
  d = length(y)
  ymean = mean(y)
  z = y - ymean
  tmp = 1
  for (i in 1:d){
    tmp = tmp + (z[i]*z[i])/(2*(d)*(d+1)) + (z[i]*z[i]*z[i])/
      (3*(d)*(d+1)*(d+2));
  }
  for (i in 1:(d-1)) {
    tmp = tmp / i;
  }
  tmp = tmp*exp(ymean);
  return(tmp);
};

# Computes (an approximation of) the integral of exp(tent function).
# This is based on the formula from Robeva, Sturmfels, Uhler
integral <- function(y, eps) {
  y = sort(y);
  d = length(y);
  tmp = 0;
  if(y[d] - y[1] < eps) {
    return(integral_appr(y));
  } else {
    if(d == 2) {
      tmp = (exp(y[1]) - exp(y[2]))/(y[1] - y[2]);
      return(tmp);
    } else {
      e1 = integral(y[2:d], eps);
      e2 = integral(y[1:(d-1)], eps);
      return((e1-e2)/(y[d] - y[1]));
    }
  }
};

# Computes the derivative of the integral of exp(tent function)
integral_deriv <- function(y, i, eps) {
  z = c(y, y[i])
  return(integral(z, eps));
};

# Computes the upper convex hull of points (x_i,y_i)
# Output are the normals with offsets for the upper convex hull
# This has been tested on the test example
upper_chull_normals <- function(x,y,return_non_triangulated_facets) {
  output <- c()
  Q = cbind(x,y)
  affine_hull <- Null(t(Q) - Q[1,])
  if (ncol(affine_hull) > 0) {
    #Q = cbind(x,y+rnorm(length(y), mean=0, sd=1))
    Cupper <- matrix(cbind(t(affine_hull),-Q[1,,drop=FALSE]%*%affine_hull),nrow=1)
    if (Cupper[ncol(Q)] < 0 ) {
      Cupper <- -Cupper
    }
    Cupper_simplices <- matrix(seq(1:nrow(x)),nrow=1)
    output$normals <- Cupper
    output$simplices <- Cupper_simplices
  } else {
    # compute convex hull including normals to facets
    C2 = convhullnNew(Q,output.options = "n",return.non.triangulated.facets = return_non_triangulated_facets)
    Cupper = 0
    Cupper_simplices = 0
    # return normals whose y coordinate is positive
    for (i in 1:nrow(C2$normals)) {
      K = C2$normals[i,]
      if (K[ncol(Q)] > 1e-08) {
        Cupper = rbind(Cupper, K)
        Cupper_simplices = rbind(Cupper_simplices, C2$hull[i,])
      }
    }
    output$normals <- Cupper[2:nrow(Cupper),,drop=FALSE]
    output$simplices <- Cupper_simplices[2:nrow(Cupper_simplices),,drop=FALSE]
  }
  return(output);
};

# lift normals on a clique to normals in the full space
lift_normals <- function(n,clique,normals) {
  normals_lifted <- matrix(0, nrow(normals), n+1)
  normals_lifted[,c(clique,n+1)] <- normals
  return(normals_lifted)
};

# upper convex hull
sum_of_tent_functions_normals<-function(cliques,x,y_cliques) {
  n <- ncol(x)
  normals_for_all_cliques <- list()
  simplices_for_all_cliques <- list()
  # for each of the cliques find facet normals
  for (i in 1:length(cliques)) {
    clique <- cliques[[i]]
    clique_size <- length(clique)
    xc <- x[,clique,drop = FALSE]
    upper_chull <- upper_chull_normals(xc,y_cliques[((i-1)*nrow(x)+1):(i*nrow(x))],TRUE)
    normals <- upper_chull$normals
    # remove the y coefficient
    normals <- normals / normals[,clique_size+1]
    normals <- normals[,-(clique_size+1)]
    # the previous converts a one-row matrix to a vector
    if(!is.matrix(normals)) {
      normals <- matrix(normals,nrow = 1)
    }
    # lift the normal from clique coordinates to n+1 coordinates
    normals_lifted <- lift_normals(n,clique,normals)
    normals_for_all_cliques[[i]] <- normals_lifted
    simplices_for_all_cliques[[i]] <- upper_chull$simplices
  }
  # construct the normals for S_{G,X}
  # by taking all possible sums of normals for grids
  n_facets <- 1
  for (i in 1:length(cliques)) {
    n_facets <- n_facets * nrow(normals_for_all_cliques[[i]])
  }
  normals_SGX <- matrix(0,n_facets,n+1)
  facet_seq <- list()
  for (i in 1:length(cliques)){
    facet_seq[[i]] <- seq(1,nrow(normals_for_all_cliques[[i]]))
  }
  grid <- expand.grid(facet_seq)
  for (i in 1:nrow(grid)) {
    for (j in 1:length(cliques)) {
      normals_SGX[i,] <- normals_SGX[i,] + normals_for_all_cliques[[j]][grid[i,j],]
    }
  }
  output <- c()
  output$normals <- normals_SGX
  # record which heights contribute to which facets on cliques
  output$simplices <- list()
  for (i in 1:ncol(grid)) {
    output$simplices[[i]] <- t(sapply(grid[,i],function(z) simplices_for_all_cliques[[i]][z,]))
  }
  return(output)
};

# compute the support S_{G,X}
compute_SGX_vertices <- function(cliques,x) {
  n <- ncol(x)
  normals <- matrix(0 , nrow = 1, ncol = n+1)
  for (i in 1:length(cliques)) {
    clique <- cliques[[i]]
    xc <- x[,clique,drop=FALSE]    
    if(!is.matrix(xc)) {
      xc <- matrix(xc,ncol=1)
    }
    if (ncol(xc) == 1) {
      xc_normals <- matrix(c(-1,1,min(xc),-max(xc)), ncol = 2)
    } else {
      # the support for clique
      Sc <- convhullnNew(xc, output.options = "n",return.non.triangulated.facets = TRUE)
      xc_normals <- Sc$normals
    }
    # lifted normals of the support for clique
    normals <- rbind(normals,lift_normals(n,clique,xc_normals))
  }
  normals <- normals[2:nrow(normals),,drop=FALSE]
  SGX_vertices <- halfspacenNew(normals,(1/nrow(x))*colSums(x))
  return(SGX_vertices)
};

# list that gives vertices in a facet
vertices_in_facets <- function(upper_hull_normals,vertices) {
  # add a column of ones
  vertices <- cbind(vertices,rep(1,nrow(vertices)))
  eval_facets <- upper_hull_normals %*% t(vertices)
  vertices_in_facets_mat <- (abs(eval_facets) < 1e-06)
  vertices_in_facets_list <- list()
  for (i in 1:nrow(vertices_in_facets_mat)) {
    vertices_in_facet <- c()
    for (j in 1:ncol(vertices_in_facets_mat)) {
      if(vertices_in_facets_mat[i,j]) {
        vertices_in_facet <- append(vertices_in_facet,j)
      }
    }
    vertices_in_facets_list[[i]] <- vertices_in_facet
  }
  return(vertices_in_facets_list)
};

# function that computes the vertices of the upper convex hull on S_{G,X}
# from the heights on each of the cliques
upper_chull_vertices_from_cliques<-function(cliques,x,y_cliques) {
  n <- ncol(x)
  SGX <- compute_SGX_vertices(cliques,x)
  # equations for S_{G,X}
  SGX_normals <- convhullnNew(SGX,output.options = "n",return.non.triangulated.facets = TRUE)$normals
  SGX_normals_lifted <-matrix(0,nrow=nrow(SGX_normals),ncol=n+2)
  SGX_normals_lifted[,c(seq(1:n),n+2)]<-SGX_normals
  # equations for upper hull on S_{G,X}
  sum_of_tent_functions <- sum_of_tent_functions_normals(cliques,x,y_cliques)
  upper_hull_normals <- sum_of_tent_functions$normals
  upper_hull_normals_lifted <- matrix(0,nrow=nrow(upper_hull_normals),ncol=n+2)
  upper_hull_normals_lifted[,c(seq(1:n),n+2)] <- upper_hull_normals
  upper_hull_normals_lifted[,n+1]<- 1
  # equation to bound the polytope from below
  y_min <- min(0,y_cliques)
  y_bound <- y_min * length(cliques)
  y_ineq <- matrix(0,nrow=1,ncol=n+2)
  y_ineq[1,n+1] <- -1
  y_ineq[1,n+2] <- y_bound - 10
  lifted_SGX_normals <- rbind(SGX_normals_lifted,upper_hull_normals_lifted,y_ineq)
  lifted_SGX <- halfspacenNew(lifted_SGX_normals,c((1/nrow(x))*colSums(x),y_bound),options="Q12")
  lifted_and_cropped_SGX <- lifted_SGX[lifted_SGX[,ncol(lifted_SGX)]>(y_bound - 0.5),,drop=FALSE]
  vert_in_facets <-vertices_in_facets(upper_hull_normals_lifted,lifted_and_cropped_SGX)
  output <- c()
  output$vertices <- lifted_and_cropped_SGX
  output$vertices_in_facets <- vert_in_facets
  output$decomposition <- sum_of_tent_functions$simplices
  return(output)
};



# Input: cliques of the graph, data points x, heights on each clique
# Output: list that consists of tentpoles and heights at tentpoles
compute_xy <- function(cliques,x,y_cliques) {
  output <- c()
  lifted_SGX <- upper_chull_vertices_from_cliques(cliques,x,y_cliques)
  lifted_vertices <- lifted_SGX$vertices
  cols <- ncol(lifted_vertices)
  output$tentpoles <- lifted_vertices[,-cols,drop=FALSE]
  output$y <- matrix(lifted_vertices[,cols],ncol=1)
  output$vertices_in_facets <- lifted_SGX$vertices_in_facets
  output$decomposition <- lifted_SGX$decomposition
  return(output)
};

# Input: tentpoles x_i and heights at the tentpoles y_i
# Output: matrix, where each row corresponds to a simplex of 
# the regular triangulation given by the
# upper convex hull of points (x_i,y_i)
# the row entries are indices of the vertices in the simplex
upper_chull_triang <- function(x,y) {
  Q = cbind(x,y)
  # if Q is not full dimensional, then add noise to the heights
  if (ncol(Null(t(Q) - Q[1,])) > 0) {
    Q = cbind(x,y+rnorm(length(y), mean=0, sd=1))
  }
  C = convhullnNew(Q,output.options = "n",return.non.triangulated.facets = FALSE)
  Cupper = 0
  # return simplices whose y coordinate is positive
  for (i in 1:nrow(C$normals)) {
    K = C$normals[i,]
    if (K[ncol(Q)] > 1e-08) {
      Cupper = rbind(Cupper, C$hull[i,])
    }
  }
  return(Cupper[2:nrow(Cupper),,drop=FALSE]);
};

# Computes the objective function to be minimized
sigmafunc <- function(cliques, x, y_cliques, w, eps) {
  # find the tentpoles and heights at the tentpoles
  xy <- compute_xy(cliques,x,y_cliques)
  tentpoles <- xy$tentpoles
  y <- xy$y
  # find the triangulation corresponding to the upper hull
  triang <- upper_chull_triang(tentpoles,y)
  # first part of the objective function
  fval <- 0
  fval <- fval - y_cliques%*% rep(w,length(cliques));
  #computing the integral
  for (i in 1:nrow(triang)) {
    integ <- integral(y[triang[i,]], eps);
    A <- cbind(tentpoles[triang[i,],], rep(1,ncol(x)+1))
    fval <- fval + integ*abs(det(A));
  }
  return(fval);
};

# Returns the objective function to be minimized
# for specified values of x, w, epsilon
sigmafuncS <- function(cliques, x, w, eps) {
  s <- function(y_cliques){
    return(sigmafunc(cliques, x, y_cliques, w, eps));
  }
  return(s);
};


# Computes the derivative of the objective funcion
gradsigmafunc <- function(cliques, x, y_cliques, w, eps) {
  n <- ncol(x)
  # find the tentpoles and heights at the tentpoles
  xy <- compute_xy(cliques,x,y_cliques)
  tentpoles <- xy$tentpoles
  y <- xy$y
  vert_in_facets <- xy$vertices_in_facets
  decomp <- xy$decomposition
  # find the triangulation corresponding to the upper hull
  triang <- upper_chull_triang(tentpoles,y)
  gradval = -rep(w,length(cliques));
  # Gradients of the integral
  for (simplex_ix in 1:nrow(triang)) {
    # Find which original facet/cell the simplex corresponds to
    facet_ix <- 0 
    for (i in 1:length(vert_in_facets)) {
      if (all(triang[simplex_ix,] %in% vert_in_facets[[i]])) {
        facet_ix <- i
      }
    }
    if(facet_ix > 0) {
      # Loop over all vertices in the simplex 
      # (the main part of the gradient will depend on it)
      for (i in 1:ncol(triang)) {
        tentpole_ix <- triang[simplex_ix,i]
        integ <- integral_deriv(y[triang[simplex_ix,]], i, eps);
        A <- cbind(tentpoles[triang[simplex_ix,],], rep(1,ncol(tentpoles)+1))
        for (clique_ix in 1:length(cliques)) {
          # Find the cell on a clique that contains the projected small simplex
          small_cell <- decomp[[clique_ix]][facet_ix,]
          # Since this might not be a simplex, we might have to triangulate it further (there might be a mistake here)
          # But hopefully this will almost never happen
          if(length(small_cell)>length(cliques[[clique_ix]])+1) {
            triang_cell <- upper_chull_normals(x[small_cell,cliques[[clique_ix]]],rnorm(length(small_cell), mean=0, sd=1),FALSE)$simplices
            for(j in 1:nrow(triang_cell)) {
              if(inhulln(convhulln(x[small_cell[triang_cell[j,]],cliques[[clique_ix]],drop = FALSE]),matrix(tentpoles[tentpole_ix,cliques[[clique_ix]]],nrow=1))) {
                small_simplex <- small_cell[triang_cell[j,]]
                print(small_simplex)
                print(small_cell)
              }
            } 
          } else {
            small_simplex <- small_cell
          }
          # barycentric coordinates for the projection of the tent pole in a simplex on a clique
          beta <- cart2baryNew(x[small_simplex,cliques[[clique_ix]],drop = FALSE],matrix(tentpoles[tentpole_ix,cliques[[clique_ix]]],nrow=1))
          for (j in 1:length(small_simplex)) {
            gradval[(clique_ix-1)*nrow(x)+small_simplex[j]] = gradval[(clique_ix-1)*nrow(x)+small_simplex[j]] + integ*abs(det(A))*beta[j];
          }
        }
      }
    }
  }
  return(gradval);
};

# Returns the integral of the objective function
# for specified values of x, w, epsilon
gradsigmafuncS <- function(cliques, x, w, eps){
  s <- function(y_cliques){
    return(gradsigmafunc(cliques, x, y_cliques, w, eps));
  }
  return(s);
};


# Computes the graphical and log-concave MLE using the conditional gradient method
graphical_cgm <- function(cliques, x, w, method = "BFGS", control=list(maxit=800,trace=TRUE,REPORT=1)) {
  eps <- 0.001
  # LogConcDEAD initialization of y_C
  y_cliques <- c()
  for (i in 1:length(cliques)) {
    xc <- x[,cliques[[i]]]
    y_cliques <- append(y_cliques, mlelcd(xc,w)$logMLE)
  }
  # random initialization of y_C (comment if not used)
  # y_cliques <- rnorm(length(cliques)*nrow(x), mean=0, sd=3)
  S <- optim(y_cliques, sigmafuncS(cliques, x, w, eps), gradsigmafuncS(cliques, x, w, eps), method = method ,control = control)
  output <- c()
  output$optimization <- S
  return(output)
};

# Prepares the computed MLE for plotting
interp_tent_mle <- function(tentpoles, S, meshpar) {
  min_coord = rep(0, ncol(tentpoles));
  max_coord = rep(0, ncol(tentpoles));
  for(i in 1:ncol(tentpoles)){
    min_coord[i] = min(tentpoles[,i]);
    max_coord[i] = max(tentpoles[,i]);
  }
  x = seq(min_coord[1]-0.1, max_coord[1]+0.1, length.out = meshpar);
  y = seq(min_coord[2]-0.1, max_coord[2]+0.1, length.out = meshpar);
  z = matrix(0, ncol=meshpar, nrow=meshpar);
  triangs = upper_chull_triang(tentpoles, S$par);
  for (i in 1:length(x)) {
    for (j in 1:length(y)) {
      triang_num = 0;
      ker = c(0,0,0);
      flag = 0;
      for (k in 1:nrow(triangs)){
        ker = Null(rbind(tentpoles[triangs[k,1],] - c(x[i], y[j]), tentpoles[triangs[k,2],]- c(x[i], y[j]), tentpoles[triangs[k,3],]- c(x[i], y[j])))[,1];
        ker = ker/sum(ker);
        if(ker[1] >=-1e-06 && ker[2] >=-1e-06 && ker[3] >=-1e-06) {
          triang_num = k;
          flag = 1;
          break;
        }
      }
      if (flag == 0) {
        z[i,j] = -Inf;
      } else {
        z[i,j] = ker %*% c(S$par[triangs[triang_num,1]], S$par[triangs[triang_num,2]], S$par[triangs[triang_num,3]]);
      }
    }
  }
  return(list(x = x, y = y, z = z));
}


#####################################################################
# Functions for finding the squared L_2 distance between densities: #
#####################################################################
hellinger <- function(p1, p2, L, R, N) {
  xcoord = seq(L, R, length.out = N);
  ycoord = seq(L, R, length.out = N);
  distance = 0;
  for (i in 1:N) {
    for (j in 1:N) {
      difference = (sqrt(p1(c(xcoord[i], ycoord[j]))) - sqrt(p2(c(xcoord[i], ycoord[j]))));
      distance = distance + difference*difference;
    }
  }
  return(distance*(R-L)*(R-L)/(N*N));
}

squared_L2 <- function(p1, p2, L, R, N) {
  xcoord = seq(L, R, length.out = N);
  ycoord = seq(L, R, length.out = N);
  distance = 0;
  for (i in 1:N) {
    for (j in 1:N) {
      difference = (p1(c(xcoord[i], ycoord[j])) - p2(c(xcoord[i], ycoord[j])));
      distance = distance + difference*difference;
    }
  }
  return(distance*(R-L)*(R-L)/(N*N));
}

standard2dGaussian <- function(x) {
  return(1/(2*pi)*exp(-(x[1]*x[1]+x[2]*x[2])/2));
}

evaluateTentFunction <- function(x, y, z) {
  val = 0;
  z1 = c(z, 1);
  C = upper_chull_triang(x,y);
  for (i in 1:nrow(C)) {
    M = t(x[C[i,],]);
    M = rbind(M, rep(1,ncol(C)));
    coords = solve(M) %*% z1;
    flag = TRUE;
    for (j in 1:ncol(C)) {
      if (coords[j] < 0) {
        flag = FALSE;
        break;
      }
    }
    if (flag == TRUE) {
      return(exp(t(coords[1:ncol(C)]) %*% y[C[i,]]));
    }
  }
  return(0);
}

tentFunction <- function(x, y) {
  s <- function(z) {
    return(evaluateTentFunction(x, y, z));
  }
  return(s);
}
