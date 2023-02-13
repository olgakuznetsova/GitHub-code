# Log-concave density estimation in undirected graphical models
#### Abstract
We study the problem of maximum likelihood estimation of densities that are log-concave and lie in the graphical model corresponding to a given 
undirected graph $G$. We show that the maximum likelihood estimate (MLE) is the product of the exponentials of several tent functions, one for each 
maximal clique of $G$. While the set of log-concave densities in a graphical model is infinite-dimensional, our results imply that the MLE can be 
found by solving a finite-dimensional convex optimization problem. We provide an implementation and a few examples. Furthermore, we show that the 
MLE exists and is unique with probability 1 as long as the number of sample points is larger than the size of the largest clique of $G$ when $G$ is 
chordal. We show that the MLE is consistent when the graph G is a disjoint union of cliques. 
Finally, we discuss the conditions under which a log-concave density in the graphical model of $G$ has a log-concave factorization according to $G$.

#### Contents of the repository
The repository contains two folders

1. "Compute log-concave graphical MLE in R", which contains the implementation of the optmization algorithm presented in Section 6 of the corresponding paper. The necessary functions are contained in the file "functions.R", which are then called by "log-concave_graphical_MLE.R". One needs to add the path to "functions.R" on line 18 of "log-concave_graphical_MLE.R". To facilitate replication, we also include the data that was used in Example 6.1 (folder "Graph_12") and Example 6.2 ("Graph_12-23"). The computations are performed in [R](https://www.r-project.org).  
2. "Compute support of log-concave graphical MLE in Macaulay2", which contain the code for studying whether the special sets $D_G^{(i)}$ converge to the support of the MLE for non-chordal graphs as discussed in Conjecture 2.11. The file "support_MLE_functions.m2" includes the declarations of all relevant functions that are called inside "support_MLE.m2" for the computation of the support of a 4-cycle. Both files should be placed in the same folder.  The computations are performed in the [Macaulay2](https://www.r-project.org).   

#### Reference
Kubjas, K., Kuznetsova, O., Robeva, E., Semnani, P., & Sodomaco, L. (2022). Log-concave density estimation in undirected graphical models. arXiv preprint [arXiv:2206.05227](https://arxiv.org/abs/2206.05227). 
