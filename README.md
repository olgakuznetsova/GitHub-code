# Log-concave density estimation in undirected graphical models
#### Abstract
We study the problem of maximum likelihood estimation of densities that are log-concave and lie in the graphical model corresponding to a given 
undirected graph G. We show that the maximum likelihood estimate (MLE) is the product of the exponentials of several tent functions, one for each 
maximal clique of G. While the set of log-concave densities in a graphical model is infinite-dimensional, our results imply that the MLE can be 
found by solving a finite-dimensional convex optimization problem. We provide an implementation and a few examples. Furthermore, we show that the 
MLE exists and is unique with probability 1 as long as the number of sample points is larger than the size of the largest clique of G when G is 
chordal. We show that the MLE is consistent when the graph G is a disjoint union of cliques. 
Finally, we discuss the conditions under which a log-concave density in the graphical model of G has a log-concave factorization according to G.

#### Contents of the repository
We include the code for specific examples as well as the underlying functions for computing the following objects in log-concave 
density estimation in undirected graphical models:
1. the MLE 
2. the support of the MLE

#### Reference
Kubjas, K., Kuznetsova, O., Robeva, E., Semnani, P., & Sodomaco, L. (2022). Log-concave density estimation in undirected graphical models. arXiv preprint [arXiv:2206.05227](https://arxiv.org/abs/2206.05227). 
