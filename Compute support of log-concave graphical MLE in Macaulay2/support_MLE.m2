restart
load "support_MLE_functions.m2"
-- input vertices of the graph
d = 4
V = toList(1..d)
-- input edges of the graph 
E = {}
E = {{1,2},{2,3},{3,4},{4,1}}
E = {{1,2},{2,3},{3,4},{4,5},{5,1}}
E = {{1,2},{2,3},{3,4},{4,5},{5,6},{6,1}}
E = {{1,2},{2,3},{3,1}}
E = {{1,2},{2,3},{3,4}}
E = {{1,2},{2,3},{3,1},{2,4},{2,5}}
-- compute maximal cliques and global Markov statements that form a partition of the vertex set
(C,globalG) = cliquesGlobalMarkov(V,E)

-- consider a random sample of n integer points in RR^d
-- (points are given by columns of V)
n = 2
KK = QQ
KK = ZZ
M = random(KK^d,KK^n)

-- specific examples that converge for the 4-cycle:
M = matrix {{7, 3, 7, 8}, {2, 7, 9, 0}, {8, 9, 8, 1}, {0, 3, 4, 8}}
M = transpose matrix {{1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 1},{0,0,0,0}}
M = transpose matrix {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}
M = transpose matrix {{1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 1}}
M = transpose matrix {{1, 0, 0, 0, 0}, {0, 1, 0, 0, 0}, {0, 0, 1, 0, 0}, {0, 0, 0, 1, 0}, {0, 0, 0, 0, 1}}
M = transpose matrix {{1, 0, 0, 0, 0, 0}, {0, 1, 0, 0, 0, 0}, {0, 0, 1, 0, 0, 0}, {0, 0, 0, 1, 0, 0}, {0, 0, 0, 0, 1, 0}, {0, 0, 0, 0, 0, 1}}
M= transpose matrix{{0,0,0,0},{1,0,0,0},{1,1,0,0},{1,1,1,0},{0,0,0,1},{0,0,1,1},{0,1,1,1},{1,1,1,1}}
-- compute the support of the MLE
time supp = suppMLE(M,C);
-- compute the volume of the support of the MLE
(volume supp)_RR

-- compute recursively the sets D_G^i
-- the command prints whether the algorithm converges or not after N iterations
N = 1;
elapsedTime (O,I,D) = DGrecursive(N,M,globalG);
D == supp -- if the algorithm converged, check that the limit coincides with the support of MLE
