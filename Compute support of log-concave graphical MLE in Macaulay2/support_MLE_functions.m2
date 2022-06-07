needsPackage "GraphicalModels";
needsPackage "Polyhedra";

-- Function "cliquesAndGlobalMarkov":
-- INPUT:
-- 1) vertSet: vertex set, a list of vertices
-- 2) edgeSet: edge set, a list of lists of elements in vertSet
-- OUTPUT:
-- 1) the list of maximal cliques of the graph G = (vertSet, edgeSet)
-- 2) the list of global Markov statements of G that form a partition of vertSet
cliquesGlobalMarkov = (vertSet,edgeSet) -> (
    G := graph(vertSet,edgeSet);
    numVertices := #(vertSet);
    cl := facets cliqueComplex(G);
    C := apply(cl, c-> indices(c)+toList((degree(c))#0:1));
    gMG := globalMarkov(G);
    gMGpartition := {};
    for g in gMG do if set((g#0)|(g#1)|(g#2)) === set(toList(1..numVertices)) then gMGpartition = append(gMGpartition, g);
    return (C,gMGpartition)
    );

-- Function "suppMLE":
-- INPUT:
-- 1) X: d x n matrix containing the vertices of the sample
-- 2) cliques: maximal cliques of a graph G
-- OUTPUT:
-- 1) a polytope, the support of the log-concave graphical MLE
suppMLE = (X,cliques) -> (
    prisms := {};
    dimSpace = numRows X;
    for c in cliques do (
	projcX = convexHull(X^(sort apply(c,i-> i-1)));
	F = facets projcX;
	nFacets = numRows F_0;
	Fnew0 = matrix for i in 0..nFacets-1 list for j in 1..dimSpace list if member(j,c) then (F_0)_(i,position(c, k-> k==j)) else 0;
	prisms = append(prisms,polyhedronFromHData((Fnew0,F_1)))
	);
    return intersection(prisms)
    );

-- Function "addSwapped":
-- INPUT:
-- 1) X: d x n matrix
-- 2) A,B,S: lists that form a partition (A,B,S) of the set {1,...,d}
-- OUTPUT:
-- 1) a matrix with d rows whose columns are the vertices of the polytope Q_{A,B,S}(P),
--    where P is the convex hull of the sample given by the columns of X and Q_{A,B,S} is the function defined in Lemma ??
addSwapped = (X,A,B,S) -> (
    dimSpace := numRows X;
    AS := apply(sort(A|S), i-> i-1);
    BS := apply(sort(B|S), i-> i-1);
    projAS := convexHull(X^AS);
    projBS := convexHull(X^BS);
    FAS := facets projAS;
    FBS := facets projBS;
    nFacetsAS := numRows FAS_0;
    nFacetsBS := numRows FBS_0;
    FASnew0 := matrix for i in 0..nFacetsAS-1 list for j in 0..dimSpace-1 list if member(j,AS) then (FAS_0)_(i,position(AS, k-> k==j)) else 0;
    FBSnew0 := matrix for i in 0..nFacetsBS-1 list for j in 0..dimSpace-1 list if member(j,BS) then (FBS_0)_(i,position(BS, k-> k==j)) else 0;
    prismAS := polyhedronFromHData((FASnew0,FAS_1));
    prismBS := polyhedronFromHData((FBSnew0,FBS_1));
    return vertices intersection(prismAS,prismBS)
    );

-- Function "mapDG":
-- INPUT:
-- 1) X: d x n matrix
-- 2) CI: list of lists of elements of {1,...,d}, corresponding to the list of conditional independence statements "A independent B given S" of a graph G,
--    where (A,B,S) is a partition of the set {1,...,d}.
-- OUTPUT:
-- 1) the polytope D_G(X), where D_G is the function defined in Lemma ??
DG = (X,CI) -> (if #CI==0 then return convexHull X else return convexHull apply(CI, c-> addSwapped(X,c_0,c_1,c_2)));

-- Function "DGrecursive":
-- INPUT:
-- 1) numIterations: a nonnegative integer
-- 2) X: d x n matrix
-- 3) CI: list of lists of elements of {1,...,d}, corresponding to the list of conditional independence statements "A independent B given S" of a graph G,
--    where (A,B,S) is a partition of the set {1,...,d}.
-- OUTPUT:
-- 1) the polytope D_G^i(X), where i = numIterations and D_G^i is the function defined in Lemma ?? composed i times
DGrecursive = (numIterations,X,CI) -> (
    ind := 0;
    Xold := X;
    Pold := convexHull(X);
    volold := (volume(Pold))_RR;
    dimensionold := dim(Pold);
    print concatenate{"-- index i = ",toString(ind),", dimension = ",toString(dimensionold),", relative volume = ",toString(volold)};
    s := false;
    out := {};
    while ind < numIterations+1 and s==false do (
	ind = ind + 1;
	if ind > 1 then Xold = Xnew;
	if ind > 1 then Pold = Pnew;
	Pnew = DG(Xold,globalG);
	Xnew = vertices(Pnew);
	if Pnew == Pold then s = true;
	volnew = (volume(Pnew))_RR;
	dimensionnew = dim(Pnew);
	out = append(out, (ind,dimensionnew,volnew));
	if s==false then print concatenate{"-- index i = ",toString(ind),", dimension = ",toString(dimensionnew),", number of vertices = ",toString(numColumns(Xnew)),", relative volume = ",toString(volnew)};
	);
    if s==true then print concatenate{"-- converged at index i = ",toString(ind-1)};
    if s==false then print concatenate{"-- not converged after ",toString(numIterations)," iterations"};
    return (out,ind,Pnew)
    );
