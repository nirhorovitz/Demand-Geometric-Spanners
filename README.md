# Demand-Geometric-Spanners

In this project, we study the problem of constructing sparse geometric networks with weighted stretch. 
Given a set $P$ of $n$ points in the Euclidean plane, a stretch parameter t > 1, 
and an edge-demand function
w : P * P -> (0,1]
our goal is to compute a graph G = (P, E) satisfying the following properties:

1. ### Edge-demand t-spanner.
   For every pair of points p, q ∈ P , the graph contains a path whose
   Euclidean length is at most t · |pq|/w(p, q), where |pq| denotes the Euclidean distance between p and q.
2. ### Sparsity.
   The number of edges is linear: |E| = O(n), where the hidden constant depends only on t and the demand function w.
3. ### Lightness.
   The total weight of the graph is O(wt(M ST (P ))), where the hidden constant again depends on t and on the demand function w.

We study constructions that satisfy these constraints simultaneously and analyze their depen-
dence on the stretch parameter t and the edge-demand function w. Our goal is to compute a graph
whose size and weight are within a constant factor of those of the optimal graph satisfying these
constraints.
