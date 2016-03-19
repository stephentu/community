The C code SDPclustering.c implements the SDP-based clustering algorithm developed in “Phase Transitions in Semidefinite Relaxations”, Adel Javanmard, Andrea Montanari and Federico Ricci-Tersenghi, arxiv:1511.08769

CODE COMPILING
The C code is self-contained and only need standard libraries (stdio.h, stdlib.h, ctype.h and math.h). It can be compiled with a standard simple command like

gcc -O3 SDPclustering.c —o SDPclustering -lm

The optimization flag -O3 is not strictly required, but strongly suggested, since it largely improves code performances. The executable file name (SDPclustering in the example above) can be set to any desired value.

INPUT
At execution time, the code requires some input parameter at the command line.
The usage is

./SDPclustering <graphFile> <m> <maxIter> <randomSeed> [<gamma>]

In <graphFile> a file should be provided containing the information on the edge set of undirected and unweighted graph (further details on the graph file format are given below).

<m> : The number of components of the vertex variables. Equivalently, the algorithm optimizes over
positive-semidefinite matrices of rank m.

<maxIter>: Maximum number of iterations (each iterations corresponds on average to one update per variable)

<randomSeed>: Seed for  the pseudo-random number generator is required.

<gamma>: Regularization parameter (optional). Large gamma enforces more balanced solutions.  If not provided it is set by default  to d/n, where d is the graph average degree and n the number of nodes.

The code uses 2 tolerance levels (for checking convergence to the maximum of the objective function and to compute the eigensystem of the empirical covariance matrix) which are set via #define directives in the first lines of the code.

GRAPH FILE FORMAT
The input is an undirected unweighted graph, which is given in a file with the following format.
Each line contain an edge, specified as a pair of nodes.
Furthermore:
- Any line beginning with a # is considered as a comment and ignored;
- Any line with at least 2 numbers is considered as an undirected  edge connecting the nodes identified by the integer part of the first 2 numbers in the line;
- Any line with less than 2 numbers returns an error and stops the execution.

OUTPUT
The program writes only on stdout.

OUTPUT HEADER AND INFO
The header of the output contains a summary of the input and of of parameters used.
During the coordinate ascent loop, at each time step, the maximum variation in the variables is provided (as a check for convergence).

Once coordinate ascent converges within the assigned tolerance, the program computes the empirical covariance matrix
(an m x m matrix) at the maximizer, and the corresponding eigenvectors and eigenvalues.
Eigenvectors corresponding to eigenvalues smaller than than 1e-6 are discarded.
Eigenvalues larger than this threshold and their number (denoted below by k)  are reported.

OUTPUT EMBEDDING:
The output of the optimization algorithm is printed as a list of tuples

vertex  degree x1 x2 ... xk

One such line is present per each vertex. Vertices of degree 0 are also reported, but the corresponding
row only contains the pair vertex-degree. For the other rows, (x1,...,xk) provides a k-dimensional embedding
of the vertex. Coordinates correspond, in order to the leading eigenvectors.
Lower-dimensional embeddings are obtained by keeping the first k'<k coordinates. 

For updates and contact information, see
https://web.stanford.edu/~montanar/sdp_clustering/
