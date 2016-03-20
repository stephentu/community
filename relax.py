"""relax.py


Implements the block coordinate algorithm described in Section 10.1.2 of
http://web.stanford.edu/~montanar/SDPgraph/sdp_phtr.pdf


"""

import numpy as np


def solve(graph, r, eta=None, epochs=1, seed=None, conv_tol=1e-4):
    """

    """

    if r <= 0 or r > graph.n:
        raise ValueError("rank should never exceed n")

    if seed is None:
        seed = np.random.randint(0, 0xFFFFFFFF)
    rng = np.random.RandomState(seed)

    # Compute eta if necessary
    if eta is None:
        eta = 2.0 * graph.nedges / (graph.n * graph.n)

    print ("n", graph.n, "r", r, "eta", eta, "epochs", epochs, "conv_tol", conv_tol)

    # Initialize the spins
    spins = np.zeros((n, r))
    for node in range(graph.n):
        if graph.deg(node):
            spins[node] = rng.choice([-1, +1], size=r)
            spins[node] /= np.linalg.norm(spins[node])

    indices, indptr = graph.adj.indices, graph.adj.indptr
    M = spins.sum(axis=0) # \sum_i^n spin_i
    #assert M.shape == (r,)

    for epoch in range(epochs):
        maxdiff = 0.0
        for node in rng.permutation(graph.n):
            #assert np.allclose(M, spins.sum(axis=0))

            # sum the spins over the neighbors of node
            #print (node, "node", "indices", indices[indptr[node]:indptr[node+1]])
            cur = -eta * M + spins[indices[indptr[node]:indptr[node+1]]].sum(axis=0)
            curnorm = np.linalg.norm(cur)
            if np.fabs(curnorm) <= 1e-10:
                # skip
                continue
            cur /= curnorm
            diff = cur - spins[node]
            M += 2.0 * diff
            maxdiff = max(maxdiff, np.dot(diff, diff))
            spins[node] = cur

        print ("epoch", epoch, "maxdiff", maxdiff)

        if maxdiff <= conv_tol:
            break

    #print (spins)

    # now compute the rounding
    print (spins)
    cov = 1.0/graph.n * spins.T.dot(spins) # empirical covariance

    evals, evecs = np.linalg.eigh(cov)
    print ("evals", evals)
    return np.sign(spins.dot(evecs[:, np.argmax(evals)]))


if __name__ == '__main__':
    from randomgraph import generate

    n = 1000
    a = 5
    b = 3
    g = generate(n, a, b)
    #print (g.adj.todense())
    g.save("test.txt")

    print ("d", g.d, "lam", g.lam)
    print (solve(g, r=20, epochs=100))
