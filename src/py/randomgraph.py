"""randomgraph.py


"""

from scipy.sparse import coo_matrix
from numba import jit


import numpy as np


class Instance(object):

    def __init__(self, n, a, b, x0, coo, seed):
        self.n = n
        self.a = a
        self.b = b
        self.x0 = x0
        self.coo = coo
        self.adj = coo.tocsr()
        self.seed = seed

    def save(self, filename):
        """Saves the adjacency matrix in an output parsable
        by the Montanari code.

        """

        coo = self.coo
        with open(filename, 'w') as fp:
            for src, dst in zip(coo.row, coo.col):
                if src <= dst:
                    print(src, dst, file=fp)

    @property
    def d(self):
        return (self.a + self.b) / 2

    @property
    def lam(self):
        return (self.a - self.b)/np.sqrt(2.0 * (self.a + self.b))

    def deg(self, node):
        return len(self.adj.indices[self.adj.indptr[node]:self.adj.indptr[node+1]])

    @property
    def nedges(self):
        coo = self.coo
        s = 0
        for src, dst in zip(coo.row, coo.col):
            if src <= dst:
                s += 1
        return s

@jit(nopython=True)
def _generate(x0, a, b, seed):
    np.random.seed(seed)
    n = x0.shape[0]
    inp, outp = a/n, b/n
    rows, cols = [], []
    for src in range(n):
        for dst in range(src, n):
            r = np.random.random()
            x0ij = x0[src] * x0[dst]
            if ((x0ij == +1 and r <= inp) or (x0ij == -1 and r <= outp)):
                rows.append(src)
                cols.append(dst)
                if src != dst:
                    rows.append(dst)
                    cols.append(src)
    return rows, cols


def generate(n, a, b, seed=None):
    """Generate a random undirected graph via the following distribution.
    First, an n length +/- 1 vector called x0 is generated with each entry iid rademacher.

    Then, for i <= j, the edges (i, j) and (j, i) are included with
    probability a/n if x0[i]*x0[j] = 1, and with probability b/n if x0[i]*x0[j] = -1

    """

    if a < 0.0 or a > n:
        raise ValueError("a needs to be in [0, n]")
    if b < 0.0 or b > n:
        raise ValueError("b needs to be in [0, n]")

    if seed is None:
        seed = np.random.randint(0, 0xFFFFFFFF)
    rng = np.random.RandomState(seed)

    x0 = rng.choice([-1, +1], size=n)
    #print(x0)
    #print(x0.dtype)

    rows, cols = _generate(x0, a, b, rng.randint(0, 0xFFFFFFFF))

    coo = coo_matrix(
        (np.ones(len(rows), dtype=np.int8), (np.array(rows), np.array(cols))),
        shape=(n, n))

    return Instance(n, a, b, x0, coo, seed)



if __name__ == '__main__':
    n = 5
    a = 3
    b = 2
    g = generate(n, a, b)
    print (g.adj.todense())
    g.save("test.txt")
    print (g.deg(0), g.deg(1), g.deg(2), g.deg(3), g.deg(4))
