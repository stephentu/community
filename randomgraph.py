"""randomgraph.py


"""

from scipy.sparse import coo_matrix


import numpy as np


class Instance(object):

    def __init__(self, n, a, b, x0, adj, seed):
        self.n = n
        self.a = a
        self.b = b
        self.x0 = x0
        self.adj = adj
        self.seed = seed


    def save(self, filename):
        """Saves the adjacency matrix in an output parsable
        by the Montanari code.

        """

        coo = self.adj.tocoo()
        with open(filename, 'w') as fp:
            for src, dst in zip(coo.row, coo.col):
                if src <= dst:
                    print(src, dst, file=fp)


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
    print(x0)
    print(x0.dtype)

    inp, outp = a/n, b/n

    rows, cols = [], []
    for src in range(n):
        for dst in range(src, n):
            r = rng.uniform()
            if ((x0[src] * x0[dst] == +1 and r <= inp) or
                (x0[src] * x0[dst] == -1 and r <= outp)):
                rows.append(src)
                cols.append(dst)
                if src != dst:
                    rows.append(dst)
                    cols.append(src)

    adj = coo_matrix(
        (np.ones(len(rows), dtype=np.int8), (np.array(rows), np.array(cols))),
        shape=(n, n)).tocsr()

    return Instance(n, a, b, x0, adj, seed)



if __name__ == '__main__':
    n = 5
    a = 3
    b = 2
    g = generate(n, a, b)
    print (g.adj.todense())
    g.save("test.txt")
