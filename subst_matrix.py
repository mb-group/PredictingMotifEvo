__author__ = "Greg Slodkowicz <gslodko@mrc-lmb.cam.ac.uk>"

import numpy as np

nts = ['U', 'C', 'A', 'G']

def eigen(Q):
    evals, evecs = np.linalg.eig(Q)
    order = np.argsort(evals)
    evals = evals[order]
    evecs = evecs[:, order]
    ivecs = np.linalg.inv(evecs)

    return (np.ascontiguousarray(evecs),
            np.ascontiguousarray(evals),
            np.asfortranarray(ivecs))


class NamedMat(object):
    def __init__(self, M, rownames, colnames):
        self.M = M
        
        assert len(rownames) == M.shape[0]
        assert len(colnames) == M.shape[1]

        self.rownames = dict(zip(rownames, range(M.shape[0])))
        self.colnames = dict(zip(colnames, range(M.shape[1])))

    def __getitem__(self, i):
        return self.M[self.rownames[i[1]],self.colnames[i[0]]]


class SubstMatrix(object):
    def __init__(self, nstates):
        self.Q = np.zeros((nstates, nstates))

    def initQ(self, M, freqs):
        self.Q = M.dot(np.diag(freqs))

    def getP(self, t):
        U, l, Ut = eigen(self.Q)

        l = np.exp(l*t)

        P = U.dot(np.diag(l)).dot(Ut) 
        P = NamedMat(P, nts, nts)

        return P


class GTRSubstMatrix(SubstMatrix):
    """
    Assuming the nt order of T, C, A, G 
    """
    def __init__(self, a, b, c, d, e, f, freqs):
        super(GTRSubstMatrix, self).__init__(4)

        M = np.ascontiguousarray([[ -(a+b+c), a, b, c ],
                                  [ a, -(a+d+e), d, e ],
                                  [ b, d, -(b+d+f), f ],
                                  [ c, e, f, -(c+e+f) ]])

        self.freqs = freqs

        self.initQ(M, freqs)


# print substmat.getP(0.1)
# print substmat.getP(100)

# substmat = GTRSubstMatrix(a=1,
#                           b=1,
#                           c=1,
#                           d=1,
#                           e=1,
#                           f=1,
#                           freqs=[0.25, 0.25, 0.25, 0.25])
# print substmat.getP(0.05)
# print substmat.getP(0.1)
# print substmat.getP(100)
