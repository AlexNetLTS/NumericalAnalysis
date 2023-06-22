# Lagrange插值法


import numpy as np


def lagrange(xx, X, Y): 
    dim = X.shape[0]
    primary = []
    for i in range(dim):
        if i == 0:
            num = np.prod(xx - X[i + 1: ]) / np.prod(X[i] - X[i + 1: ])
            primary.append(num)
        elif i == dim - 1:
            num = np.prod(xx - X[: i]) / np.prod(X[i] - X[: i])
            primary.append(num)
        else:
            num = (np.prod(xx - X[: i]) * np.prod(xx - X[i + 1: ])) / (np.prod(X[i] - X[: i]) * np.prod(X[i] - X[i + 1: ]))
            primary.append(num)
    return np.dot(np.array(primary), Y.T)
