# Newton插值法


import numpy as np


def diff_coef(X, Y):
    dim = X.shape[0]
    primary = []
    for i in range(dim):
        if i == 0:
            num = 1 / np.prod(X[i] - X[i + 1:])
            primary.append(num)
        elif i == dim - 1:
            num = 1 / np.prod(X[i] - X[: i])
            primary.append(num)
        else:
            num = 1 / (np.prod(X[i] - X[: i]) * np.prod(X[i] - X[i + 1: ]))
            primary.append(num)
    return np.dot(np.array(primary), Y.T)


def newton(xx, X, Y):
    dim = X.shape[0]
    primary, prot = [], []
    for i in range(1, dim):
        num = diff_coef(X=X[: i + 1], Y=Y[: i + 1])
        prot.append(np.prod(xx - X[: i]))
        primary.append(num)
    # print(prot)
    # print(primary)
    return np.dot(np.array(prot), np.array(primary)) + Y[0]