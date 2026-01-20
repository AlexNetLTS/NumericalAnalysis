import numpy as np


def fit(X, y, reg: int):
    x = np.vander(x=X, increasing=True, N=reg+1)
    coef = np.linalg.inv(x.T @ x) @ x.T @ y
    return coef, np.poly1d(coef[::-1])
