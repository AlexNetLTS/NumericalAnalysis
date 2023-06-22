import numpy as np


def lu(A):
    n = A.shape[0]
    L = np.eye(n)
    U = np.zeros_like(A)

    for k in range(n):
        U[k, k:] = A[k, k:] - L[k, :k] @ U[:k, k:]
        L[(k+1):, k] = (A[(k+1):, k] - L[(k+1):, :k] @ U[:k, k]) / U[k, k]

    return L, U



def lu_solve(L, U, b):
    n = L.shape[0]
    y = np.zeros(n)
    x = np.zeros(n)

    # Ly = b
    for i in range(n):
        y[i] = b[i] - L[i, :i] @ y[:i]

    # Ux = y
    for i in range(n - 1, -1, -1):
        x[i] = (y[i] - U[i, i + 1:] @ x[i + 1:]) / U[i, i]

    return x
