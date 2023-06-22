import numpy as np

''' 基本参量的求解 '''


def base_coef(X, y):
    h = np.diff(X)
    lam = h[1:] / (h[1:] + h[: -1])
    mu = 1 - lam
    delta_y = np.diff(y)
    g = 3 * (mu * delta_y[1:] / h[1:] + lam * delta_y[: -1] / h[: -1])
    return h, delta_y, lam, mu, g


def s(xx, index, X, y, cj1, cj2, cj3):
    return y[index] + cj1[index] * (xx - X[index]) + cj2[index] * (xx - X[index]) ** 2 + cj3[index] * (xx - X[index]) ** 3


def mspline(xx, X, y, kind=1, df1=None, dfn=None, ddf1=None, ddfn=None):
    n = X.shape[0] - 1
    index = np.searchsorted(X, xx, side='left') - 1    # 二分搜索
    h, deltay, lam, mu, g = base_coef(X, y)
    match kind:
        case (1):

            if df1 is None or dfn is None:
                raise ValueError('未设置参数df1, dfn')

            mat = 2 * np.eye(n - 1)
            np.fill_diagonal(mat[: -1, 1:], mu[: -1])  # 修改上对角线
            np.fill_diagonal(mat[1:, : -1], lam[1:])   # 修改下对角线
            m = np.zeros(n + 1)
            m[0], m[-1] = df1, dfn
            g[0], g[-1] = g[0] - lam[0] * df1, g[-1] - mu[-1] * dfn
            m1 = np.linalg.solve(mat, g)               # 此处是可以修改的
            m[1: n] = m1
            cj1, cj2, cj3 = m[: -1], (3 * deltay / h - 2 * m[: -1] - m[1:]) / h, (m[1:] + m[: -1] - 2 * (deltay / h)) / h ** 2
            cj = np.array([cj3, cj2, cj1, y[: -1]])
            return s(xx, index, X, y, cj1, cj2, cj3), cj

        case (2):

            if ddf1 is None or ddfn is None:
                raise ValueError('未设置参数ddf1, ddfn')

            mu, lam = np.concatenate((np.array([1]), mu)), np.concatenate((lam, np.array([1])))
            mat = 2 * np.eye(n + 1)
            np.fill_diagonal(mat[: -1, 1:], mu)        # 修改上对角线
            np.fill_diagonal(mat[1:, : -1], lam)       # 修改下对角线
            g = np.concatenate((np.array([3 * (deltay[0] / h[0] + 0.5 * h[0] * ddf1)]), g,
                                np.array([3 * (deltay[-1] / h[-1] + 0.5 * h[-1] * ddfn)])))
            m = np.linalg.solve(mat, g)                # 此处是可以修改的

            cj1, cj2, cj3 = m[: -1], (3 * deltay / h - 2 * m[: -1] - m[1:]) / h, (m[1:] + m[: -1] - 2 * (deltay / h)) / h ** 2
            cj = np.array([cj3, cj2, cj1, y[: -1]])
            return s(xx, index, X, y, cj1, cj2, cj3), cj
