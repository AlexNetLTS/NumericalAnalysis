""" 基于大M法的三次样条插值 """


from MSSpline import s
import numpy as np


def base_coef(X, y):
    h = np.diff(X)
    lam = h[1:] / (h[1:] + h[: -1])
    mu = 1 - lam
    delta_y = np.diff(y)
    d = 6 * (delta_y[1:] / h[1:] + delta_y[: -1] / h[: -1]) / (h[1:] + h[: -1])
    return h, delta_y, lam, mu, d


def mspline(xx, X, y, kind=1, df1=None, dfn=None, ddf1=None, ddfn=None):
    """
        Args:
            xx: 代求值
            X: 自变量(是numpy.ndarray)
            y: 因变量(是numpy.ndarray)
            kind: 只取 1 或 2【1是利用两侧一阶导数的信息; 2是利用两侧二阶导数的信息】
            df1: x0处的一阶导数
            dfn: xn处的一阶导数
            ddf1: x0处的二阶导数
            ddfn: xn处的二阶导数
        """

    n = X.shape[0] - 1
    index = np.searchsorted(X, xx, side='left') - 1  # 二分搜索
    h, deltay, lam, mu, d = base_coef(X, y)
    match kind:
        case 1:

            if df1 is None or dfn is None:
                raise ValueError('未设置参数df1, dfn')

            mu = np.concatenate((mu, np.array([1])))
            lam = np.concatenate((np.array([1]), lam))
            mat = 2 * np.eye(n + 1)
            np.fill_diagonal(mat[: -1, 1:], lam)
            np.fill_diagonal(mat[1:, : -1], mu)
            d1, dn = np.array([6 * ((y[1] - y[0])/h[0] - df1) / h[0]]), np.array([6 * (dfn - (y[-1] - y[-2]) / h[-1])])
            d = np.concatenate((d1, d, dn))
            M = np.linalg.solve(mat, d)
            cj1, cj2, cj3 = deltay / h - (2 * M[: -1] + M[1:]) * h / 6, M[: -1] / 2, (-M[: -1] + M[1:]) / h / 6
            return s(xx, index, X, y, cj1, cj2, cj3)

        case 2:

            if ddf1 is None or ddfn is None:
                raise ValueError('未设置参数df1, dfn')

            mat = 2 * np.eye(n - 1)
            np.fill_diagonal(mat[: -1, 1:], lam[: -1])
            np.fill_diagonal(mat[1:, : -1], mu[1:])
            d[0], d[-1] = d[0] - mu[0] * ddf1, d[-1] - lam[-1] * ddfn
            M = np.linalg.solve(mat, d)
            cj1, cj2, cj3 = deltay / h - (2 * M[: -1] + M[1:]) * h / 6, M[: -1] / 2, (-M[: -1] + M[1:]) / h / 6
            return s(xx, index, X, y, cj1, cj2, cj3)