import numpy as np


def gauss_quadrature(f, a, b, n):
    # 预先计算节点和权重
    x, w = np.polynomial.legendre.leggauss(n)

    # 将积分区间从[-1, 1]映射到[a, b]
    t = 0.5 * (b - a) * x + 0.5 * (b + a)

    # 计算积分结果
    integral = 0.5 * (b - a) * np.sum(w * f(t))

    return integral