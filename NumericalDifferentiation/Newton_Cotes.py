import numpy as np


def trapezoidal_rule(f, a, b, n):
    h = (b - a) / n
    x = np.linspace(a, b, n + 1)
    y = f(x)
    area = h * (np.sum(y) - 0.5 * (y[0] + y[-1]))
    return area


def simpson_rule(f, a, b, n):
    if n % 2 != 0:
        raise ValueError("n must be an even number for Simpson's rule.")
    
    h = (b - a) / n
    x = np.linspace(a, b, n+1)
    y = f(x)
    area = h / 3 * (np.sum(y[::2]) + 4 * np.sum(y[1::2]) + y[-1])
    return area


def adaptive_composite_trapezoidal(f, a, b, tol):
    # 自适应变步长复化的梯形公式
    n = 1  # 初始分段数
    integral_old = trapezoidal_rule(f, a, b, n)  # 初始积分值
    while True:
        n *= 2  # 增加分段数
        integral_new = trapezoidal_rule(f, a, b, n)  # 新的积分值
        if abs(integral_new - integral_old) < tol:
            return integral_new
        integral_old = integral_new