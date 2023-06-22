import numpy as np


def romberg_integration(f, a, b, max_iters, tolerance):
    """
    Args:
        f: 积分的函数
        a: 积分区间上限
        b: 积分区间下限
        max_iters: 最大迭代次数
        tolerance: 容差
    """
    R = np.zeros((max_iters, max_iters), dtype=float)
    h = b - a

    R[0, 0] = 0.5 * h * (f(a) + f(b))
    iter_count = 0

    while iter_count < max_iters - 1:
        iter_count += 1
        h *= 0.5
        sum_term = 0

        # 计算新的R[i, 0]
        for k in range(1, 2**iter_count, 2):
            sum_term += f(a + k * h)

        R[iter_count, 0] = 0.5 * (R[iter_count-1, 0] + h * sum_term)

        # 计算其他列的值
        for j in range(1, iter_count + 1):
            R[iter_count, j] = (4**j * R[iter_count, j-1] - R[iter_count-1, j-1]) / (4**j - 1)

        # 检查收敛条件
        if np.abs(R[iter_count, iter_count] - R[iter_count-1, iter_count-1]) < tolerance:
            break

    return R[iter_count, iter_count]

# 定义要积分的函数
def f(x):
    return np.sin(x)

# 定义积分区间
a = 0
b = np.pi

# 定义最大迭代次数和容差
max_iters = 10
tolerance = 1e-6

# 使用龙贝格算法计算积分结果
result = romberg_integration(f, a, b, max_iters, tolerance)

# 打印结果
print("积分结果:", result)
