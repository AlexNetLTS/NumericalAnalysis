import numpy as np


class Interpolation(object):
    def __init__(self, X, Y):
        """
            X: 自变量
            Y: 因变量 
            coef: 相应的系数
        """
        self.X = X
        self.Y = Y
        self.coef = []

    def lagrange(self, xx):
        self.coef.clear()
        dim = self.X.shape[0]
        for i in range(dim):
            if i == 0:
                num = np.prod(xx - self.X[i + 1: ]) / np.prod(self.X[i] - self.X[i + 1: ])
                self.coef.append(num)
            elif i == dim - 1:
                num = np.prod(xx - self.X[: i]) / np.prod(self.X[i] - self.X[: i])
                self.coef.append(num)
            else:
                num = (np.prod(xx - self.X[: i]) * np.prod(xx - self.X[i + 1: ])) / (np.prod(self.X[i] - self.X[: i]) * np.prod(self.X[i] - self.X[i + 1: ]))
                self.coef.append(num)
        return np.dot(np.array(self.coef), self.Y.T)
    
    def diff_coef(self, X, Y):
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


    def newton(self, xx):
        dim = self.X.shape[0]
        primary, prot = [], []
        for i in range(1, dim):
            num = self.diff_coef(self.X[: i + 1], self.Y[: i + 1])
            prot.append(np.prod(xx - self.X[: i]))
            primary.append(num)
        return np.dot(np.array(prot), np.array(primary)) + self.Y[0]