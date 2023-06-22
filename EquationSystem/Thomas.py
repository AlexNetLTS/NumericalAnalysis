import numpy as np

def thomas(a, b, c, d):
    n = len(d)
    c_ = np.zeros(n - 1)
    d_ = np.zeros(n)
    x = np.zeros(n)
    
    c_[0] = c[0] / b[0]
    d_[0] = d[0] / b[0]
    
    for i in range(1, n-1):
        c_[i] = c[i] / (b[i] - a[i-1] * c_[i-1])
    
    for i in range(1, n):
        d_[i] = (d[i] - a[i-1] * d_[i-1]) / (b[i] - a[i-1] * c_[i-1])
    
    x[n-1] = d_[n-1]
    
    for i in range(n-2, -1, -1):
        x[i] = d_[i] - c_[i] * x[i+1]
    
    return x


if __name__ == '__main__': 
    a = np.array([[4, 1, 0], [1, 4, 1], [0, 1, 4]])  # Coefficient matrix
    b = np.array([1, 2, 3])  # Right-hand side
    x = thomas(a.diagonal(-1), a.diagonal(0), a.diagonal(1), b)
    print("Thomas Algorithm Solution:", x)
    