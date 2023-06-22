import numpy as np


def gaussian_elimination(A, b):
    augmented_matrix = np.concatenate((A, b.reshape(-1, 1)), axis=1)

    rows, cols = augmented_matrix.shape

    for i in range(rows):
        max_row = i
        for j in range(i + 1, rows):
            if abs(augmented_matrix[j, i]) > abs(augmented_matrix[max_row, i]):
                max_row = j

        if max_row != i:
            augmented_matrix[[i, max_row]] = augmented_matrix[[max_row, i]]

        augmented_matrix[i, :] /= augmented_matrix[i, i]

        for j in range(i + 1, rows):
            augmented_matrix[j, :] -= augmented_matrix[j, i] * augmented_matrix[i, :]

    # 回代求解
    x = np.zeros(rows)
    x[rows - 1] = augmented_matrix[rows - 1, cols - 1]
    for i in range(rows - 2, -1, -1):
        x[i] = augmented_matrix[i, cols - 1] - np.dot(augmented_matrix[i, i + 1:cols - 1], x[i + 1:])

    return x
