#include <iostream>
using namespace std;

#define MAX_SIZE 100

void gaussian_elimination(double A[MAX_SIZE][MAX_SIZE], double b[MAX_SIZE], int n) {
    int i, j, k;
    
    // 高斯消元
    for (i = 0; i < n - 1; i++) {
        // 如果主对角线元素为0，则进行行交换以避免除以零的错误
        if (A[i][i] == 0) {
            for (j = i + 1; j < n; j++) {
                if (A[j][i] != 0) {
                    // 交换第i行和第j行
                    for (k = 0; k < n; k++) {
                        double temp = A[i][k];
                        A[i][k] = A[j][k];
                        A[j][k] = temp;
                    }
                    double temp = b[i];
                    b[i] = b[j];
                    b[j] = temp;
                    break;
                }
            }
        }
        
        // 将主对角线元素缩放为1
        double pivot = A[i][i];
        for (j = i; j < n; j++) {
            A[i][j] /= pivot;
        }
        b[i] /= pivot;
        
        // 将主对角线以下的元素消为零
        for (j = i + 1; j < n; j++) {
            double factor = A[j][i];
            for (k = i; k < n; k++) {
                A[j][k] -= factor * A[i][k];
            }
            b[j] -= factor * b[i];
        }
    }
    
    // 回代求解
    double x[MAX_SIZE];
    x[n - 1] = b[n - 1];
    for (i = n - 2; i >= 0; i--) {
        double sum = 0;
        for (j = i + 1; j < n; j++) {
            sum += A[i][j] * x[j];
        }
        x[i] = b[i] - sum;
    }
    
    // 输出解向量
    printf("Solution vector:\n");
    for (i = 0; i < n; i++) {
        printf("x[%d] = %f\n", i, x[i]);
    }
}

int main() {
    int n;
    printf("Enter the size of the system: ");
    scanf("%d", &n);
    
    double A[MAX_SIZE][MAX_SIZE];
    double b[MAX_SIZE];
    
    printf("Enter the coefficients of the system:\n");
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            scanf("%lf", &A[i][j]);
        }
        scanf("%lf", &b[i]);
    }
    
    gaussian_elimination(A, b, n);
    
    return 0;
}
