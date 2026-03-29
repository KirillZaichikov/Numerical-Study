#include <cmath>
#include <cstring>
#include <iostream>

const int N = 6;

void multiply6x6(double A[N][N], double B[N][N], double C[N][N]) {
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            C[i][j] = 0.0;
            for (int k = 0; k < N; ++k) {
                C[i][j] += A[i][k] * B[k][j];
            }
        }
    }
}

void multiply3x3(double A[3][3], double B[3][3], double C[3][3]) {
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            C[i][j] = 0.0;
            for (int k = 0; k < 3; ++k) {
                C[i][j] += A[i][k] * B[k][j];
            }
        }
    }
}


int invert_matrix(const double A[6][6], double invA[6][6]) {
    double aug[6][2 * 6]; // Расширенная матрица [A | I]
    int i, j, k;

    // Формируем расширенную матрицу [A | I]
    for (i = 0; i < 6; i++) {
        for (j = 0; j < 6; j++) {
            aug[i][j] = A[i][j];
            aug[i][j + 6] = (i == j) ? 1.0 : 0.0; // единичная матрица справа
        }
    }

    // Прямой ход метода Гаусса–Жордана
    for (i = 0; i < 6; i++) {
        // Находим ведущий элемент
        double pivot = aug[i][i];
        if (fabs(pivot) < 1e-12) {
            // ищем строку для обмена
            int swap_row = -1;
            for (int r = i + 1; r < 6; r++) {
                if (fabs(aug[r][i]) > 1e-12) {
                    swap_row = r;
                    break;
                }
            }
            if (swap_row == -1)
                return 0; // вырожденная матрица — обратной нет

            // обмен строк
            for (j = 0; j < 2 * 6; j++) {
                double tmp = aug[i][j];
                aug[i][j] = aug[swap_row][j];
                aug[swap_row][j] = tmp;
            }
            pivot = aug[i][i];
        }

        // Нормализуем строку (делим на ведущий элемент)
        for (j = 0; j < 2 * 6; j++)
            aug[i][j] /= pivot;

        // Обнуляем столбцы выше и ниже ведущего элемента
        for (k = 0; k < 6; k++) {
            if (k == i) continue;
            double factor = aug[k][i];
            for (j = 0; j < 2 * 6; j++)
                aug[k][j] -= factor * aug[i][j];
        }
    }

    // Копируем правую часть (обратную матрицу)
    for (i = 0; i < 6; i++)
        for (j = 0; j < 6; j++)
            invA[i][j] = aug[i][j + 6];

    return 1;
}

int inverse3x3(double A[3][3], double inv[3][3]) {
    double det =
          A[0][0]*(A[1][1]*A[2][2] - A[1][2]*A[2][1])
        - A[0][1]*(A[1][0]*A[2][2] - A[1][2]*A[2][0])
        + A[0][2]*(A[1][0]*A[2][1] - A[1][1]*A[2][0]);

    if (det == 0) {
        return 0; // матрица необратима
    }

    inv[0][0] =  (A[1][1]*A[2][2] - A[1][2]*A[2][1]) / det;
    inv[0][1] = -(A[0][1]*A[2][2] - A[0][2]*A[2][1]) / det;
    inv[0][2] =  (A[0][1]*A[1][2] - A[0][2]*A[1][1]) / det;

    inv[1][0] = -(A[1][0]*A[2][2] - A[1][2]*A[2][0]) / det;
    inv[1][1] =  (A[0][0]*A[2][2] - A[0][2]*A[2][0]) / det;
    inv[1][2] = -(A[0][0]*A[1][2] - A[0][2]*A[1][0]) / det;

    inv[2][0] =  (A[1][0]*A[2][1] - A[1][1]*A[2][0]) / det;
    inv[2][1] = -(A[0][0]*A[2][1] - A[0][1]*A[2][0]) / det;
    inv[2][2] =  (A[0][0]*A[1][1] - A[0][1]*A[1][0]) / det;

    return 1;
}

double dotProduct(const double* firstVec, const double* secondVec, const int length) {
    double res = 0;
    for (int i = 0; i < length; i++)
        res += firstVec[i] * secondVec[i];
    return res;
}

double vecNorm(const double* vec, const int length) {
    double norm = 0;
    for (int i = 0; i < length; i++)
        norm += vec[i] * vec[i];
    return sqrt(norm);
}

void normalizeVecs(double* vecs, int dimension, int numVec, double eps) {
    double norm = 0;
    for (int i = 0; i < numVec; i++) {
        norm = vecNorm(&vecs[i * dimension], dimension);
        for (int j = 0; j < dimension; j++)
            vecs[i * dimension + j] = (vecs[i * dimension + j] / norm) * eps;
    }
}

double getAngle(const double* vec1, const double* vec2, const double dimension) {
    return acos((dotProduct(vec1, vec2, dimension)) / (vecNorm(vec1, dimension) * vecNorm(vec2, dimension)));
}

void ortVecs(double* vecs, int dimension, int numVecs, double* projSum) {
    double projCoeff = 0;

    for (int i = 1; i < numVecs; i++) {
        for (int j = 0; j < i; j++) {
            projCoeff = dotProduct(&vecs[j * dimension], &vecs[i * dimension], dimension) / dotProduct(&vecs[j * dimension], &vecs[j * dimension], dimension);
            for (int k = 0; k < dimension; k++)
                projSum[k] += projCoeff * vecs[j * dimension + k];
        }
        for (int j = 0; j < dimension; j++) {
            vecs[i * dimension + j] -= projSum[j];
            projSum[j] = 0;
        }
    }
}

void matrixTranspose(double* matr, int rows, int cols) {
    double temp;

    for (int i = 0; i < rows; i++)
        for (int j = 0; j < cols; j++) {
            if (j > i) {
                temp = matr[i * cols + j];
                matr[i * cols + j] = matr[j * cols + i];
                matr[j * cols + i] = temp;
            }
        }
}

void matrixMult(double* res, double* A, int rowsA, int colsA, double* B, int rowsB, int colsB) {
    for (int i = 0; i < rowsA; i++)
        for (int j = 0; j < colsB; j++) {
            res[i * colsB + j] = 0;
            for (int k = 0; k < colsA; k++)
                res[i * colsB + j] += (A[i * colsA + k] * B[k * colsB + j]);
        }
}

double dotProductColMaj(double* firstVec, double* secondVec, int rows, int cols) {
    double res = 0;
    for (int i = 0; i < rows; i++)
        res += firstVec[i * cols] * secondVec[i * cols];
    return res;
}

double vecNormColMaj(double* vec, int rows, int cols) {
    double norm = 0;
    for (int i = 0; i < rows; i++)
        norm += vec[i * cols] * vec[i * cols];
    return sqrt(norm);
}

void normalizeVecsColMaj(double* vecs, int rows, int cols, double eps) {
    double norm = 0;
    for (int i = 0; i < rows; i++) {
        norm = vecNormColMaj(&vecs[i], rows, cols);
        for (int j = 0; j < cols; j++)
            vecs[i + j * cols] = (vecs[i + j * cols] / norm) * eps;
    }
}

void ortVecsColMaj(double* vecs, int rows, int cols, double* projSum) {
    double projCoeff = 0;

    for (int i = 1; i < cols; i++) {
        for (int j = 0; j < i; j++) {
            projCoeff = dotProductColMaj(&vecs[j], &vecs[i], rows, cols) / dotProductColMaj(&vecs[j], &vecs[j], rows, cols);
            for (int k = 0; k < rows; k++)
                projSum[k] += projCoeff * vecs[j + k * cols];
        }
        for (int j = 0; j < rows; j++) {
            vecs[i + j * cols] -= projSum[j];
            projSum[j] = 0;
        }
    }
}

void qrDecomposition(double* R, double* Q, double* A, int rows, int cols, double* projSum) {
    memcpy(Q, A, cols * cols * sizeof(double));
    ortVecsColMaj(Q, rows, cols, projSum);
    normalizeVecsColMaj(Q, rows, cols, 1.0);
    matrixTranspose(Q, cols, rows);
    matrixMult(R, Q, rows, cols, A, rows, cols);
    matrixTranspose(Q, cols, rows);
}

void compare(double * a, double n1){
    int n = n1;
    for (int i = 0; i < n-1; i++) {
        for (int j = 0; j < n - i - 1; j++) {
            if (a[j] > a[j + 1]) {
                double tmp = a[j];
                a[j] = a[j + 1];
                a[j + 1] = tmp;
            }
        }
    }
}

void findEigenValues(double* A, double *R, double *Q, double* projSum, int rows, int cols, int numberOfIterations, double *eigenValues) {
    for (int i = 0; i < numberOfIterations; i++) {
        qrDecomposition(R, Q, A, rows, cols, projSum);
        matrixMult(A, R, cols, cols, Q, cols, cols);
    }
    for (int i = 0; i < rows; i++) {
        eigenValues[i] = A[i * cols + i];
    }
    compare(eigenValues, rows);
}