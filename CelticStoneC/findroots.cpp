#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "model.cpp"
#include "integrator.cpp"
// #include "linal.cpp"

#define N 6
#define EPS 1e-2
#define MAX_ITER 100
#define H 1e-5   // шаг для численного дифференцирования
#define dimension 6



// ---------- Отображение пуанкаре ---------- ВРОДЕ БЫ РАБОТАЕТ КОРРЕКТНО
int G(const double x[N],  double fx[N], void(*diffFunc)(const double*, double*, const double*), void(*diffFuncP)(const double*, double*, const double*),  
       double step, double* params,double* arg,double* k1,double* k2,double* k3,double* k4,double* k5,double* k6,double* k7,double* k8) {
    double old_ps;
    double new_ps;

    double mainCopy[6];
    fx[0] = x[0];
    fx[1] = x[1];
    fx[2] = x[2];
    fx[3] = x[3];
    fx[4] = x[4];
    fx[5] = x[5]; 
    
    std::cout << "**** start iteraions" << std::endl;
    for (int i=0;i<6;i++)
            std::cout << fx[i] << ' ';
    while (1){
        old_ps = fx[1] * fx[3] - fx[0] * fx[4];
        dverkStep(fx, dimension, diffFunc, params, step, arg, k1, k2, k3, k4, k5, k6, k7, k8);
        new_ps = fx[1] * fx[3] - fx[0] * fx[4];
        // std::cout<< "new_ps " << new_ps << "old_ps" << old_ps << std::endl;
        if ((new_ps < 0) and (old_ps > 0)){
            dverkStep(fx, dimension, diffFuncP, params, -new_ps, arg, k1, k2, k3, k4, k5, k6, k7, k8); break;}
    }
    std::cout << "poincare point" << std::endl;
    for (int i=0;i<6;i++){
        std::cout << fx[i] << ' ';
    }std::cout<<std::endl;
    fx[0] = x[0]-fx[0];
    fx[1] = x[0]-fx[1];
    fx[2] = x[0]-fx[2];
    fx[3] = x[0]-fx[3];
    fx[4] = x[0]-fx[4];
    fx[5] = x[0]-fx[5]; 
}

// ---------- Численное вычисление Якобиана ----------
void numerical_J(const double x[N], double J[N][N], void(*diffFunc)(const double*, double*, const double*), void(*diffFuncP)(const double*, double*, const double*), double step, double* params,double* arg,double* k1,double* k2,double* k3,double* k4,double* k5,double* k6,double* k7,double* k8) {
    double fx0[N], fxh[N], xh[N];
    std::cout << "start numerical J" << std::endl;
    G(x, fx0, diffFunc, diffFuncP, step, params, arg, k1, k2, k3, k4, k5, k6, k7, k8);

    std::cout << "x0-fx0" << std::endl;
    for (int i=0;i<6;i++){
        std::cout << fx0[i] << ' ';
    }std::cout<<std::endl;

    for (int j = 0; j < N; j++) {
        for (int i = 0; i < N; i++)
            xh[i] = x[i];
        xh[j] += H;
        G(xh, fxh, diffFunc, diffFuncP, step, params, arg, k1, k2, k3, k4, k5, k6, k7, k8);
        for (int i = 0; i < N; i++)
            J[i][j] = (fxh[i] - fx0[i]) / H;
    }
}


void mat_vec_mul(const double A[N][N], const double x[N], double result[N]) {
    for (int i = 0; i < N; i++) {
        result[i] = 0.0;
        for (int j = 0; j < N; j++) {
            result[i] += A[i][j] * x[j];
        }
    }
}

// ---------- Многомерный метод Ньютона ----------
int newton_numeric(double x[N], void(*diffFunc)(const double*, double*, const double*), void(*diffFuncP)(const double*, double*, const double*),
                   double step, double* params,double* arg,double* k1,double* k2,double* k3,double* k4,double* k5,double* k6,double* k7,double* k8) {
    double fx[N], Jm[N][N];
    for (int iter = 0; iter < MAX_ITER; iter++) {
        
        G(x, fx , diffFunc, diffFuncP, step, params, arg, k1, k2, k3, k4, k5, k6, k7, k8);
        
        std::cout << "x - fx point" << std::endl;
        for (int i=0;i<6;i++){
            std::cout << fx[i] << ' ';
        }std::cout<<std::endl;

        // Проверка невязки
        double norm = 0;
        for (int i = 0; i < N; i++) norm += fx[i]*fx[i];
        // std::cout<<sqrt(norm)<<std::endl;
        if (sqrt(norm) < EPS) {
            printf("get point for %d ineration.\n", iter);
            return 1;
        }

        // Численный Якобиан
        numerical_J(x, Jm, diffFunc, diffFuncP, step, params, arg, k1, k2, k3, k4, k5, k6, k7, k8);

        std::cout << "Jacobi matrix" << std::endl;
        for (int i=0;i<6;i++){
            for (int j=0;j<6;j++){
                std::cout << Jm[i][j] << ' ';
            }
            std::cout<<std::endl;
        }

        // Решаем xn+1 = xn- J-1*f(xn)
        double Jm_inv[N][N];
        // double test[N][N];
        double tmp[N];
        invert_matrix(Jm, Jm_inv);
        // multiply6x6(Jm_inv, Jm, test);
        // for (int i=0;i<6;i++){
        //     for (int j=0;j<6;j++){
        //         std::cout << test[i][j] << ' ';
        //     }
        //     std::cout<<std::endl;
        // }
        mat_vec_mul(Jm_inv, fx, tmp);
        for (int i=0;i<N;i++){
            x[i] = x[i] - tmp[i];
        }

        // Обновляем x
        double delta_norm = 0;
        for (int i = 0; i < N; i++) {
            delta_norm += x[i]*x[i];
        }

        if (sqrt(delta_norm) < EPS) {
            printf("Достигнута точность на %d итерации.\n", iter);
            return 1;
        }
    }
    printf("Method not find fixed point.\n");
    return 1; // Поставил тут 1 для дебага
}

// ---------- Пример ----------
int main() {
    void (*diffFunc)(const double*, double*, const double*) {CelticStone_6D_flow};
    void (*diffFuncP)(const double*, double*, const double*) {CelticStone_6D_flow_poincare};
    double projSum[6], arg[6];
    for (int i=0;i<6;i++)
        projSum[i]=0;
    double k1[6], k2[6], k3[6], k4[6], k5[6], k6[6], k7[6], k8[6];

// [-59.75327948 -38.86436877  64.18356172] [0.27811238 0.18088818 0.94336259]
    double params[] = {0.485, 2, 6, 7, 9, 4, 1, 752, 100};
    double M[] = {-60.22089994, -30.39261507,  62.56014739};
    double gamma[] = {0.30853354, 0.1557124,  0.93838196}; // Тут вбита неподвижная точка при 0.485 752
    double step = 0.0001;
    double initial_point[] = {M[0], M[1], M[2], gamma[0], gamma[1], gamma[2]};
    if (newton_numeric(initial_point, diffFunc, diffFuncP, step, params, arg, k1, k2, k3, k4, k5, k6, k7, k8)) {
        printf("Solution founded:\n");
        for (int i = 0; i < N; i++)
            printf("x[%d] = %.10f\n", i,  initial_point[i]);
    }
    return 0;
}
