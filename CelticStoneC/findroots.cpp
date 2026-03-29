#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "model.cpp"
#include "integrator.cpp"
// #include "linal.cpp"

#define N 6
#define EPS 1e-5
#define MAX_ITER 50
#define H 1e-5   // шаг для численного дифференцирования
#define dimension 6

void check_geom_integral(double* gamma){
    double a = pow(gamma[0], 2)+pow(gamma[1], 2)+pow(gamma[2], 2);
    if ((a > 0.999) && (a < 1.001)) {
        std::cout << "**** Geom SUCCESS " << a << std::endl;
    } else {
        std::cout << "**** Geom FAIL " << a << std::endl;
    }
}

void check_energy_integral(double* M, double* gamma, double* r, double* omega, double energy, double g0){
    double a = 0.5 * (M[0] * omega[0] + M[1] * omega[1] + M[2] * omega[2]) - g0 * (r[0] * gamma[0] + r[1] * gamma[1] + r[2] * gamma[2]);
    if ((a > energy-0.1) && (a < energy+0.1)) {
        std::cout << "**** Energy SUCCESS " << a << std::endl;
    } else {
        std::cout << "**** Energy FAIL " << a << std::endl;
    }
}

// ---------- Отображение пуанкаре ---------- ВРОДЕ БЫ РАБОТАЕТ КОРРЕКТНО
int F(const double x[N],  double fx[N], void(*diffFunc)(const double*, double*, const double*), void(*diffFuncP)(const double*, double*, const double*),  
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
    
    // std::cout << "**** start iteraions" << std::endl;
    // for (int i=0;i<6;i++)
    //         std::cout << fx[i] << ' ';
    int c = 0;
    while (1){
        old_ps = fx[1] * fx[3] - fx[0] * fx[4];
        dverkStep(fx, dimension, diffFunc, params, step, arg, k1, k2, k3, k4, k5, k6, k7, k8);
        new_ps = fx[1] * fx[3] - fx[0] * fx[4];
        c++;
        // std::cout<< "new_ps " << new_ps << "old_ps" << old_ps << std::endl;
        if ((new_ps < 0) && (old_ps > 0)&&(c>1)){
            dverkStep(fx, dimension, diffFuncP, params, -new_ps, arg, k1, k2, k3, k4, k5, k6, k7, k8); 
            // std::cout<<fx[1] * fx[3] - fx[0] * fx[4]<<std::endl;
            // c++;
            // if (c>1){
            // break;}
            break;}
            if (c>100000){
                std::cout<<"****ERROR****"<<std::endl;
            }
    }
    // std::cout << "poincare point" << std::endl;
    // for (int i=0;i<6;i++){
    //     std::cout << fx[i] << ' ';
    // }std::cout<<std::endl;
    // fx[0] = x[0]-fx[0];
    // fx[1] = x[1]-fx[1];
    // fx[2] = x[2]-fx[2];
    // fx[3] = x[3]-fx[3];
    // fx[4] = x[4]-fx[4];
    // fx[5] = x[5]-fx[5]; 
}

// ---------- Отображение пуанкаре в переменных Андуайе Депри ------------
void G(const double x[3],  double gx[3], void(*diffFunc)(const double*, double*, const double*), void(*diffFuncP)(const double*, double*, const double*),  
       double step, double* params,double* arg,double* k1,double* k2,double* k3,double* k4,double* k5,double* k6,double* k7,double* k8) {
    double initial[6], fx[6];

    double l  = x[0];
    double LG = x[1];
    double HG = x[2];
    initial[0] = sqrt(1-LG*LG) * sin(l);
    initial[1] = sqrt(1-LG*LG) * cos(l);
    initial[2] = LG;
    initial[3] = -(HG * sqrt(1-LG*LG) + LG*sqrt(1-HG*HG)) * sin(l);
    initial[4] = -(HG * sqrt(1-LG*LG) + LG*sqrt(1-HG*HG)) * cos(l);
    initial[5] = -(HG * LG - sqrt(1-LG*LG) * sqrt(1-HG*HG)) * 1; // Поставил минус для gamma

    // std::cout<<"point after varchange"<<std::endl;
    // for (int i =0;i<6;i++)
    //     std::cout<<initial[i]<<" ";
    // std::cout<<std::endl;
    // Q @ M
    double tmp[6];
    tmp[0] =  cos(params[0]) * initial[0] + sin(params[0]) * initial[1];
    tmp[1] = -sin(params[0]) * initial[0] + cos(params[0]) * initial[1];
    tmp[2] = initial[2];

    initial[0] = tmp[0];
    initial[1] = tmp[1];
    initial[2] = tmp[2];

    // Q @ gamma
    tmp[3] =  cos(params[0]) * initial[3] + sin(params[0]) * initial[4];
    tmp[4] = -sin(params[0]) * initial[3] + cos(params[0]) * initial[4];
    tmp[5] = initial[5];
    
    initial[3] = tmp[3];
    initial[4] = tmp[4];
    initial[5] = tmp[5];

    for (int i =0;i<6;i++)
        std::cout<<initial[i]<<" ";
    std::cout<<std::endl;

    double r[3];
    double omega[3];
    double tmpI[6];
    calc_vector_r(r,initial,params);
    std::cout<<"r"<<std::endl;
    for (int i =0;i<6;i++)
        std::cout<<r[i]<<" ";
    std::cout<<std::endl;
    calc_vector_omega(r,omega,initial,params);
    std::cout<<"omega"<<std::endl;
    for (int i =0;i<6;i++)
        std::cout<<omega[i]<<" ";
    std::cout<<std::endl;
    tmpI[0] = initial[0];
    tmpI[1] = initial[1];
    tmpI[2] = initial[3];
    tmpI[3] = initial[3];
    tmpI[4] = initial[3];
    tmpI[5] = initial[3];
    for (int k = 0; k < 3; k++)
        initial[k] = tmp[k] * pow((2 * (params[7] + params[8] * (r[0] * initial[3] + r[1] * initial[4] + r[2] * initial[5])) / (tmp[0] * omega[0] + tmp[1] * omega[1] + tmp[2] * omega[2])), 0.5);
    // double gamma[] = {x[3],x[4],x[5]};
    // std::cout<<"point before poincare"<<std::endl;
    // for (int i =0;i<6;i++)
    //     std::cout<<initial[i]<<" ";
    // std::cout<<std::endl;
    F(initial, fx, diffFunc, diffFuncP, step, params, arg, k1, k2, k3, k4, k5, k6, k7, k8);
    // std::cout<<"point after poincare"<<std::endl;
    // for (int i =0;i<6;i++)
    //     std::cout<<fx[i]<<" ";
    // std::cout<<std::endl;
    fx[3] = -fx[3];
    fx[4] = -fx[4];
    fx[5] = -fx[5];
    tmp[0] =  cos(-params[0]) * fx[0] + sin(-params[0]) * fx[1];
    tmp[1] = -sin(-params[0]) * fx[0] + cos(-params[0]) * fx[1];
    tmp[2] = fx[2];

    fx[0] = tmp[0];
    fx[1] = tmp[1];
    fx[2] = tmp[2];  

    // Q @ gamma
    tmp[3] =  cos(-params[0]) * fx[3] + sin(-params[0]) * fx[4];
    tmp[4] = -sin(-params[0]) * fx[3] + cos(-params[0]) * fx[4];
    tmp[5] = fx[5];
    
    fx[3] = tmp[3];
    fx[4] = tmp[4];
    fx[5] = tmp[5];

    double L = fx[2];
    double G = sqrt(fx[0]*fx[0]+fx[1]*fx[1]+fx[2]*fx[2]);
    double h = fx[0]*fx[3] + fx[1]*fx[4] + fx[2]*fx[5];
    gx[0] = atan(fx[0]/fx[1])+3.1415926;
    gx[1] = L/G;
    gx[2] = h/G;
    // std::cout<<"point after backchange"<<std::endl;
    // for (int i =0;i<3;i++)
    //     std::cout<<gx[i]<<" ";
    // std::cout<<std::endl;
    gx[0] = x[0]-gx[0];
    gx[1] = x[1]-gx[1];
    gx[2] = x[2]-gx[2];
}

// ---------- Численное вычисление Якобиана ----------
void numerical_J(const double x[3], double J[3][3], void(*diffFunc)(const double*, double*, const double*), void(*diffFuncP)(const double*, double*, const double*), double step, double* params,double* arg,double* k1,double* k2,double* k3,double* k4,double* k5,double* k6,double* k7,double* k8) {
    double gxh[3], gx0[3], xh[3];
    // std::cout << "start numerical J" << std::endl;
    G(x, gx0, diffFunc, diffFuncP, step, params, arg, k1, k2, k3, k4, k5, k6, k7, k8);
    // std::cout << "x0-fx0" << std::endl;
    // for (int i=0;i<6;i++){
    //     std::cout << fx0[i] << ' ';
    // }std::cout<<std::endl;

    for (int j = 0; j < 3; j++) {
        for (int i = 0; i <3; i++)
            xh[i] = x[i];
        xh[j] += H;
        
        G(xh, gxh, diffFunc, diffFuncP, step, params, arg, k1, k2, k3, k4, k5, k6, k7, k8);

        if (j==0){
            std::cout << "xh - fxh" << std::endl;
            for (int k = 0; k < 3; k++)
                std::cout<<gxh[k]<<' ';
            std::cout<<std::endl;
        }
        for (int i = 0; i < 3; i++){
            J[i][j] = (gxh[i] - gx0[i]) / H;
            // if (j==0){
            // std::cout << J[i][j] << ' ';}
        }std::cout<<std::endl;
    }
}


// void mat_vec_mul(const double A[N][N], const double x[N], double result[N]) {
//     for (int i = 0; i < N; i++) {
//         result[i] = 0.0;
//         for (int j = 0; j < N; j++) {
//             result[i] += A[i][j] * x[j];
//         }
//     }
// }
void mat_vec_mul(const double A[3][3], const double x[3], double result[3]) {
    for (int i = 0; i < 3; i++) {
        result[i] = 0.0;
        for (int j = 0; j < 3; j++) {
            result[i] += A[i][j] * x[j];
        }
    }
}

// ---------- Многомерный метод Ньютона ----------
int newton_numeric(double x[3],void(*diffFunc)(const double*, double*, const double*), void(*diffFuncP)(const double*, double*, const double*),
                   double step, double* params,double* arg,double* k1,double* k2,double* k3,double* k4,double* k5,double* k6,double* k7,double* k8) {
    double gx[3], Jm[3][3];
    for (int iter = 0; iter < MAX_ITER; iter++) {
        // std::cout << "start x" << std::endl;
        // for (int i=0;i<3;i++){
        //     std::cout << x[i] << ' ';
        // }std::cout<<std::endl;

        //G(x)
        G(x, gx, diffFunc, diffFuncP, step, params, arg, k1, k2, k3, k4, k5, k6, k7, k8);
        std::cout << "x - fx point" << std::endl;
        for (int i=0;i<3;i++){
            std::cout << gx[i] << ' ';
        }std::cout<<std::endl;

        // Проверка невязки
        double norm = 0;
        for (int i = 0; i < 3; i++) norm += gx[i]*gx[i];
        std::cout<<sqrt(norm)<<std::endl;
        if ((sqrt(norm) < EPS) && (iter != 0)){
            printf("get point for %d ineration.\n", iter);
            return 1;
        }

        // Численный Якобиан
        // std::cout<<gx[0]<<" "<<gx[1]<<' '<<gx[2]<<std::endl;
        numerical_J(x, Jm, diffFunc, diffFuncP, step, params, arg, k1, k2, k3, k4, k5, k6, k7, k8);
        // std::cout << "Jacobi matrix" << std::endl;
        // for (int i=0;i<3;i++){
        //     for (int j=0;j<3;j++){
        //         std::cout << Jm[i][j] << ' ';
        //     }
        //     std::cout<<std::endl;
        // }

        // Решаем xn+1 = xn- J^{-1}*g(xn)
        double Jm_inv[3][3];
        double test[3][3];
        double tmp[3];
        inverse3x3(Jm, Jm_inv);
        multiply3x3(Jm_inv, Jm, test);
        // for (int i=0;i<3;i++){
        //     for (int j=0;j<3;j++){
        //         std::cout << test[i][j] << ' ';
        //     }
        //     std::cout<<std::endl;
        // }
        // std::cout<<"invert matrix"<<std::endl;
        // for (int i=0;i<3;i++){
        //     for (int j=0;j<3;j++){
        //         std::cout << Jm_inv[i][j] << ' ';
        //     }
        //     std::cout<<std::endl;
        // }
        // std::cout<<gx[0]<<" "<<gx[1]<<' '<<gx[2]<<std::endl;
        mat_vec_mul(Jm_inv, gx, tmp);
        // std::cout<<tmp[0]<<" "<<tmp[1]<<' '<<tmp[2]<<std::endl;

        //Думаю что перед тем как вычесть точки нужно обе привести к одному уровню энергии
        // Нормировка интегралов
        // double gamma_norm = sqrt(x[3]*x[3]+x[4]*x[4]+x[5]*x[5]);
        // for (int i=3;i<6;i++)
        //     x[i] = x[i] / gamma_norm;
        // double r[3];
        // double omega[3];
        // double tmpI[6];
        // calc_vector_r(r,x,params);
        // calc_vector_omega(r,omega,x,params);
        // tmpI[0] = x[0];
        // tmpI[1] = x[1];
        // tmpI[2] = x[2];
        // tmpI[3] = x[3];
        // tmpI[4] = x[4];
        // tmpI[5] = x[5];
        // for (int k = 0; k < 3; k++)
        //     x[k] = tmp[k] * pow((2 * (params[7] + params[8] * (r[0] * x[3] + r[1] * x[4] + r[2] * x[5])) / (tmp[0] * omega[0] + tmp[1] * omega[1] + tmp[2] * omega[2])), 0.5);
        // double gamma[] = {x[3],x[4],x[5]};
        // double M[] = {x[0],x[1],x[2]};
        // check_geom_integral(gamma);
        // check_energy_integral(M, gamma, r, omega, params[7], params[8]);

        for (int i=0;i<3;i++){
            x[i] = x[i] - tmp[i];
        }

        std::cout << "result x" << std::endl;
        for (int i=0;i<3;i++){
            std::cout << x[i] << ' ';
        }std::cout<<std::endl;

        // Проверим норму
        double delta_norm = 0;
        for (int i = 0; i < 3; i++) {
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
    // double M[] = {-59.75329018, -38.86436555,  64.18351007};
    // double gamma[] = {0.27811238, 0.18088813, 0.9433626}; // Тут вбита неподвижная точка при 0.485 752 достигается за 0 итераций
    double M[] = {-60.63246519, -33.62842703,  62.43066679};
    double gamma[] = {0.29012771, 0.16091278, 0.9433626}; // Тут вбита близкая к неподвижной
    double step = 0.001;
    double initial_point[] = {M[0], M[1], M[2], gamma[0], gamma[1], gamma[2]};
    double initial_pointAD[] = {3.721,0.6691441646,-0.3847013110};
    if (newton_numeric(initial_pointAD, diffFunc, diffFuncP, step, params, arg, k1, k2, k3, k4, k5, k6, k7, k8)) {
        printf("Solution founded:\n");
        for (int i = 0; i < 3; i++)
            printf("x[%d] = %.10f\n", i,  initial_pointAD[i]);
    }
    return 0;
}
