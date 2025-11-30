// index 0  1  2  3      4      5
// Var - M1 M2 M3 gamma1 gamma2 gamma3
// index - 0 1  2  3  4  5  6 7 8
// Param - d I1 I2 I3 a1 a2 h E g0
// type - flow
#include <cmath>
#include "integrator.cpp"
// #include "linal.cpp"
#include "model.cpp"
#include <iostream>
#include <vector>
#include <fstream>


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

void copyMas(double* a, double* b){
    a[0] = b[0];
    a[1] = b[1];
    a[2] = b[2];
    a[3] = b[3];
    a[4] = b[4];
    a[5] = b[5];
}


int main(){
    void (*diffFunc)(const double*, double*, const double*) {CelticStone_6D_flow};
    double projSum[6], arg[6];
    for (int i=0;i<6;i++)
        projSum[i]=0;
    double k1[6], k2[6], k3[6], k4[6], k5[6], k6[6], k7[6], k8[6];
    

    // double M[] = {176.39671645238903,
    //             -168.85257112196294,
    //             -94.877918535603939};
    // double gamma[] = {-0.69524430378189617,
    //                 -0.40429993079272591, 
    //                 0.59428864887075838}; // Поставил минус
    // Протянул к рождению аттрактора
    double M[] = {-54.06849765, -41.5933255,   76.16298962};
    double gamma[] = {0.26638367, 0.20492122, 0.94183174}; // Тут минус уже не надо

    double tmp[6];
    
    const double step = 0.0025;
    const int time = 2000;
    const int dimension = 6;

    double params[] = {0.423, 2, 6, 7, 9, 4, 1, 744, 100};
    double delta = params[0];
    
    // Q @ M
    // tmp[0] = cos(delta) * M[0] + sin(delta) * M[1];
    // tmp[1] = -sin(delta) * M[0] + cos(delta) * M[1];
    // tmp[2] = M[2];
    
    // M[0] = tmp[0];
    // M[1] = tmp[1];
    // M[2] = tmp[2];

    // // Q @ gamma
    // tmp[3] = cos(delta) * gamma[0] + sin(delta) * gamma[1];
    // tmp[4] = -sin(delta) * gamma[0] + cos(delta) * gamma[1];
    // tmp[5] = gamma[2];
    
    // gamma[0] = tmp[3];
    // gamma[1] = tmp[4];
    // gamma[2] = tmp[5];
        
    double r[3];
    double omega[3];
    double mainTrajectory[] = {M[0], M[1], M[2], gamma[0], gamma[1], gamma[2]};
    
    // ТУТ ДЕЛАЕТСЯ ПРИВЕДЕНИЕ К УРОВНЮ ЭНЕРГИИ ЕСЛИ НАДО
    // Вектора вычисляются верно
    // calc_vector_r(r,mainTrajectory,params);
    // calc_vector_omega(r,omega,mainTrajectory,params);
    // tmp[0] = mainTrajectory[0];
    // tmp[1] = mainTrajectory[1];
    // tmp[2] = mainTrajectory[2];
    // tmp[3] = mainTrajectory[3];
    // tmp[4] = mainTrajectory[4];
    // tmp[5] = mainTrajectory[5];
    // for (int k = 0; k < 3; k++){
    //     M[k] = tmp[k] * pow(2 * (params[7] + params[8] * (r[0] * gamma[0] + r[1] * gamma[1] + r[2] * gamma[2])) 
    //     / (tmp[0] * omega[0] + tmp[1] * omega[1] + tmp[2] * omega[2]), 0.5);
    // }
    // mainTrajectory[0] = M[0];
    // mainTrajectory[1] = M[1];
    // mainTrajectory[2] = M[2];

    calc_vector_r(r,mainTrajectory,params);
    calc_vector_omega(r,omega,mainTrajectory,params);
    check_geom_integral(gamma);
    check_energy_integral(M, gamma, r, omega, params[7], params[8]);
<<<<<<< HEAD
    
    double step = 0.0025;
    int iterSkip = 50000;
    int iterNum = 5000;
    int dimension = 6;
=======

>>>>>>> 6a93a972f79559d189a0a4cfac72198c51d67a95
    double old_ps;
    double new_ps, H;
    double mainCopy[6];
    std::vector<std::vector<double>> matrix(time/step, std::vector<double>(dimension));

    for (double t=0;t<time;t+=step)
        dverkStep(mainTrajectory, dimension, diffFunc, params, step, arg, k1, k2, k3, k4, k5, k6, k7, k8);

    copyMas(mainCopy, mainTrajectory);
    std::cout << "LAST POINT" << std::endl;
    std::cout << mainCopy[0] << " " << mainCopy[1] << " " <<mainCopy[2] << " " << mainCopy[3]<< " " <<mainCopy[4] << " " << mainCopy[5] << std::endl;

    const int intervals = 500;
    const double param_step = (0.485 -0.423) / intervals;
    params[0] = params[0]+param_step;
    calc_vector_r(r,mainTrajectory,params);
    std::cout << "R" << std::endl;
    std::cout << r[0] << " " << r[1] << " " <<r[2] << std::endl;
    calc_vector_omega(r,omega,mainTrajectory,params);
    std::cout << "omega" << std::endl;
    std::cout << omega[0] << " " << omega[1] << " " << omega[2] << std::endl;

    M[0] = mainTrajectory[0];
    M[1] = mainTrajectory[1];
    M[2] = mainTrajectory[2];
    gamma[0] = mainTrajectory[3];
    gamma[1] = mainTrajectory[4];
    gamma[2] = mainTrajectory[5];
    check_energy_integral(M, gamma, r, omega, params[7], params[8]);
    check_geom_integral(gamma);

    tmp[0] = M[0];
    tmp[1] = M[1];
    tmp[2] = M[2];
    tmp[3] = gamma[0];
    tmp[4] = gamma[1];
    tmp[5] = gamma[2];
    std::cout << "NORM COEF "<< pow((2 * (params[7] + params[8] * (r[0] * gamma[0] + r[1] * gamma[1] + r[2] * gamma[2])) / (tmp[0] * omega[0] + tmp[1] * omega[1] + tmp[2] * omega[2])), 0.5) << std::endl;


    // Вывод в файл
    std::ofstream fout("matrix.txt");
    for (const auto& row : matrix) {
        for (double x : row) {
            fout << x << " ";
        }
        fout << "\n";
    }
}