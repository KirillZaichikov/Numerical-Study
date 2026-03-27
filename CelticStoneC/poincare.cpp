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
    void (*diffFuncP)(const double*, double*, const double*) {CelticStone_6D_flow_poincare};
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
    // 120.63233985265303, -58.58808404718062, -107.37176969353844, -0.7753054865740857, -0.4022390300169346, 0.4869709205288837
    double M[] = {-59.75329104, -38.86436549,  64.18350573};
    double gamma[] = {0.27811237, 0.18088813, 0.9433626}; // Тут минус уже не надо
    // double M[] = {-60.22089994, -30.39261507,  62.56014739};
    // double gamma[] = {0.30853354, 0.1557124,  0.93838196}; // Тут вбита неподвижная точка при 0.485 752
    double tmp[6];
    
    const double step = 0.0025;
    const int iterSkip = 1000;
    const int iterNum = 2000;
    const int dimension = 6;

    double params[] = {0.485, 2, 6, 7, 9, 4, 1, 752, 100};
    // double params[] = {0.485, 2, 6, 7, 9, 4, 1, 752, 100};
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
    calc_vector_r(r,mainTrajectory,params);
    calc_vector_omega(r,omega,mainTrajectory,params);
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

    double old_ps;
    double new_ps, H;
    double count =0;
    double iteration=0;
    double result_time = 0;
    double mainCopy[6];
    std::vector<std::vector<double>> matrix(iterNum+1, std::vector<double>(dimension));

    for (int i=0;i<dimension;i++)
        matrix[count][i] = mainTrajectory[i];
    while (count < iterNum){
        old_ps = mainTrajectory[1] * mainTrajectory[3] - mainTrajectory[0] * mainTrajectory[4];
        dverkStep(mainTrajectory, dimension, diffFunc, params, step, arg, k1, k2, k3, k4, k5, k6, k7, k8);
        result_time = result_time + step;

        new_ps = mainTrajectory[1] * mainTrajectory[3] - mainTrajectory[0] * mainTrajectory[4];
        if ((new_ps < 0) and (old_ps > 0)){
            copyMas(mainCopy, mainTrajectory);
            dverkStep(mainCopy, dimension, diffFuncP, params, -new_ps, arg, k1, k2, k3, k4, k5, k6, k7, k8);
            // std::cout<<mainCopy[1] * mainCopy[3] - mainCopy[0] * mainCopy[4]<<std::endl;
            if (iterSkip < iteration) {
                for (int i=0;i<dimension;i++)
                    matrix[count+1][i] = mainCopy[i];
                count++;
            }
            // std::cout << mainCopy[0] << " " << mainCopy[1] << " " <<mainCopy[2] << " " << mainCopy[3]<< " " <<mainCopy[4] << " " << mainCopy[5] << std::endl;
            iteration++;
            // std::cout << mainCopy[3] << " " << mainCopy[4] << std::endl;
        }
    }
    std::cout << "LAST POINT" << std::endl;
    std::cout << mainCopy[0] << " " << mainCopy[1] << " " <<mainCopy[2] << " " << mainCopy[3]<< " " <<mainCopy[4] << " " << mainCopy[5] << std::endl;
    std::cout << result_time << std::endl;

    calc_vector_r(r,mainTrajectory,params);
    calc_vector_omega(r,omega,mainTrajectory,params);

    M[0] = mainTrajectory[0];
    M[1] = mainTrajectory[1];
    M[2] = mainTrajectory[2];
    gamma[0] = mainTrajectory[3];
    gamma[1] = mainTrajectory[4];
    gamma[2] = mainTrajectory[5];
    check_energy_integral(M, gamma, r, omega, params[7], params[8]);
    check_geom_integral(gamma);

    // Вывод в файл
    std::ofstream fout("matrix.txt");
    for (const auto& row : matrix) {
        for (double x : row) {
            fout << x << " ";
        }
        fout << "\n";
    }
}