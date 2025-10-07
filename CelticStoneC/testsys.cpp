// index 0  1  2  3      4      5
// Var - M1 M2 M3 gamma1 gamma2 gamma3
// index - 0 1  2  3  4  5  6 7 8
// Param - d I1 I2 I3 a1 a2 h E g0
// type - flow
#include <cmath>
#include "integrator.cpp"
// #include "linal.cpp"
// #include "model.cpp"
#include <iostream>
#include <vector>
#include <fstream>




int inverse3x3(double A[3][3], double inv[3][3]) {
    double det =
        A[0][0] * (A[1][1] * A[2][2] - A[1][2] * A[2][1])
        - A[0][1] * (A[1][0] * A[2][2] - A[1][2] * A[2][0])
        + A[0][2] * (A[1][0] * A[2][1] - A[1][1] * A[2][0]);

    if (det == 0) {
        return 0; // матрица необратима
    }

    inv[0][0] = (A[1][1] * A[2][2] - A[1][2] * A[2][1]) / det;
    inv[0][1] = -(A[0][1] * A[2][2] - A[0][2] * A[2][1]) / det;
    inv[0][2] = (A[0][1] * A[1][2] - A[0][2] * A[1][1]) / det;

    inv[1][0] = -(A[1][0] * A[2][2] - A[1][2] * A[2][0]) / det;
    inv[1][1] = (A[0][0] * A[2][2] - A[0][2] * A[2][0]) / det;
    inv[1][2] = -(A[0][0] * A[1][2] - A[0][2] * A[1][0]) / det;

    inv[2][0] = (A[1][0] * A[2][1] - A[1][1] * A[2][0]) / det;
    inv[2][1] = -(A[0][0] * A[2][1] - A[0][1] * A[2][0]) / det;
    inv[2][2] = (A[0][0] * A[1][1] - A[0][1] * A[1][0]) / det;

    return 1;
}

void calc_vector_r(double* r, const double* state, const double* params) {
    // Вычисление вектора r
    double gamma[3];
    gamma[0] = state[3];
    gamma[1] = state[4];
    gamma[2] = state[5];
    r[0] = -params[4] * gamma[0] / gamma[2];
    r[1] = -params[5] * gamma[1] / gamma[2];
    r[2] = -params[6] + (params[4] * pow(gamma[0], 0.2e1) + params[5] * pow(gamma[1], 0.2e1)) * pow(gamma[2], -0.2e1) / 0.2e1;
    // std::cout << "vector r " << r[0] << " " << r[1] << " " << r[2] << std::endl;
}

void calc_vector_omega(const double* r, double* omega, const double* state, const double* params) {
        // Вычисление вектора r
    double J[3];
    J[0] = params[1];
    J[1] = params[2];
    J[2] = params[3];
    double a1 = params[4];
    double a2 = params[5];
    double delta = params[0];
    double g0 = params[8];

    double M[3];
    double JQ_kuz[3][3];
    double JQ_rev[3][3];

    M[0] = state[0];
    M[1] = state[1];
    M[2] = state[2];

    JQ_kuz[0][0] = pow(cos(delta), 0.2e1) * J[0] + pow(sin(delta), 0.2e1) * J[1] + pow(r[1], 0.2e1) + pow(r[2], 0.2e1);
    JQ_kuz[0][1] = -cos(delta) * J[0] * sin(delta) + sin(delta) * J[1] * cos(delta) - r[0] * r[1];
    JQ_kuz[0][2] = -r[0] * r[2];
    JQ_kuz[1][0] = -cos(delta) * J[0] * sin(delta) + sin(delta) * J[1] * cos(delta) - r[0] * r[1];
    JQ_kuz[1][1] = pow(sin(delta), 0.2e1) * J[0] + pow(cos(delta), 0.2e1) * J[1] + pow(r[0], 0.2e1) + pow(r[2], 0.2e1);
    JQ_kuz[1][2] = -r[1] * r[2];
    JQ_kuz[2][0] = -r[0] * r[2];
    JQ_kuz[2][1] = -r[1] * r[2];
    JQ_kuz[2][2] = pow(r[0], 0.2e1) + pow(r[1], 0.2e1) + J[2];
    // std::cout << "vector r " << r[0] << " " << r[1] << " " << r[2] << std::endl;

    inverse3x3(JQ_kuz, JQ_rev);

    omega[0] = JQ_rev[0][0] * M[0] + JQ_rev[0][1] * M[1] + JQ_rev[0][2] * M[2];
    omega[1] = JQ_rev[1][0] * M[0] + JQ_rev[1][1] * M[1] + JQ_rev[1][2] * M[2];
    omega[2] = JQ_rev[2][0] * M[0] + JQ_rev[2][1] * M[1] + JQ_rev[2][2] * M[2];
}

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

void CelticStone_6D_flow(const double* state, double* res, const double* params) {
    double r[3];
    double r1, r2, r3;
    double omega[3];

    double JQ_kuz[3][3];
    double JQ_rev[3][3];
    double Jr[3][3];

    double J[3];
    J[0] = params[1];
    J[1] = params[2];
    J[2] = params[3];
    double a1 = params[4];
    double a2 = params[5];
    double delta = params[0];
    double g0 = params[8];

    // Начальные условия и фазовые
    double M[3], gamma[3];
    double tmp[6];
    M[0] = state[0];
    M[1] = state[1];
    M[2] = state[2];
    gamma[0] = state[3];
    gamma[1] = state[4];
    gamma[2] = state[5];
    double dr[3], dM[3], dgamma[3];

    // Вычисление вектора r
    r[0] = -params[4] * gamma[0] / gamma[2];
    r[1] = -params[5] * gamma[1] / gamma[2];
    r[2] = -params[6] + (params[4] * pow(gamma[0], 0.2e1) + params[5] * pow(gamma[1], 0.2e1)) * pow(gamma[2], -0.2e1) / 0.2e1;

    // матрица М от омега
    JQ_kuz[0][0] = pow(cos(delta), 0.2e1) * J[0] + pow(sin(delta), 0.2e1) * J[1] + pow(r[1], 0.2e1) + pow(r[2], 0.2e1);
    JQ_kuz[0][1] = -cos(delta) * J[0] * sin(delta) + sin(delta) * J[1] * cos(delta) - r[0] * r[1];
    JQ_kuz[0][2] = -r[0] * r[2];
    JQ_kuz[1][0] = -cos(delta) * J[0] * sin(delta) + sin(delta) * J[1] * cos(delta) - r[0] * r[1];
    JQ_kuz[1][1] = pow(sin(delta), 0.2e1) * J[0] + pow(cos(delta), 0.2e1) * J[1] + pow(r[0], 0.2e1) + pow(r[2], 0.2e1);
    JQ_kuz[1][2] = -r[1] * r[2];
    JQ_kuz[2][0] = -r[0] * r[2];
    JQ_kuz[2][1] = -r[1] * r[2];
    JQ_kuz[2][2] = pow(r[0], 0.2e1) + pow(r[1], 0.2e1) + J[2];

    // Обратная матрица М от омега
    inverse3x3(JQ_kuz, JQ_rev);

    // Посчитали вектор omega
    omega[0] = JQ_rev[0][0] * M[0] + JQ_rev[0][1] * M[1] + JQ_rev[0][2] * M[2];
    omega[1] = JQ_rev[1][0] * M[0] + JQ_rev[1][1] * M[1] + JQ_rev[1][2] * M[2];
    omega[2] = JQ_rev[2][0] * M[0] + JQ_rev[2][1] * M[1] + JQ_rev[2][2] * M[2];

    // Отнормируем к необходимому уровню энергии
    tmp[0] = M[0];
    tmp[1] = M[1];
    tmp[2] = M[2];
    tmp[3] = gamma[0];
    tmp[4] = gamma[1];
    tmp[5] = gamma[2];
    for (int k = 0; k < 3; k++)
        M[k] = tmp[k] * pow((2 * (params[7] + params[8] * (r[0] * gamma[0] + r[1] * gamma[1] + r[2] * gamma[2])) / (tmp[0] * omega[0] + tmp[1] * omega[1] + tmp[2] * omega[2])), 0.5);

    // Посчитали вектор dgamma
    dgamma[0] = gamma[1] * omega[2] - gamma[2] * omega[1];
    dgamma[1] = -gamma[0] * omega[2] + gamma[2] * omega[0];
    dgamma[2] = gamma[0] * omega[1] - gamma[1] * omega[0];

    // Посчитали якобиан r
    Jr[0][0] = -a1 / gamma[2];
    Jr[0][1] = 0;
    Jr[0][2] = a1 * gamma[0] * pow(gamma[2], -0.2e1);
    Jr[1][0] = 0;
    Jr[1][1] = -a2 / gamma[2];
    Jr[1][2] = a2 * gamma[1] * pow(gamma[2], -0.2e1);
    Jr[2][0] = a1 * gamma[0] * pow(gamma[2], -0.2e1);
    Jr[2][1] = a2 * gamma[1] * pow(gamma[2], -0.2e1);
    Jr[2][2] = -(a1 * pow(gamma[0], 0.2e1) + a2 * pow(gamma[1], 0.2e1)) * pow(gamma[2], -0.3e1);

    // Посчитали dr
    dr[0] = Jr[0][0] * dgamma[0] + Jr[0][1] * dgamma[1] + Jr[0][2] * dgamma[2];
    dr[1] = Jr[1][0] * dgamma[0] + Jr[1][1] * dgamma[1] + Jr[1][2] * dgamma[2];
    dr[2] = Jr[2][0] * dgamma[0] + Jr[2][1] * dgamma[1] + Jr[2][2] * dgamma[2];

    // Посчитали dM
    dM[0] = M[1] * omega[2] - M[2] * omega[1] + dr[1] * (omega[0] * r[1] - omega[1] * r[0]) - dr[2] * (-omega[0] * r[2] + omega[2] * r[0]) + g0 * (-r[2] * gamma[1] + r[1] * gamma[2]);
    dM[1] = -M[0] * omega[2] + M[2] * omega[0] - dr[0] * (omega[0] * r[1] - omega[1] * r[0]) + dr[2] * (omega[1] * r[2] - omega[2] * r[1]) + g0 * (r[2] * gamma[0] - r[0] * gamma[2]);
    dM[2] = M[0] * omega[1] - M[1] * omega[0] + dr[0] * (-omega[0] * r[2] + omega[2] * r[0]) - dr[1] * (omega[1] * r[2] - omega[2] * r[1]) + g0 * (-r[1] * gamma[0] + r[0] * gamma[1]);

    res[0] = dM[0];
    res[1] = dM[1];
    res[2] = dM[2];
    res[3] = dgamma[0];
    res[4] = dgamma[1];
    res[5] = dgamma[2];
}

void CelticStone_6D_flow_poincare(const double* state, double* res, const double* params) {
    double r[3];
    double r1, r2, r3;
    double omega[3];

    double JQ_kuz[3][3];
    double JQ_rev[3][3];
    double Jr[3][3];

    double J[3];
    J[0] = params[1];
    J[1] = params[2];
    J[2] = params[3];
    double a1 = params[4];
    double a2 = params[5];
    double delta = params[0];
    double g0 = params[8];

    // Начальные условия и фазовые
    double M[3], gamma[3];
    double tmp[6];
    M[0] = state[0];
    M[1] = state[1];
    M[2] = state[2];
    gamma[0] = state[3];
    gamma[1] = state[4];
    gamma[2] = state[5];
    double dr[3], dM[3], dgamma[3];

    // Вычисление вектора r
    r[0] = -params[4] * gamma[0] / gamma[2];
    r[1] = -params[5] * gamma[1] / gamma[2];
    r[2] = -params[6] + (params[4] * pow(gamma[0], 0.2e1) + params[5] * pow(gamma[1], 0.2e1)) * pow(gamma[2], -0.2e1) / 0.2e1;

    // матрица М от омега
    JQ_kuz[0][0] = pow(cos(delta), 0.2e1) * J[0] + pow(sin(delta), 0.2e1) * J[1] + pow(r[1], 0.2e1) + pow(r[2], 0.2e1);
    JQ_kuz[0][1] = -cos(delta) * J[0] * sin(delta) + sin(delta) * J[1] * cos(delta) - r[0] * r[1];
    JQ_kuz[0][2] = -r[0] * r[2];
    JQ_kuz[1][0] = -cos(delta) * J[0] * sin(delta) + sin(delta) * J[1] * cos(delta) - r[0] * r[1];
    JQ_kuz[1][1] = pow(sin(delta), 0.2e1) * J[0] + pow(cos(delta), 0.2e1) * J[1] + pow(r[0], 0.2e1) + pow(r[2], 0.2e1);
    JQ_kuz[1][2] = -r[1] * r[2];
    JQ_kuz[2][0] = -r[0] * r[2];
    JQ_kuz[2][1] = -r[1] * r[2];
    JQ_kuz[2][2] = pow(r[0], 0.2e1) + pow(r[1], 0.2e1) + J[2];

    // Обратная матрица М от омега
    inverse3x3(JQ_kuz, JQ_rev);

    // Посчитали вектор omega
    omega[0] = JQ_rev[0][0] * M[0] + JQ_rev[0][1] * M[1] + JQ_rev[0][2] * M[2];
    omega[1] = JQ_rev[1][0] * M[0] + JQ_rev[1][1] * M[1] + JQ_rev[1][2] * M[2];
    omega[2] = JQ_rev[2][0] * M[0] + JQ_rev[2][1] * M[1] + JQ_rev[2][2] * M[2];

    // Отнормируем к необходимому уровню энергии
    tmp[0] = M[0];
    tmp[1] = M[1];
    tmp[2] = M[2];
    tmp[3] = gamma[0];
    tmp[4] = gamma[1];
    tmp[5] = gamma[2];
    for (int k = 0; k < 3; k++)
        M[k] = tmp[k] * pow((2 * (params[7] + params[8] * (r[0] * gamma[0] + r[1] * gamma[1] + r[2] * gamma[2])) / (tmp[0] * omega[0] + tmp[1] * omega[1] + tmp[2] * omega[2])), 0.5);

    // Посчитали вектор dgamma
    dgamma[0] = gamma[1] * omega[2] - gamma[2] * omega[1];
    dgamma[1] = -gamma[0] * omega[2] + gamma[2] * omega[0];
    dgamma[2] = gamma[0] * omega[1] - gamma[1] * omega[0];

    // Посчитали якобиан r
    Jr[0][0] = -a1 / gamma[2];
    Jr[0][1] = 0;
    Jr[0][2] = a1 * gamma[0] * pow(gamma[2], -0.2e1);
    Jr[1][0] = 0;
    Jr[1][1] = -a2 / gamma[2];
    Jr[1][2] = a2 * gamma[1] * pow(gamma[2], -0.2e1);
    Jr[2][0] = a1 * gamma[0] * pow(gamma[2], -0.2e1);
    Jr[2][1] = a2 * gamma[1] * pow(gamma[2], -0.2e1);
    Jr[2][2] = -(a1 * pow(gamma[0], 0.2e1) + a2 * pow(gamma[1], 0.2e1)) * pow(gamma[2], -0.3e1);

    // Посчитали dr
    dr[0] = Jr[0][0] * dgamma[0] + Jr[0][1] * dgamma[1] + Jr[0][2] * dgamma[2];
    dr[1] = Jr[1][0] * dgamma[0] + Jr[1][1] * dgamma[1] + Jr[1][2] * dgamma[2];
    dr[2] = Jr[2][0] * dgamma[0] + Jr[2][1] * dgamma[1] + Jr[2][2] * dgamma[2];

    // Посчитали dM
    dM[0] = M[1] * omega[2] - M[2] * omega[1] + dr[1] * (omega[0] * r[1] - omega[1] * r[0]) - dr[2] * (-omega[0] * r[2] + omega[2] * r[0]) + g0 * (-r[2] * gamma[1] + r[1] * gamma[2]);
    dM[1] = -M[0] * omega[2] + M[2] * omega[0] - dr[0] * (omega[0] * r[1] - omega[1] * r[0]) + dr[2] * (omega[1] * r[2] - omega[2] * r[1]) + g0 * (r[2] * gamma[0] - r[0] * gamma[2]);
    dM[2] = M[0] * omega[1] - M[1] * omega[0] + dr[0] * (-omega[0] * r[2] + omega[2] * r[0]) - dr[1] * (omega[1] * r[2] - omega[2] * r[1]) + g0 * (-r[1] * gamma[0] + r[0] * gamma[1]);

    double H;
    H = dM[1] * gamma[0] + M[1] * dgamma[0] - dM[0] * gamma[1] - M[0] * dgamma[1];//M2' * gamma1 + M2 * gamma1' - M1' * gamma2 - M1 * gamma2'

    res[0] = dM[0] / H;
    res[1] = dM[1] / H;
    res[2] = dM[2] / H;
    res[3] = dgamma[0] / H;
    res[4] = dgamma[1] / H;
    res[5] = dgamma[2] / H;
}


int main(){
    void (*diffFunc)(const double*, double*, const double*) {CelticStone_6D_flow};
    void (*diffFuncP)(const double*, double*, const double*) {CelticStone_6D_flow_poincare};
    double projSum[6], arg[6];
    for (int i=0;i<6;i++)
    projSum[i]=0;
    double k1[6], k2[6], k3[6], k4[6], k5[6], k6[6], k7[6], k8[6];
    
    double params[] = {0.485, 2, 6, 7, 9, 4, 1, 752, 100};
    double delta = params[0];
    double M[] = {176.39671645238903,
                -168.85257112196294,
                -94.877918535603939};
    double gamma[] = {-0.69524430378189617,
                    -0.40429993079272591, 
                    0.59428864887075838}; // Поставил минус
    double tmp[6];
    
    
    // Q @ M
    tmp[0] = cos(delta) * M[0] + sin(delta) * M[1];
    tmp[1] = -sin(delta) * M[0] + cos(delta) * M[1];
    tmp[2] = M[2];
    
    M[0] = tmp[0];
    M[1] = tmp[1];
    M[2] = tmp[2];

    // Q @ gamma
    tmp[3] = cos(delta) * gamma[0] + sin(delta) * gamma[1];
    tmp[4] = -sin(delta) * gamma[0] + cos(delta) * gamma[1];
    tmp[5] = gamma[2];
    
    gamma[0] = tmp[3];
    gamma[1] = tmp[4];
    gamma[2] = tmp[5];
        
    double r[3];
    double omega[3];
    for (int i =0;i<3;i++)
        std::cout << M[i] << std::endl;
    for (int i =0;i<3;i++)
        std::cout << gamma[i] << std::endl;
    double mainTrajectory[] = {M[0], M[1], M[2], gamma[0], gamma[1], gamma[2]};

    calc_vector_r(r,mainTrajectory,params);
    calc_vector_omega(r,omega,mainTrajectory,params);
    check_geom_integral(gamma);
    check_energy_integral(M, gamma, r, omega, params[7], params[8]);
    
    double step = 0.0025;
    int iterSkip = 0;
    int iterNum = 200000;
    int dimension = 6;
    double old_ps;
    double new_ps, H;
    int count =0;
    double iteration=0;
    
    double result_time = 0;
    double mainCopy[6];
    std::vector<std::vector<double>> matrix(iterNum, std::vector<double>(dimension));

    while (count < iterNum){
        old_ps = mainTrajectory[1] * mainTrajectory[3] - mainTrajectory[0] * mainTrajectory[4];
        dverkStep(mainTrajectory, dimension, diffFunc, params, step, arg, k1, k2, k3, k4, k5, k6, k7, k8);
        result_time = result_time + step;

        new_ps = mainTrajectory[1] * mainTrajectory[3] - mainTrajectory[0] * mainTrajectory[4];
        if ((new_ps > 0) and (old_ps < 0)){
            copyMas(mainCopy, mainTrajectory);
            dverkStep(mainCopy, dimension, diffFuncP, params, -new_ps, arg, k1, k2, k3, k4, k5, k6, k7, k8);
            // std::cout<<mainCopy[1] * mainCopy[3] - mainCopy[0] * mainCopy[4]<<std::endl;
            if (iterSkip < iteration) {
                for (int i=0;i<dimension;i++)
                    matrix[count][i] = mainCopy[i];
                count++;
            }
            if (count % 2000 == 0)
                std::cout<<count<<std::endl;
            iteration++;
        }
    }
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