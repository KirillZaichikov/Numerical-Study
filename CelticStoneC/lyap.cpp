// index 0  1  2  3      4      5
// Var - M1 M2 M3 gamma1 gamma2 gamma3
// index - 0 1  2  3  4  5  6 7 8
// Param - d I1 I2 I3 a1 a2 h E g0
// type - flow
#include <cmath>
// #include "integrator.cpp"
#include "linal.cpp"
#include "model.cpp"
#include <iostream>


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
    if ((a > energy-0.01) && (a < energy+0.01)) {
        std::cout << "**** Energy SUCCESS " << a << std::endl;
    } else {
        std::cout << "**** Energy FAIL " << a << std::endl;
    }
}

void dverkStep(double* val, const int dimension, void(*diffFunc)(const double*, double*, const double*), const double* params, double step,
    double* arg, double* k1, double* k2, double* k3, double* k4, double* k5, double* k6, double* k7, double* k8) {

    diffFunc(val, k1, params);
    for (int j = 0; j < dimension; j++) {
        arg[j] = val[j] + step / 6. * k1[j];
    }
    diffFunc(arg, k2, params);
    for (int j = 0; j < dimension; j++) {
        arg[j] = val[j] + step * (4. / 75. * k1[j] + 16. / 75. * k2[j]);
    }
    diffFunc(arg, k3, params);
    for (int j = 0; j < dimension; j++) {
        arg[j] = val[j] + step * (5. / 6. * k1[j] - 8. / 3. * k2[j] + 5. / 2. * k3[j]);
    }
    diffFunc(arg, k4, params);
    for (int j = 0; j < dimension; j++) {
        arg[j] = val[j] + step * (-165. / 64. * k1[j] + 55. / 6. * k2[j] - 425. / 64. * k3[j] + 85. / 96. * k4[j]);
    }
    diffFunc(arg, k5, params);
    for (int j = 0; j < dimension; j++) {
        arg[j] = val[j] + step * (12. / 5. * k1[j] - 8. * k2[j] + 4015. / 612. * k3[j] - 11. / 36. * k4[j] + 88. / 255. * k5[j]);
    }
    diffFunc(arg, k6, params);
    for (int j = 0; j < dimension; j++) {
        arg[j] = val[j] + step * (-8263. / 15000. * k1[j] + 124. / 75. * k2[j] - 643. / 680. * k3[j] - 81. / 250. * k4[j] + 2484. / 10625. * k5[j]);
    }
    diffFunc(arg, k7, params);
    for (int j = 0; j < dimension; j++) {
        arg[j] = val[j] + step * (3501. / 1720. * k1[j] - 300. / 43. * k2[j] + 297275. / 52632. * k3[j] - 319. / 2322. * k4[j] + 24068. / 84065. * k5[j] + 3850. / 26703. * k7[j]);
    }
    diffFunc(arg, k8, params);
    for (int j = 0; j < dimension; j++) {
        val[j] = val[j] + step * (3. / 40. * k1[j] + 875. / 2244. * k3[j] + 23. / 72. * k4[j] + 264. / 1955. * k5[j] + 125. / 11592. * k7[j] + 43. / 616. * k8[j]);
    }
}

int main(){
    std::cout<<"****LYAP CALC"<<std::endl;
    // Технические переменные
    void (*diffFunc)(const double*, double*, const double*) {CelticStone_6D_flow};
    double projSum[6], arg[6];
    for (int i=0;i<6;i++)
        projSum[i]=0;
    double k1[6], k2[6], k3[6], k4[6], k5[6], k6[6], k7[6], k8[6];
    double mainTrajectory[6];
    double slaveTrajectories[6*6];
    double lyapExp[6];
    double tmp[6];

    // Начальная точка в хаосовских координатах
    double M[] = {-59.70209864, -40.90167698,  62.21294067};
    double gamma[] = {0.27785914, 0.19036022, 0.94157171}; // Поставил мину
    

    double step = 0.0025;
    double timeSkip = 0;
    double timeSkipClP = 500;
    double calcTime = 200;
    int dimension = 6;
    int expsNum = 6;
    double eps = 1e-7;
    
    double params[] = {0.4892, 2, 6, 7, 9, 4, 1, 752, 100};
    double delta = params[0];
    // 77.3333, -231.617, -94.8779, -0.803553, -0.0335453, 0.594289
    // // Q @ M
    // mainTrajectory[0] = cos(delta) * M[0] + sin(delta) * M[1];
    // mainTrajectory[1] = -sin(delta) * M[0] + cos(delta) * M[1];
    // mainTrajectory[2] = M[2];
    // // Q @ gamma
    // mainTrajectory[3] = cos(delta) * gamma[0] + sin(delta) * gamma[1];
    // mainTrajectory[4] = -sin(delta) * gamma[0] + cos(delta) * gamma[1];
    // mainTrajectory[5] = gamma[2];
    // gamma[0] = mainTrajectory[3];
    // gamma[1] = mainTrajectory[4];
    // gamma[2] = mainTrajectory[5];

    mainTrajectory[0] = M[0];
    mainTrajectory[1] = M[1];
    mainTrajectory[2] = M[2];
    mainTrajectory[3] = gamma[0];
    mainTrajectory[4] = gamma[1];
    mainTrajectory[5] = gamma[2];

    double r[3];
    double omega[3];
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
    

    calc_vector_r(r,mainTrajectory,params);
    calc_vector_omega(r,omega,mainTrajectory,params);
    check_geom_integral(gamma);
    check_energy_integral(M, gamma, r, omega, params[7], params[8]);

    for (int i=0; i<6;i++)
        std::cout << mainTrajectory[i] << std::endl;

    for (double t = 0; t < timeSkip-step; t += step) // Skip points to reach the attractor
        dverkStep(mainTrajectory, dimension, diffFunc, params, step, arg, k1, k2, k3, k4, k5, k6, k7, k8);
    std::cout<<"****skip end" << std::endl;

    for (int i = 0; i < expsNum; i++){
        for (int j = 0; j < dimension; j++)
            slaveTrajectories[i * dimension + j] = mainTrajectory[j];
        slaveTrajectories[i * dimension + i] += eps;}

    for (double t = 0; t < timeSkipClP-step; t += step) {
        dverkStep(mainTrajectory, dimension, diffFunc, params, step, arg, k1, k2, k3, k4, k5, k6, k7, k8);

        for (int i = 0; i < expsNum; i++) {
            dverkStep(&slaveTrajectories[i * dimension], dimension, diffFunc, params, step, arg, k1, k2, k3, k4, k5, k6, k7, k8);
            for (int j = 0; j < dimension; j++)
                slaveTrajectories[i * dimension + j] -= mainTrajectory[j];
        }

        ortVecs(slaveTrajectories, dimension, expsNum, projSum);
        normalizeVecs(slaveTrajectories, dimension, expsNum, eps);

        for (int32_t i = 0; i < expsNum; i++) {
            for (int32_t j = 0; j < dimension; j++)
                slaveTrajectories[i * dimension + j] += mainTrajectory[j];
        }
    } std::cout<<"****warm end" << std::endl;
    
    //////////////////////////////////////////////////////////
    // std::cout << "main" << std::endl;
    // for (int32_t j = 0; j < dimension; j++)
    //     std::cout << mainTrajectory[j] << " "<< std::endl;
    // std::cout << std::endl;

    for (int i = 0; i < expsNum; i++)
        for (int j = 0; j < dimension; j++)
            slaveTrajectories[i * dimension + j] -= mainTrajectory[j];

    // // std::cout << "vectors" << std::endl;
    // // for (int32_t i = 0; i < expsNum; i++){
    // //     for (int32_t j = 0; j < dimension; j++)
    // //         std::cout << slaveTrajectories[i*dimension+j] << " ";
    // //     std::cout << std::endl;}

    for (int i = 0; i < expsNum; i++)
        for (int j = 0; j < dimension; j++)
            slaveTrajectories[i * dimension + j] += mainTrajectory[j];
    ///////////////////////////////////////////////////////////////////
    for (double t = 0; t < calcTime; t += step) {
        dverkStep(mainTrajectory, dimension, diffFunc, params, step, arg, k1, k2, k3, k4, k5, k6, k7, k8);
        for (int i = 0; i < expsNum; i++){
            dverkStep(&slaveTrajectories[i * dimension], dimension, diffFunc, params, step, arg, k1, k2, k3, k4, k5, k6, k7, k8);
            for (int j = 0; j < dimension; j++)
                slaveTrajectories[i * dimension + j] -= mainTrajectory[j];
            }
            
        ortVecs(slaveTrajectories, dimension, expsNum, projSum);
        for (int i = 0; i < expsNum; i++)
            lyapExp[i] += log(vecNorm(&slaveTrajectories[i * dimension], dimension) / eps);
        normalizeVecs(slaveTrajectories, dimension, expsNum, eps);
        // for (int32_t i = 0; i < expsNum; i++)
        //     std::cout << i << " " << vecNorm(&slaveTrajectories[i * dimension], dimension) << std::endl;
            
        for (int32_t i = 0; i < expsNum; i++) {
            for (int32_t j = 0; j < dimension; j++)
                slaveTrajectories[i * dimension + j] += mainTrajectory[j];
        }
    }std::cout<<"****calc end" << std::endl;
    ///////////////////////////////////////////////////////
    // std::cout << "main" << std::endl;
    // for (int32_t j = 0; j < dimension; j++)
    //     std::cout << mainTrajectory[j] << " "<< std::endl;
    // std::cout << std::endl;
    ////////////////////////////////////////////////////////

    for (int i=0;i<expsNum;i++)
        std::cout << "L" << i << ": " << lyapExp[i] / calcTime << std::endl;
}