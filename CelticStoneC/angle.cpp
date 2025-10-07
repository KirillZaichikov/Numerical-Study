// index 0  1  2  3      4      5
// Var - M1 M2 M3 gamma1 gamma2 gamma3
// index - 0 1  2  3  4  5  6 7 8
// Param - d I1 I2 I3 a1 a2 h E g0
// type - flow
#include <cmath>
#include "integrator.cpp"
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
    if ((a > energy-0.1) && (a < energy+0.1)) {
        std::cout << "**** Energy SUCCESS " << a << std::endl;
    } else {
        std::cout << "**** Energy FAIL " << a << std::endl;
    }
}

int main(){
    std::cout<<"****BACKWARD LYAP CALC"<<std::endl;
    // Технические переменные
    void (*diffFunc)(const double*, double*, const double*) {CelticStone_6D_flow};
    double projSum[6], arg[6];
    for (int i=0;i<6;i++)
        projSum[i]=0;
    double k1[6], k2[6], k3[6], k4[6], k5[6], k6[6], k7[6], k8[6];
    double mainTrajectory[6];
    double slaveTrajectories[6*6];
    double lyapExp[6];
    double lyapExpBw[6];
    double tmp[6];

    // Начальная точка в хаосовских координатах
    double M[] = {176.39671645238903,
                  -168.85257112196294,
                  -94.877918535603939};
    double gamma[] = {-0.69524430378189617,
                      -0.40429993079272591, 
                      0.59428864887075838}; // Поставил мину

    double step = 0.001;
    double timeSkip = 0;
    double timeSkipClP = 100;
    double calcTime = 200;
    int dimension = 6;
    int expsNum = 6;
    int essDim = 1;
    double eps = 0.001;
    int samplesNum = 1;
    double params[] = {0.485, 2, 6, 7, 9, 4, 1, 752, 100};

    size_t ncu_dots_cnt = (size_t) (calcTime / samplesNum / step);
    size_t traj_dots_cnt = (size_t) ((calcTime / samplesNum + timeSkipClP) / step);
    std::cout<<traj_dots_cnt<<std::endl;
    
    double delta = params[0];
    // 77.3333, -231.617, -94.8779, -0.803553, -0.0335453, 0.594289
    // Q @ M
    mainTrajectory[0] = cos(delta) * M[0] + sin(delta) * M[1];
    mainTrajectory[1] = -sin(delta) * M[0] + cos(delta) * M[1];
    mainTrajectory[2] = M[2];
    // Q @ gamma
    mainTrajectory[3] = cos(delta) * gamma[0] + sin(delta) * gamma[1];
    mainTrajectory[4] = -sin(delta) * gamma[0] + cos(delta) * gamma[1];
    mainTrajectory[5] = gamma[2];

    gamma[0] = mainTrajectory[3];
    gamma[1] = mainTrajectory[4];
    gamma[2] = mainTrajectory[5];

    double* ncu_s = (double*)malloc(ncu_dots_cnt * dimension * essDim * sizeof(double));
    double* traj_dots = (double*)malloc(traj_dots_cnt * dimension * sizeof(double));
    double r[3];
    double omega[3];
    calc_vector_r(r,mainTrajectory,params);
    calc_vector_omega(r,omega,mainTrajectory,params);
    
    tmp[0] = mainTrajectory[0];
    tmp[1] = mainTrajectory[1];
    tmp[2] = mainTrajectory[2];
    tmp[3] = mainTrajectory[3];
    tmp[4] = mainTrajectory[4];
    tmp[5] = mainTrajectory[5];
    
    for (int k = 0; k < 3; k++){
        M[k] = tmp[k] * pow(2 * (params[7] + params[8] * (r[0] * gamma[0] + r[1] * gamma[1] + r[2] * gamma[2])) 
        / (tmp[0] * omega[0] + tmp[1] * omega[1] + tmp[2] * omega[2]), 0.5);
    }
    mainTrajectory[0] = M[0];
    mainTrajectory[1] = M[1];
    mainTrajectory[2] = M[2];

    calc_vector_r(r,mainTrajectory,params);
    calc_vector_omega(r,omega,mainTrajectory,params);
    check_geom_integral(gamma);
    check_energy_integral(M, gamma, r, omega, params[7], params[8]);

    for (int i=0; i<6;i++)
        std::cout << mainTrajectory[i] << std::endl;

    // Пропуск до аттрактора
    for (double t = 0; t < timeSkip-step; t += step)
        dverkStep(mainTrajectory, dimension, diffFunc, params, step, arg, k1, k2, k3, k4, k5, k6, k7, k8);

    // Задаем вектора возмущения
    for (int32_t i = 0; i < expsNum; i++){
        for (int32_t j = 0; j < dimension; j++)
            slaveTrajectories[i * dimension + j] = mainTrajectory[j];
        slaveTrajectories[i * dimension + i] += eps;
    }

    // Прогреваем вектора
    for (double t = 0; t < timeSkipClP-step; t += step) {
        dverkStep(mainTrajectory, dimension, diffFunc, params, step, arg, k1, k2, k3, k4, k5, k6, k7, k8);

        for (int32_t i = 0; i < expsNum; i++) {
            dverkStep(&slaveTrajectories[i * dimension], dimension, diffFunc, params, step, arg, k1, k2, k3, k4, k5, k6, k7, k8);
            
            for (int32_t j = 0; j < dimension; j++)
                slaveTrajectories[i * dimension + j] -= mainTrajectory[j];
        }

        ortVecs(slaveTrajectories, dimension, expsNum, projSum);
        normalizeVecs(slaveTrajectories, dimension, expsNum, eps);

        for (int32_t i = 0; i < expsNum; i++) {
            for (int32_t j = 0; j < dimension; j++)
                slaveTrajectories[i * dimension + j] += mainTrajectory[j];
        }
    }std::cout<<"****forward warm end"<<std::endl;

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

    // Считаем в прямом времени
    int stepCount = 0;
    for (double t = 0; t < calcTime; t += step) {
        dverkStep(mainTrajectory, dimension, diffFunc, params, step, arg, k1, k2, k3, k4, k5, k6, k7, k8);
        // std::cout<<stepCount<<std::endl;
        
        for (int32_t i = 0; i < expsNum; i++){
            dverkStep(&slaveTrajectories[i * dimension], dimension, diffFunc, params, step, arg, k1, k2, k3, k4, k5, k6, k7, k8);
            for (int32_t j = 0; j < dimension; j++)
                slaveTrajectories[i * dimension + j] -= mainTrajectory[j];
        }
        
        ortVecs(slaveTrajectories, dimension, expsNum, projSum);
        for (int32_t i = 0; i < expsNum; i++)
            lyapExp[i] += log(vecNorm(&slaveTrajectories[i * dimension], dimension) / eps);
        normalizeVecs(slaveTrajectories, dimension, expsNum, eps);

        memcpy(ncu_s + stepCount * dimension * essDim, &slaveTrajectories[5*dimension], sizeof(double) * dimension * essDim);
        memcpy(traj_dots + stepCount * dimension, mainTrajectory, sizeof(double) * dimension);
        stepCount++;

        for (int32_t i = 0; i < expsNum; i++) {
            for (int32_t j = 0; j < dimension; j++)
                slaveTrajectories[i * dimension + j] += mainTrajectory[j];
        }
    } std::cout<<"****forward calc end"<<std::endl;

    // Запоминаем точки для прогрева в обратном времени
    for (double t = 0; t < timeSkipClP; t += step) {
        dverkStep(mainTrajectory, dimension, diffFunc, params, step, arg, k1, k2, k3, k4, k5, k6, k7, k8);
        memcpy(traj_dots + stepCount * dimension, mainTrajectory, sizeof(double) * dimension);
        // std::cout<<stepCount<<std::endl;
        // std::cout<<(traj_dots + stepCount * dimension)[0]<<std::endl;
        stepCount++;
    } std::cout<<"****forward remember points end"<<std::endl;

    // Пошли в обратном времени
    // Создадим вектора
    double slaveTrajBw[6*6];
    memset(slaveTrajBw, 0., dimension * expsNum * sizeof(double));

    for (int i = 0; i < expsNum; i++){
        for (int j = 0; j < dimension; j++)
            slaveTrajBw[i * dimension + j] += (traj_dots + (stepCount-1) * dimension)[j];
        slaveTrajBw[i * dimension + i] += eps;
    }

    // Прогрев векторов в обратном времени
    for (double t = 0; t < timeSkipClP; t += step) {
        stepCount--;
        for (int i = 0; i < expsNum; i++) {
            dverkStep(&slaveTrajBw[i * dimension], dimension, diffFunc, params, -step, arg, k1, k2, k3, k4, k5, k6, k7, k8);
            for (int32_t j = 0; j < dimension; j++)
                slaveTrajBw[i * dimension + j] -= (traj_dots + stepCount * dimension)[j];
        }
        
        
        ortVecs(slaveTrajBw, dimension, expsNum, projSum);
        normalizeVecs(slaveTrajBw, dimension, expsNum, eps);

        for (int i = 0; i < expsNum; i++)
            for (int j = 0; j < dimension; j++)
                slaveTrajBw[i * dimension + j] += (traj_dots + stepCount * dimension)[j];

    } std::cout<<"****backtime warm end"<<std::endl;

    double minang = 10;
    double ang = 20;
    for (double t = 0; t < calcTime; t += step) {
        stepCount--;
        for (int i = 0; i < expsNum; i++) {
            dverkStep(&slaveTrajBw[i * dimension], dimension, diffFunc, params, -step, arg, k1, k2, k3, k4, k5, k6, k7, k8);
            for (int32_t j = 0; j < dimension; j++)
                slaveTrajBw[i * dimension + j] -= (traj_dots + stepCount * dimension)[j];
        }
        // std::cout<<t<<std::endl;
        // std::cout<<stepCount<<std::endl;
        // std::cout<<(traj_dots + stepCount * dimension)[0]<<std::endl;
        ortVecs(slaveTrajBw, dimension, expsNum, projSum);

        for (int i = 0; i < expsNum; i++)
            lyapExpBw[i] += log(vecNorm(&slaveTrajBw[i * dimension], dimension) / eps);

        normalizeVecs(slaveTrajBw, dimension, expsNum, eps);
        ang = abs(3.14159265358979323846 / 2 - getAngle(slaveTrajBw, ncu_s + stepCount * dimension, dimension)) ;
        if (ang < minang)
            minang = ang;
        for (int i = 0; i < expsNum; i++)
            for (int j = 0; j < dimension; j++)
                slaveTrajBw[i * dimension + j] += (traj_dots + stepCount * dimension)[j];
    } std::cout<<"****backtime calc end"<<std::endl;

    for (int i=0;i<expsNum;i++)
        std::cout << "L" << i << ": " << lyapExp[i] / calcTime << std::endl;

    for (int i=0;i<expsNum;i++)
        std::cout << "LB" << i << ": " << lyapExpBw[i] / calcTime << std::endl;

    std::cout << "min angle: " << minang << std::endl; 
}