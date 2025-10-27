#include <cmath>
#include "linal.cpp"

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

    // Нормировка gamma
    // double gamma_norm = sqrt(gamma[0]*gamma[0]+gamma[1]*gamma[1]+gamma[2]*gamma[2]);
    // for (int i=0;i<3;i++)
    //     gamma[i] = gamma[i] / gamma_norm;

    // Вычисление вектора r
    r[0] = -params[4] * gamma[0] / gamma[2];
    r[1] = -params[5] * gamma[1] / gamma[2];
    r[2] = -params[6] + (params[4] * pow(gamma[0], 0.2e1) + params[5] * pow(gamma[1], 0.2e1)) * pow(gamma[2], -0.2e1) / 0.2e1;
    // std::cout << "vector r " << r[0] << " " << r[1] << " " << r[2] << std::endl;

    // r1 = r[0];
    // r2 = r[1];
    // r3 = r[2];
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

    // JQ_kuz[0][0] = pow(cos(delta), 0.2e1) * J[0] + pow(sin(delta), 0.2e1) * J[1] + r2 * r2 + r3 * r3;
    // JQ_kuz[0][1] = -cos(delta) * J[0] * sin(delta) + sin(delta) * J[1] * cos(delta) - r1 * r2;
    // JQ_kuz[0][2] = -r1 * r3;
    // JQ_kuz[1][0] = -cos(delta) * J[0] * sin(delta) + sin(delta) * J[1] * cos(delta) - r1 * r2;
    // JQ_kuz[1][1] = pow(sin(delta), 0.2e1) * J[0] + pow(cos(delta), 0.2e1) * J[1] + r1 * r1 + r3 * r3;
    // JQ_kuz[1][2] = -r2 * r3;
    // JQ_kuz[2][0] = -r1 * r3;
    // JQ_kuz[2][1] = -r2 * r3;
    // JQ_kuz[2][2] = r1 * r1 + r2 * r2 + J[2];

    // Обратная матрица М от омега
    inverse3x3(JQ_kuz, JQ_rev);

    // Посчитали вектор omega
    omega[0] = JQ_rev[0][0] * M[0] + JQ_rev[0][1] * M[1] + JQ_rev[0][2] * M[2];
    omega[1] = JQ_rev[1][0] * M[0] + JQ_rev[1][1] * M[1] + JQ_rev[1][2] * M[2];
    omega[2] = JQ_rev[2][0] * M[0] + JQ_rev[2][1] * M[1] + JQ_rev[2][2] * M[2];
    // std::cout << "vector omega " << omega[0] << " " << omega[1] << " " << omega[2] << std::endl;


    // НОРМИРОВКА ЭНЕРГИИ
    // tmp[0] = M[0];
    // tmp[1] = M[1];
    // tmp[2] = M[2];
    // tmp[3] = gamma[0];
    // tmp[4] = gamma[1];
    // tmp[5] = gamma[2];
    // for (int k = 0; k < 3; k++)
    //     M[k] = tmp[k] * pow((2 * (params[7] + params[8] * (r[0] * gamma[0] + r[1] * gamma[1] + r[2] * gamma[2])) / (tmp[0] * omega[0] + tmp[1] * omega[1] + tmp[2] * omega[2])), 0.5);


    // Посчитали вектор dgamma
    dgamma[0] =  gamma[1] * omega[2] - gamma[2] * omega[1];
    dgamma[1] = -gamma[0] * omega[2] + gamma[2] * omega[0];
    dgamma[2] =  gamma[0] * omega[1] - gamma[1] * omega[0];
    // std::cout << "vector dgamma " << dgamma[0] << " " << dgamma[1] << " " << dgamma[2] << std::endl;

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
    // std::cout << "vector dr " << dr[0] << " " << dr[1] << " " << dr[2] << std::endl;

    // Посчитали dM
    dM[0] =  M[1] * omega[2] - M[2] * omega[1] + dr[1] * ( omega[0] * r[1] - omega[1] * r[0]) - dr[2] * (-omega[0] * r[2] + omega[2] * r[0]) + g0 * (-r[2] * gamma[1] + r[1] * gamma[2]);
    dM[1] = -M[0] * omega[2] + M[2] * omega[0] - dr[0] * ( omega[0] * r[1] - omega[1] * r[0]) + dr[2] * ( omega[1] * r[2] - omega[2] * r[1]) + g0 * ( r[2] * gamma[0] - r[0] * gamma[2]);
    dM[2] =  M[0] * omega[1] - M[1] * omega[0] + dr[0] * (-omega[0] * r[2] + omega[2] * r[0]) - dr[1] * ( omega[1] * r[2] - omega[2] * r[1]) + g0 * (-r[1] * gamma[0] + r[0] * gamma[1]);
    // std::cout << "vector dM " << dM[0] << " " << dM[1] << " " << dM[2] << std::endl;

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
    double H;

    // Начальные условия и фазовые
    double M[3], gamma[3];
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
    // std::cout << "vector r " << r[0] << " " << r[1] << " " << r[2] << std::endl;

    r1 = r[0];
    r2 = r[1];
    r3 = r[2];
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

    // JQ_kuz[0][0] = pow(cos(delta), 0.2e1) * J[0] + pow(sin(delta), 0.2e1) * J[1] + r2 * r2 + r3 * r3;
    // JQ_kuz[0][1] = -cos(delta) * J[0] * sin(delta) + sin(delta) * J[1] * cos(delta) - r1 * r2;
    // JQ_kuz[0][2] = -r1 * r3;
    // JQ_kuz[1][0] = -cos(delta) * J[0] * sin(delta) + sin(delta) * J[1] * cos(delta) - r1 * r2;
    // JQ_kuz[1][1] = pow(sin(delta), 0.2e1) * J[0] + pow(cos(delta), 0.2e1) * J[1] + r1 * r1 + r3 * r3;
    // JQ_kuz[1][2] = -r2 * r3;
    // JQ_kuz[2][0] = -r1 * r3;
    // JQ_kuz[2][1] = -r2 * r3;
    // JQ_kuz[2][2] = r1 * r1 + r2 * r2 + J[2];

    // Обратная матрица М от омега
    inverse3x3(JQ_kuz, JQ_rev);

    // Посчитали вектор omega
    omega[0] = JQ_rev[0][0] * M[0] + JQ_rev[0][1] * M[1] + JQ_rev[0][2] * M[2];
    omega[1] = JQ_rev[1][0] * M[0] + JQ_rev[1][1] * M[1] + JQ_rev[1][2] * M[2];
    omega[2] = JQ_rev[2][0] * M[0] + JQ_rev[2][1] * M[1] + JQ_rev[2][2] * M[2];
    // std::cout << "vector omega " << omega[0] << " " << omega[1] << " " << omega[2] << std::endl;

    // Посчитали вектор dgamma
    dgamma[0] =  gamma[1] * omega[2] - gamma[2] * omega[1];
    dgamma[1] = -gamma[0] * omega[2] + gamma[2] * omega[0];
    dgamma[2] =  gamma[0] * omega[1] - gamma[1] * omega[0];
    // std::cout << "vector dgamma " << dgamma[0] << " " << dgamma[1] << " " << dgamma[2] << std::endl;

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
    // std::cout << "vector dr " << dr[0] << " " << dr[1] << " " << dr[2] << std::endl;

    // Посчитали dM
    dM[0] =  M[1] * omega[2] - M[2] * omega[1] + dr[1] * ( omega[0] * r[1] - omega[1] * r[0]) - dr[2] * (-omega[0] * r[2] + omega[2] * r[0]) + g0 * (-r[2] * gamma[1] + r[1] * gamma[2]);
    dM[1] = -M[0] * omega[2] + M[2] * omega[0] - dr[0] * ( omega[0] * r[1] - omega[1] * r[0]) + dr[2] * ( omega[1] * r[2] - omega[2] * r[1]) + g0 * ( r[2] * gamma[0] - r[0] * gamma[2]);
    dM[2] =  M[0] * omega[1] - M[1] * omega[0] + dr[0] * (-omega[0] * r[2] + omega[2] * r[0]) - dr[1] * ( omega[1] * r[2] - omega[2] * r[1]) + g0 * (-r[1] * gamma[0] + r[0] * gamma[1]);
    // std::cout << "vector dM " << dM[0] << " " << dM[1] << " " << dM[2] << std::endl;

    H = dM[1] * gamma[0] + M[1] * dgamma[0] - dM[0] * gamma[1] - M[0] * dgamma[1];//M2' * gamma1 + M2 * gamma1' - M1' * gamma2 - M1 * gamma2'

    res[0] = dM[0] / H;
    res[1] = dM[1] / H;
    res[2] = dM[2] / H;
    res[3] = dgamma[0] / H; 
    res[4] = dgamma[1] / H;
    res[5] = dgamma[2] / H;
}