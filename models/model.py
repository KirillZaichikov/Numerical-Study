import numpy as np
import math as m
# from numba import jit

def Shimizu_3D_flow(state, res, params) -> None:
    res[0] = state[1]
    res[1] = - params[0] * state[1] - state[0] * state[2] + state[0]
    res[2] = - params[1] * state[2] + pow( state[0] , 2.0 )

def Lorenz_3D_flow(state, res, params) -> None:
    res[0] = params[0] * (state[1] - state[0])
    res[1] = state[0] * (params[1] - state[2]) - state[1]
    res[2] = state[0] * state[1] - params[2] * state[2]

def Lorenz_3D_flow_var(state, res, params, stateOld):
   res[0] = state[0] * ( - params[0] ) + state[1] * ( params[0] ) 
   res[1] = state[0] * ( params[1] - stateOld[2] ) + state[1] * ( - 1 ) + state[2] * ( - stateOld[0] ) 
   res[2] = state[0] * ( stateOld[1] ) + state[1] * ( stateOld[0] ) + state[2] * ( - params[2] ) 

def Lorenz_4D_flow(state, res, params) -> None: # sigma r b mu 
    res[0] = (params[0] * ( - state[0] + state[1] ) )
    res[1] = (state[0] * ( params[1] - state[2] ) - state[1])
    res[2] = (- params[2] * state[2] + params[3] * state[3] + state[0] * state[1])
    res[3] = (- params[2] * state[3] - params[3] * state[2])

def Lorenz_4D_flow_var(state, res, params, stateOld) -> None: # sigma r b mu 
    res[0] = state[0] * (-params[0]) + state[1] * (params[0]) + state[2] * (0) + state[3] * (0)
    res[1] = state[0] * (params[1] - stateOld[2]) + state[1] * (-1) + state[2] * (-stateOld[0]) + state[3] * (0)
    res[2] = state[0] * (stateOld[1]) + state[1] * (stateOld[0]) + state[2] * (-params[2]) + state[3] * (params[3])
    res[3] = state[0] * (0) + state[1] * (0) + state[2] * (-params[3]) + state[3] * (-params[2])

def Lorenz4D_cross(state):
    return state[2] - 40

def Gonchenko_3D_map(state, res, params) -> None:
    res[0] = state[1] 
    res[1] = ( - params[4] * state[0] + params[3] * state[1] * ( params[0] * state[2] + params[2] ) ) / params[3] 
    res[2] = params[1] * pow ( state[1] , 2.0 ) + params[3] * state[2] 

def ShimizuX3_3D_flow(state, res, params, H) -> None:
    # Param - alpha lamda B 
    res[0] = state[1]
    res[1] = params[2] * pow ( state[0] , 3.0 ) - params[0] * state[1] - state[0] * state[2] + state[0]
    res[2] = - params[1] * state[2] + pow ( state[0] , 2.0 )

def ShimizuX3_3D_flow_var(state, res, params, stateOld):
   res[0] = state[1] *  1  
   res[1] = state[0] * ( 3.0 * params[2] * pow ( stateOld[0] , 2.0 ) - stateOld[2] + 1 ) + state[1] * ( - params[0] ) + state[2] * ( - stateOld[0] ) 
   res[2] = state[0] * ( 2.0 * pow ( stateOld[0] , 1.0 ) ) + state[2] * ( - params[1] ) 

# def ShimizuX3_3D_flow_var_trans(state, res, params, stateOld):
#    res[0] = state[1] * ( 3.0 * params[2] * pow ( stateOld[0] , 2.0 ) - stateOld[2] + 1 ) + state[2] * ( 2.0 * pow ( stateOld[0] , 1.0 ) ) 
#    res[1] = state[0] *  1  + state[1] * ( - params[0] ) 
#    res[2] = state[1] * ( - stateOld[0] ) + state[2] * ( - params[1] ) 


def ShimizuX3_3D_flow_var_rev(state, res, params, stateOld):
   res[0] = -1 * ( state[0] * ( 0 ) +state[1] * ( 1 ) +state[2] * ( 0 )  )
   res[1] = -1 * ( state[0] * ( 3.0 * params[2] * pow ( stateOld[0] , 2.0 ) - stateOld[2] + 1 ) +state[1] * ( - params[0] ) +state[2] * ( - stateOld[0] )  )
   res[2] = -1 * ( state[0] * ( 2.0 * pow ( stateOld[0] , 1.0 ) ) +state[1] * ( 0 ) +state[2] * ( - params[1] )  )

def ShimizuX3_3D_flow_var_trans_rev(state, res, params, stateOld):
   res[0] = -1 * (state[1] * ( 3.0 * params[2] * pow ( stateOld[0] , 2.0 ) - stateOld[2] + 1 ) + state[2] * ( 2.0 * pow ( stateOld[0] , 1.0 ) ) )
   res[1] = -1 * (state[0] *  1  + state[1] * ( - params[0] ) )
   res[2] = -1 * (state[1] * ( - stateOld[0] ) + state[2] * ( - params[1] ) )


def Hw_2D_flow(state, res, params) -> None: # a b
    res[0] = params[0] - (params[1]-1)*state[0] + state[0] ** 2 * state[1]
    res[1] = params[1] * state[0] - state[0] ** 2 * state[1]

def Hw_3D_flow(state, res, params) -> None: # l w m
    res[0] = state[1]
    res[1] = (params[0] + state[2]+state[0]**2-state[0]**4 / 18)*state[1]+params[1]*state[0]
    res[2] = params[2] - state[0] ** 2

# @jit(nopython=True, cache=True)
def lerFRM_3D_map(state, res, params):
    eps, alpha, beta, p = params
    # print(eps, alpha, p_s, p_u)
    ksi0, eta0, teta0 = state

    # LOCAL 1
    # r_factor = (p**2) / np.sqrt(ksi0**2 + eta0**2)
    # scale = (r_factor) ** ((-2 * eps) / (alpha - eps))

    # Phi0 = np.arctan2(eta0, ksi0)
    # # Phi1 = (-(2 * eps) / (alpha - eps)) * np.log(r_factor) + Phi0
    # # if Phi1>np.pi:
    # #     Phi1 -= 2 * np.pi
    # # elif Phi1<-np.pi:
    # #     Phi1 += 2 * np.pi

    # # ksi1 = p_s * p_u * scale * np.cos(Phi1)
    # # eta1 = p_s * p_u * scale * np.sin(Phi1)
    # # phi1 = Phi1 - teta0
    # ksi1 = scale * ( ksi0 * np.cos( ((2*eps)/(alpha-eps)) * np.log(r_factor)) + eta0 * np.sin(((2*eps)/(alpha-eps)) * np.log(r_factor)))
    # eta1 = scale * (-ksi0 * np.sin( ((2*eps)/(alpha-eps)) * np.log(r_factor)) + eta0 * np.cos(((2*eps)/(alpha-eps)) * np.log(r_factor)))
    # # phi1 = (teta0 + Phi0 + ((beta - eps)/(alpha-eps)) * np.log(r_factor)) % (2*np.pi)
    # phi1 =  teta0 + Phi0 + ((beta - eps)/(alpha-eps)) * np.log(r_factor)

    # LOCAL 2
    ro = np.sqrt(ksi0**2+eta0**2)/p
    teta0 = teta0
    phi0 = np.atan2(eta0, ksi0) + teta0

    # while (phi0 > m.pi):
    #     phi0 -= 2 * m.pi
    # while (phi0 < -m.pi):
    #     phi0 += 2 * m.pi

    # print("ro, teta0, phi0")
    # print(ro, teta0, phi0)

    r = p * (ro / p) ** ((alpha+eps)/(alpha-eps))
    teta1 = ((beta + eps) / (alpha - eps)) * np.log(p/ro) + teta0
    phi1 = (((beta - eps) / (alpha - eps)) * np.log(p/ro) + phi0)

    # print(r, teta1, phi1)

    while (phi1 > m.pi):
        phi1 = phi1 - 2 * m.pi
    while (phi1 < -m.pi):
        phi1 += 2 * m.pi

    # while (teta1 > m.pi):
    #     teta1 = teta1 - 2 * m.pi
    # while (teta1 < -m.pi):
    #     teta1 += 2 * m.pi
    # print("r, teta1, phi1")
    # print(r, teta1, phi1)
    ksi1 = r * p * np.cos(phi1-teta1)
    eta1 = r * p * np.sin(phi1-teta1)
    phi1 = phi1

    # print(ksi1, eta1, phi1)
    # LOCAL 3
    # r_factor = (p * p)
    # scale = pow(r_factor, (-2.0 * eps) / (alpha - eps))
    # Phi0 = m.atan2(eta0, ksi0)
    # # angle = (2.0 * eps / (alpha - eps)) * m.log(r_factor)
        
    # ksi1 = scale * ksi0 * (ksi0**2+eta0**2)**(eps/(alpha-eps))
    # eta1 = scale * eta0 * (ksi0**2+eta0**2)**(eps/(alpha-eps))
    # phi1 = teta0 + Phi0 + ((beta - eps) / (alpha - eps)) * m.log(r_factor/(np.sqrt(ksi0**2+eta0**2)))
    # while (phi1 > m.pi):
    #     phi1 = phi1 - 2 * m.pi
    # while (phi1 < -m.pi):
    #     phi1 += 2 * m.pi

    # ******** S3
    ksi2 = ksi1 + eps * np.cos(2*phi1)
    eta2 = eta1 + eps * np.sin(2*phi1)
    # teta2 = (phi1 + 1 + ksi1 + eta1 + eps * np.sin(phi1)) % (2 * np.pi)
    teta2 =  phi1 + 1 + ksi1 + eta1 + eps * np.sin(phi1)
    # teta2 = phi1

    # ******** S2
    # ksi2 = ksi1 + eps * np.cos(phi1)
    # eta2 = eta1 + eps * np.sin(phi1)
    # # # teta2 = (phi1 + 1 + ksi1 + eta1 + eps * np.sin(phi1)) % (2 * np.pi)
    # teta2 =  phi1 + 1 + ksi1 + eta1 + eps * np.sin(phi1)

    # ******** S1
    # ksi2 = ksi1 + eps * (0.5 + 0.25 * np.cos(phi1))
    # eta2 = eta1 + 0.75 * eps * np.sin(phi1)
    # # teta2 = (phi1 + 1 + ksi1 + eta1 + eps * np.sin(phi1)) % (2 * np.pi)
    # teta2 =  phi1 + 1 + ksi1 + eta1 + eps * np.sin(phi1)
    # teta2 =  phi1

    # print(ksi2, eta2, teta2)
    while teta2>np.pi:
        teta2 -= 2 * np.pi
    while teta2<-np.pi:
        teta2 += 2 * np.pi

    res[0] = ksi2 
    res[1] = eta2
    res[2] = teta2

def lerFRM_3D_map_rev(state, res, params):
    eps, alpha, beta, p = params
    p_u = p
    p_s = p
    # print(eps, alpha, p_s, p_u)
    ksi2, eta2, teta2 = state
    ####### Global S3
    # x = 0
    # for i in range(10):
    #     x_0 = x
    #     x = x_0 - (x_0+eps*m.sin(x_0)-eps*m.cos(2*x_0)-eps*m.sin(2*x_0)-teta2+ksi2+eta2+1) / \
    #             (1+eps*m.cos(x_0) + 2*eps*m.sin(2*x_0)-2*eps*m.cos(2*x_0))
    # phi1 = x
    # if (phi1 > m.pi):
    #     phi1 -= 2 * m.pi
    # elif (phi1 < -m.pi):
    #     phi1 += 2 * m.pi
    # ksi1 = ksi2 - eps * m.cos(2*phi1)
    # eta1 = eta2 - eps * m.sin(2*phi1)

    ####### Global S2
    # x = 0
    # for i in range(10):
    #     x_0 = x
    #     x = x_0 - (x_0 - eps * m.cos(x_0) - teta2 + ksi2 + eta2 + 1) / (1 + eps*m.sin(x_0))
    # phi1 = x
    # if (phi1 > m.pi):
    #     phi1 -= 2 * m.pi
    # elif (phi1 < -m.pi):
    #     phi1 += 2 * m.pi
    # ksi1 = ksi2 - eps * m.cos(phi1)
    # eta1 = eta2 - eps * m.sin(phi1)


    ####### Global S1
    x = 0
    for i in range(10):
        x_0 = x
        x = x_0 - (x_0+0.25*eps*m.sin(x_0)-0.25*eps*m.cos(x_0)-teta2+ksi2+eta2-0.5*eps+1) / \
                (1+0.25*eps*m.cos(x_0) + 0.25*eps*m.sin(x_0))
    phi1 = x
    if (phi1 > m.pi):
        phi1 -= 2 * m.pi
    elif (phi1 < -m.pi):
        phi1 += 2 * m.pi
    ksi1 = ksi2 - eps * (0.5 + 0.25 * m.cos(phi1))
    eta1 = eta2 - 0.75 * eps * m.sin(phi1)

    # print("(RevS3) phi1, ksi1, eta1=", phi1, ksi1,eta1)

    ####### Local 1
    eps_ = -eps
    alpha_ = -alpha
    beta_ = -beta

    r = m.sqrt(ksi1**2 + eta1**2) / p_u
    teta1 = phi1 - m.atan2(eta1, ksi1)
    while (teta1 > m.pi):
        teta1 = teta1 - 2 * m.pi
    while (teta1 < -m.pi):
        teta1 += 2 * m.pi
    phi1 = phi1

    ro = p_u * (r / p_s) ** ((alpha_-eps_)/(alpha_+eps_))
    teta0 = ((beta_ + eps_) / (-(alpha_ + eps_))) * m.log(p_s/r) + teta1
    phi0 =  ((beta_ - eps_) / (-(alpha_ + eps_))) * m.log(p_s/r) + phi1

    ksi0 = ro * p_s * m.cos(phi0-teta0)
    eta0 = ro * p_s * m.sin(phi0-teta0)
    phi0 = phi0

    # while (phi0 > m.pi):
    #     phi0 = phi0 - 2 * m.pi
    # while (phi0 < -m.pi):
    #     phi0 += 2 * m.pi

    ##### LOCAL 2
    # ksi0 = (p_s*p_u)**((2*eps)/(alpha+eps))*ksi1*(ksi1**2+eta1**2)**(-eps/(alpha+eps))
    # eta0 = (p_s*p_u)**((2*eps)/(alpha+eps))*eta1*(ksi1**2+eta1**2)**(-eps/(alpha+eps))
    # teta0 = phi1-((beta-eps)/(alpha+eps))*np.log(p_s*p_u/(np.sqrt(ksi1**2+eta1**2)))-np.atan2(eta1,ksi1)

    while (teta0 > m.pi):
        teta0 = teta0 - 2 * m.pi
    while (teta0 < -m.pi):
        teta0 += 2 * m.pi

    res[0] = ksi0
    res[1] = eta0
    res[2] = teta0

def lerFRMwithMu_3D_map(state, res, params):
    mu, eps, alpha, beta, p = params
    # print(eps, alpha, p_s, p_u)
    ksi0, eta0, teta0 = state

    r_factor = (p**2) / np.sqrt(ksi0**2 + eta0**2)
    scale = (r_factor) ** ((-2 * eps) / (alpha - eps))

    Phi0 = np.arctan2(eta0, ksi0)
    # Phi1 = (-(2 * eps) / (alpha - eps)) * np.log(r_factor) + Phi0
    # if Phi1>np.pi:
    #     Phi1 -= 2 * np.pi
    # elif Phi1<-np.pi:
    #     Phi1 += 2 * np.pi

    ksi1 = scale * ( ksi0 * np.cos( ((2*eps)/(alpha-eps)) * np.log(r_factor)) + eta0 * np.sin(((2*eps)/(alpha-eps)) * np.log(r_factor)))
    eta1 = scale * (-ksi0 * np.sin( ((2*eps)/(alpha-eps)) * np.log(r_factor)) + eta0 * np.cos(((2*eps)/(alpha-eps)) * np.log(r_factor)))
    # phi1 = (teta0 + Phi0 + ((beta - eps)/(alpha-eps)) * np.log(r_factor)) % (2*np.pi)
    phi1 =  teta0 + Phi0 + ((beta - eps)/(alpha-eps)) * np.log(r_factor)
    if phi1>np.pi:
        phi1 -= 2 * np.pi
    elif phi1<-np.pi:
        phi1 += 2 * np.pi

    # ******** S1
    # ksi2 = ksi1 + eps * (0.5 + 0.25 * np.cos(phi1))
    # eta2 = eta1 + 0.75 * eps * np.sin(phi1)
    # teta2 = (phi1 + 1 + ksi1 + eta1 + eps * np.sin(phi1)) % (2 * np.pi)
    # teta2 =  phi1 + 1 + ksi1 + eta1 + eps * np.sin(phi1)

    # ******** S2
    ksi2  = ksi1 + mu * np.cos(phi1)
    eta2  = eta1 + mu * np.sin(phi1)
    teta2 = (phi1 + 1 + ksi1 + eta1 + mu * np.sin(phi1)) % (2 * np.pi) - np.pi
    # teta2 = (phi1 + 1 + ksi1 + eta1 + eps * np.sin(phi1)) % (2 * np.pi)
    # teta2 =  phi1 + 1 + ksi1 + eta1 + mu * np.sin(phi1)
    # if teta2>np.pi:
    #     teta2 -= 2 * np.pi
    # elif teta2<-np.pi:
    #     teta2 += 2 * np.pi

    res[0] = ksi2 
    res[1] = eta2
    res[2] = teta2

def lerFRMwithMuRev_3D_map(state, res, params):
    mu, eps, alpha, beta, p = params
    # print(eps, alpha, p_s, p_u)
    ksi2, eta2, teta2 = state

    # ******** S2
    ksi1 = ksi2 + mu * np.cos(phi1)
    eta1 = eta2 + mu * np.sin(phi1)
    phi1 =  teta2 - 1 - ksi1 - eta1 - mu * np.sin(phi1)
    if phi1>np.pi:
        phi1 -= 2 * np.pi
    elif phi1<-np.pi:
        phi1 += 2 * np.pi

    r_factor = (p**2) / np.sqrt(ksi0**2 + eta0**2)
    scale = (r_factor) ** ((-2 * eps) / (alpha - eps))

    Phi0 = np.arctan2(eta0, ksi0)
    # Phi1 = (-(2 * eps) / (alpha - eps)) * np.log(r_factor) + Phi0
    # if Phi1>np.pi:
    #     Phi1 -= 2 * np.pi
    # elif Phi1<-np.pi:
    #     Phi1 += 2 * np.pi

    ksi1 = scale * ( ksi0 * np.cos( ((2*eps)/(alpha-eps)) * np.log(r_factor)) + eta0 * np.sin(((2*eps)/(alpha-eps)) * np.log(r_factor)))
    eta1 = scale * (-ksi0 * np.sin( ((2*eps)/(alpha-eps)) * np.log(r_factor)) + eta0 * np.cos(((2*eps)/(alpha-eps)) * np.log(r_factor)))
    # phi1 = (teta0 + Phi0 + ((beta - eps)/(alpha-eps)) * np.log(r_factor)) % (2*np.pi)
    phi1 =  teta0 + Phi0 + ((beta - eps)/(alpha-eps)) * np.log(r_factor)
    if phi1>np.pi:
        phi1 -= 2 * np.pi
    elif phi1<-np.pi:
        phi1 += 2 * np.pi

    # ******** S1
    # ksi2 = ksi1 + eps * (0.5 + 0.25 * np.cos(phi1))
    # eta2 = eta1 + 0.75 * eps * np.sin(phi1)
    # teta2 = (phi1 + 1 + ksi1 + eta1 + eps * np.sin(phi1)) % (2 * np.pi)
    # teta2 =  phi1 + 1 + ksi1 + eta1 + eps * np.sin(phi1)



    res[0] = ksi2 
    res[1] = eta2
    res[2] = teta2

def smaleWilliams_3D_map(state,res,params):
    phi, X, Y = state
    alpha = params[0]

    phi_ = 2 * phi
    # if (phi_ > m.pi):
    #     phi_ = phi_ - 2 * m.pi
    # elif (phi_ < -m.pi):
    #     phi_ += 2 * m.pi
    if (phi_ > 2*m.pi):
        phi_ = phi_ - 2 * m.pi
    elif (phi_ < 0):
        phi_ += 2 * m.pi
    X_ = alpha * X + (1/2)*np.cos(phi)
    Y_ = alpha * Y + (1/2)*np.sin(phi)

    res[0] = phi_
    res[1] = X_
    res[2] = Y_

def smaleWilliams_3D_map_rev(state,res,params):
    phi_, X_, Y_ = state
    alpha, phiOld = params

    phi = phi_ / 2
    if abs(phiOld - phi) > 0.2:
        phi = phi + np.pi
    # phi = phi_ / 2
    # if (phi_ > m.pi):
    #     phi_ = phi_ - 2 * m.pi
    # elif (phi_ < -m.pi):
    #     phi_ += 2 * m.pi
    # if (phi > 2*m.pi):
    #     phi = phi - 2 * m.pi
    # elif (phi < 0):
    #     phi += 2 * m.pi
    X = (X_ - (1/2)*np.cos(phi)) / alpha
    Y = (Y_ - (1/2)*np.sin(phi)) / alpha

    res[0] = phi
    res[1] = X
    res[2] = Y

def test_model_3D_map(state,res,params):
    eps, u, v = params

    x = - state[1]
    y = state[0] + 3*state[1] + 2*state[2]
    z = 2*state[1] + state[2]

    x = norm(x)
    y = norm(y)
    z = norm(z)

    x_ = u + (1/np.pi) * np.atan2((1-eps) * np.sin(np.pi*(x-v)),(1+eps) * np.cos(np.pi*(x-v)))
    y_ = u + (1/np.pi) * np.atan2((1-eps) * np.sin(np.pi*(y-v)),(1+eps) * np.cos(np.pi*(y-v)))
    z_ = u + (1/np.pi) * np.atan2((1-eps) * np.sin(np.pi*(z-v)),(1+eps) * np.cos(np.pi*(z-v)))

    res[0] = x_
    res[1] = y_
    res[2] = z_

def test_model_3D_map_rev(state,res,params):
    eps, u, v = params
    eps = -eps

    x_ = v + (1/np.pi) * np.atan2((1-eps) * np.sin(np.pi*(state[0]-u)),(1+eps) * np.cos(np.pi*(state[0]-u)))
    y_ = v + (1/np.pi) * np.atan2((1-eps) * np.sin(np.pi*(state[1]-u)),(1+eps) * np.cos(np.pi*(state[1]-u)))
    z_ = v + (1/np.pi) * np.atan2((1-eps) * np.sin(np.pi*(state[2]-u)),(1+eps) * np.cos(np.pi*(state[2]-u)))
    
    x = y_ - x_ - 2 * z_
    y = - x_
    z = z_ + 2 * x_

    x = norm(x)
    y = norm(y)
    z = norm(z)

    res[0] = x
    res[1] = y
    res[2] = z

def test_model_3D_map_var(state,res,params):
    eps, u, v = params
    eps = -eps

    x_ = v + (1/np.pi) * np.atan2((1-eps) * np.sin(np.pi*(state[0]-u)),(1+eps) * np.cos(np.pi*(state[0]-u)))
    y_ = v + (1/np.pi) * np.atan2((1-eps) * np.sin(np.pi*(state[1]-u)),(1+eps) * np.cos(np.pi*(state[1]-u)))
    z_ = v + (1/np.pi) * np.atan2((1-eps) * np.sin(np.pi*(state[2]-u)),(1+eps) * np.cos(np.pi*(state[2]-u)))
    
    x = y_ - x_ - 2 * z_
    y = - x_
    z = z_ + 2 * x_

    x = norm(x)
    y = norm(y)
    z = norm(z)

    res[0] = x
    res[1] = y
    res[2] = z

def norm(x):
    while x>(1/2):
        x = x-1
    while x<=(-1/2):
        x = x+1
    return x


