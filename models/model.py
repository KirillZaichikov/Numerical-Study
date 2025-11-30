import numpy as np

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

def lerFRM_3D_map(state, res, params):
    eps, alpha, beta, p = params
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

    # ksi1 = p_s * p_u * scale * np.cos(Phi1)
    # eta1 = p_s * p_u * scale * np.sin(Phi1)
    # phi1 = Phi1 - teta0
    ksi1 = scale * ( ksi0 * np.cos( ((2*eps)/(alpha-eps)) * np.log(r_factor)) + eta0 * np.sin(((2*eps)/(alpha-eps)) * np.log(r_factor)))
    eta1 = scale * (-ksi0 * np.sin( ((2*eps)/(alpha-eps)) * np.log(r_factor)) + eta0 * np.cos(((2*eps)/(alpha-eps)) * np.log(r_factor)))
    # phi1 = (teta0 + Phi0 + ((beta - eps)/(alpha-eps)) * np.log(r_factor)) % (2*np.pi)
    phi1 =  teta0 + Phi0 + ((beta - eps)/(alpha-eps)) * np.log(r_factor)
    if phi1>np.pi:
        phi1 -= 2 * np.pi
    elif phi1<-np.pi:
        phi1 += 2 * np.pi

    # ******** S3
    ksi2 = ksi1 + eps * np.cos(2 * phi1)
    eta2 = eta1 + eps * np.sin(2 * phi1)
    # teta2 = (phi1 + 1 + ksi1 + eta1 + eps * np.sin(phi1)) % (2 * np.pi)
    teta2 =  phi1 + 1 + ksi1 + eta1 + eps * np.sin(phi1)

    # ******** S1
    # ksi2 = ksi1 + eps * (0.5 + 0.25 * np.cos(phi1))
    # eta2 = eta1 + 0.75 * eps * np.sin(phi1)
    # teta2 = (phi1 + 1 + ksi1 + eta1 + eps * np.sin(phi1)) % (2 * np.pi)
    # teta2 =  phi1 + 1 + ksi1 + eta1 + eps * np.sin(phi1)
    if teta2>np.pi:
        teta2 -= 2 * np.pi
    elif teta2<-np.pi:
        teta2 += 2 * np.pi

    res[0] = ksi2 
    res[1] = eta2
    res[2] = teta2

def lerFRM_3D_map(state, res, params):
    eps, alpha, beta, p = params
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

    # ksi1 = p_s * p_u * scale * np.cos(Phi1)
    # eta1 = p_s * p_u * scale * np.sin(Phi1)
    # phi1 = Phi1 - teta0
    ksi1 = scale * ( ksi0 * np.cos( ((2*eps)/(alpha-eps)) * np.log(r_factor)) + eta0 * np.sin(((2*eps)/(alpha-eps)) * np.log(r_factor)))
    eta1 = scale * (-ksi0 * np.sin( ((2*eps)/(alpha-eps)) * np.log(r_factor)) + eta0 * np.cos(((2*eps)/(alpha-eps)) * np.log(r_factor)))
    # phi1 = (teta0 + Phi0 + ((beta - eps)/(alpha-eps)) * np.log(r_factor)) % (2*np.pi)
    phi1 =  teta0 + Phi0 + ((beta - eps)/(alpha-eps)) * np.log(r_factor)
    if phi1>np.pi:
        phi1 -= 2 * np.pi
    elif phi1<-np.pi:
        phi1 += 2 * np.pi

    # ******** S3
    ksi2 = ksi1 + eps * np.cos(phi1)
    eta2 = eta1 + eps * np.sin(phi1)
    # teta2 = (phi1 + 1 + ksi1 + eta1 + eps * np.sin(phi1)) % (2 * np.pi)
    teta2 =  phi1 + 1 + ksi1 + eta1 + eps * np.sin(phi1)

    # ******** S1
    # ksi2 = ksi1 + eps * (0.5 + 0.25 * np.cos(phi1))
    # eta2 = eta1 + 0.75 * eps * np.sin(phi1)
    # teta2 = (phi1 + 1 + ksi1 + eta1 + eps * np.sin(phi1)) % (2 * np.pi)
    # teta2 =  phi1 + 1 + ksi1 + eta1 + eps * np.sin(phi1)
    if teta2>np.pi:
        teta2 -= 2 * np.pi
    elif teta2<-np.pi:
        teta2 += 2 * np.pi

    res[0] = ksi2 
    res[1] = eta2
    res[2] = teta2

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




