def Shimizu_3D_flow(state, res, params) -> None:
    res[0] = state[1]
    res[1] = - params[0] * state[1] - state[0] * state[2] + state[0]
    res[2] = - params[1] * state[2] + pow( state[0] , 2.0 )

def Lorenz_3D_flow(state, res, params) -> None:
    res[0] = params[0] * (state[1] - state[0])
    res[1] = state[0] * (params[1] - state[2]) - state[1]
    res[2] = state[0] * state[1] - params[2] * state[2]

def Lorenz_4D_flow(state, res, params) -> None: # sigma r b mu 
    res[0] = params[0] * ( - state[0] + state[1] ) 
    res[1] = state[0] * ( params[1] - state[2] ) - state[1] 
    res[2] = - params[2] * state[2] + params[3] * state[3] + state[0] * state[1] 
    res[3] = - params[2] * state[3] - params[3] * state[2]

def Lorenz_4D_flow_var(state, res, params, stateOld) -> None: # sigma r b mu 
    res[0] = state[0] * (-params[0]) + state[1] * (params[0]) + state[2] * (0) + state[3] * (0)
    res[1] = state[0] * (params[1] - stateOld[2]) + state[1] * (-1) + state[2] * (-stateOld[0]) + state[3] * (0)
    res[2] = state[0] * (stateOld[1]) + state[1] * (stateOld[0]) + state[2] * (-params[2]) + state[3] * (params[3])
    res[3] = state[0] * (0) + state[1] * (0) + state[2] * (-params[3]) + state[3] * (-params[2])

def Gonchenko_3D_map(state, res, params) -> None:
    res[0] = state[1] 
    res[1] = ( - params[4] * state[0] + params[3] * state[1] * ( params[0] * state[2] + params[2] ) ) / params[3] 
    res[2] = params[1] * pow ( state[1] , 2.0 ) + params[3] * state[2] 

def ShimizuX3_3D_flow(state, res, params) -> None:
    # Param - alpha lamda B 
    res[0] = state[1]
    res[1] = params[2] * pow ( state[0] , 3.0 ) - params[0] * state[1] - state[0] * state[2] + state[0]
    res[2] = - params[1] * state[2] + pow ( state[0] , 2.0 )
