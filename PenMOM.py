from Problem import info
import numpy as np
def bracket(X):
    return min(0,X)

def penalty_mom(prob, sigma, tau,R,X):
    func, inequality_constraints, equality_constraints = info(prob)
    errG = []
    errH = []
    
    for i in range(len(inequality_constraints)):
        errG.append(bracket(bracket(inequality_constraints[i](X))+sigma[i]))
    
    for i in range(len(equality_constraints)):
        errH.append(equality_constraints[i](X)+tau[i])

    errG = np.array(errG)
    errH = np.array(errH)
    #print('errG and sigma',errG,sigma)
    
    cons_voil = 0
    if len(errG)>0:
        cons_voil += (np.dot(errG,errG) - np.dot(sigma,sigma))
    if len(errH)>0:
        cons_voil += (np.dot(errH,errH) - np.dot(tau,tau))
    #print(X,func(X)+cons_voil,cons_voil)
    return [func(X)+R*cons_voil,cons_voil]
