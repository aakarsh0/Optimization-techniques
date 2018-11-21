import math
import numpy as np
import matplotlib.pyplot as plt
from PenMOM import penalty_mom, bracket
from Problem import info

## WRITE THE PROBLEM NO., DEFAULT IS ZERO(HAMMELBLAU FUNCTION)
prob = 0

print('Problem_no',prob)
print()
## list of initial points
starting_points2 = [[1,4],[-5,16],[-1,-1],[11,-3],[-5,-5],[0.1,-0.3]]
starting_points5 = [[1,1,1,1,1],[0.5,2,-1,3,-2],[4,5,-2,2,-1]]
#starting_points = [[1,2,3,4,5,6,7,8],[600,1500,5000,300,400,500,800,700]]
#starting_points = [[1,1,1],[0.5,2,-1]]
##tolerance
e = 0.000001
R = 1

def func(X):
    global sigma
    global tau
    global R
    global prob
    return penalty_mom(prob,sigma,tau,R,X)[0]

def dotproduct(directions,newdir):
    for d in directions:
        if 1-np.dot(d,newdir)<1e-3:
            return 1
    return 0		

def new_s(X2,X1): ##gives unit vector and dist in direction of X2-X1
    #X2 = list(X2)
    #X1 = list(X1)
    #s = np.array([],dtype = 'float64')
    s = X2 - X1
    dist = np.float64(0)
    for i in range(len(s)):
        #np.append(s,X2[i]-X1[i])
        dist += (s[i])**2
    dist = dist**(0.5)
    s = s/dist
    return [s,dist]


## Bracketing
def bracketing1(X,s):
    delta=np.float64(0.01)
    X = np.array(X,dtype='float64')
    s = np.array(s,dtype='float64')
    fmin = func(X - delta*s)
    fplus = func(X + delta*s)
    fx = func(X)
    f_eval = 3
    if fplus>fx and fmin>fx:
        return [X-delta*s, X+delta*s,f_eval]
    if fplus>fx>fmin:
        delta = -1*delta
    while True:
        #alpha = delta
        delta = delta*2
        #fx = func(X+alpha*s)
        fplus = func(X+s*delta)
        f_eval += 1
        if fplus>fx :
            #print(X+alpha*s)
            #print(X+delta*s)
            #print('bracketing_done')
            return [X, X+delta*s, f_eval]
        
def bracketing(X,s):
    delta=0.01
    X = np.array(X)
    s = np.array(s)
    fmin = func(X - delta*s)
    fplus = func(X + delta*s)
    fx = func(X)
    f_eval = 3
    if fplus>fx and fmin>fx:
        return [X-delta*s, X+delta*s,f_eval]
    if fplus>fx>fmin :
        delta = -1*delta
    Xless=X
    while True:
        delta = delta*2
        Xplus = X + s*delta
        fx = fplus
        fplus = func(Xplus)
        f_eval += 1
        if fplus>fx :
            return [Xless, Xplus,f_eval]
        Xless= X
        X= Xplus

##golden section method
def golden(low,high):
    global e
    dist = new_s(low,high)[1]
    tau = np.float64(1-(pow(5,0.5)-1)/2)
    eps = 1e-12
    a = np.array(low,dtype='float64')
    b = np.array(high,dtype='float64')
    x1 = a+(b-a)*tau
    x2 = b+(a-b)*tau
    fx1 = func(x1)
    fx2 = func(x2)
    f_eval = 2
    while dist>eps:
        f_eval += 1
        if fx1>fx2:
            a = x1
            x1 = x2
            fx1 = fx2
            x2 = b + (a-b)*tau
            fx2 = func(x2)
        else:
            b = x2
            x2 = x1
            fx2 = fx1
            x1 = a +(b-a)*tau
            fx1 = func(x1)
        dist = dist*tau
    #print('golden_section_done')
    return [(x1+x2)/2, f_eval]
            
##this function calls bracketing and golden section functions
def find_minima(X,s):
    low,high, func_eval_bracket = bracketing(X,s)
    minima, func_eval_golden = golden(low,high)
    return [minima, func_eval_bracket + func_eval_golden]

##Conjugate Direction Method
def conjugate_direction(X):
    global e
    global sigma
    global tau
    global prob
    global f_calc
    global fx_now
    f_eval = 0
    n = len(X)
    X = np.array(X,dtype = 'float64')
    ## creating list of conjugate directions
    directions = [[np.float64(0) for i in range(n)] for j in range(n)]
    for i in range(n):
        directions[i][i] = np.float64(1)
    ## Iterations to find next point and direction
    iterations = 0
    while iterations<450:
        minii = find_minima(X,np.array(directions[0]))
        X_new = minii[0]
        f_eval = minii[1]
        X1 = X_new
        for s in directions[1:]:
            minii = find_minima(X_new,np.array(s))
            X_new = minii[0]
            f_eval += minii[1]
        minii = find_minima(X_new,np.array(directions[0]))
        X_new = minii[0]
        f_eval = minii[1]
        if abs(new_s(X_new,X1)[1])<=e:
            #print('break',X_new,X1)
            break
        if dotproduct(directions,new_s(X_new,X1)[0]):
            return conjugate_direction(X_new)
        f_calc.append(f_calc[-1]+f_eval)
        fx_now.append(func(list(X_new)))
        X = X_new
        #print(X,func(X))
        directions.pop(0)
        directions.append(new_s(X_new,X1)[0])
        #print(directions)
        iterations+=1
    #print(X_new,'minima',penalty(X_new,0)[0],'penality',penalty(X_new,0)[1])
    
    return X_new


def solveMOM(X):
    global R
    global prob
    global sigma
    global tau
    global fx_now
    global f_calc
    func, inequality_constraints, equality_constraints = info(prob)
    e1=1e-5
    # e2=1e-5
    R=0.1
    MAX_ITR = 800
    pen_prev= 1e9+0.1
    tau = np.zeros(len(equality_constraints))
    sigma = np.zeros(len(inequality_constraints))
    fx_now = [penalty_mom(prob, sigma, tau, R, X)[0]]
    f_calc = [0]
    it=0
    X1 = X
    #print('e1',e1)
    while True:
        it+=1
        X_now = np.array(conjugate_direction(X1))
        pen_now,cons_voil = penalty_mom(prob,sigma,tau,R,X_now)
        if(abs(pen_now-pen_prev)<e1) and cons_voil<.0001:
            #print('voilet',pen_now,pen_prev)
            return X_now
        if it>=MAX_ITR:
            print('terminated_after_max_iterations ',abs(pen_now-pen_prev))
            return X_now
        #print('pen_prev:', pen_prev,'cons_voil:',cons_voil)
##        print('pen_now:', pen_now)
        pen_prev = pen_now
##        sigma = bracket(inequality_constraints(X_now)+sigma)
##        tau = equality_constraints(X_now)+tau
        sigma1 = np.array([bracket(inequality_constraints[i](X_now)+sigma[i]) for i in range(len(sigma))])
        tau1 = np.array([(equality_constraints[i](X_now)+tau[i]) for i in range(len(tau))])
        sigma = sigma1
        tau = tau1
        X1 = X_now

sigma = []
tau = []
fx_now = []
f_calc = [0]

if prob==5:
    starting_points = starting_points5
else:
    starting_points = starting_points2
for s in starting_points:
    mini = solveMOM(s)
    #print(f_calc,fx_now)
    print('Start_Point',s)
    print('minima:',mini)
    print('f(X)= ',penalty_mom(prob,sigma,tau,R,mini)[0]-R*penalty_mom(prob,sigma,tau,R,mini)[1])
    print('constraint_voilation= ',penalty_mom(prob,sigma,tau,R,mini)[1])
    plt.plot(f_calc,fx_now)
    plt.xlabel('Function Evaluations')
    plt.ylabel('Penalty')
    plt.show()
    
    print()
    print()
       
    



