import math
import numpy as np
import matplotlib.pyplot as plt
from Problem import info
#from Penalty import penalty, set_prob
import matplotlib.pyplot as plt
#from input import starts

## WRITE THE PROBLEM NO. HERE FOR VARIABLE NAME, 'prob'
## Default problem is Hammelblau function otherwise.
prob = 0

obj, inequality_constraints, equality_constraints = info(prob) ##getting the function and contraints equations

def bracket(value):
    return min(0,value) ## check positive or negative

def penalty(X,R): ##calculate and return the penalty
    global obj
    global inequality_constraints
    global equality_constraints
    objective = obj(X)
    extra = 0
    for eqn in inequality_constraints:
        extra += bracket(eqn(X))**2
    for eqn in equality_constraints:
        extra += eqn(X)**2
    # returns the total penalty and penalty due to constraints
    return [objective + R*extra, extra]

## list of initial points
starting_points2 = [[0,0],[1,4],[-5,16],[11,-3],[-5,-5],[0.1,-0.3],[14,0.8]] #0,2,3
starting_points5 = [[1,1,1,1,1],[0.5,2,-1,3,-2],[2,4,7,-2,5]] #1,5

if prob==5:
    starting_points = starting_points5
else:
    starting_points = starting_points2
    
##tolerance
e = 0.000001
c = 2
R = 0.0001
x_current = [] ## used to store points to plot curves if needed
f_calc = [0]


def func(X):## calls penalty function
    global R
    #X = list(X)
    return penalty(X,R)[0]
    ## it calls the penalty function to get (f(x)+penalty)		

def new_s(X2,X1): ##gives unit vector and dist in direction of X2-X1
    #X2 = list(X2)
    #X1 = list(X1)
    #s = np.array([],dtype = 'float64')
    s = X2-X1 ##vector in direction of X1 to X2
    dist = np.float64(0) ## to calculate distance between the two points
    for i in range(len(s)):
        #np.append(s,X2[i]-X1[i])
        try:
            dist += (s[i])**2
        except:
            pass
    if dist<1e-7: #
        return [s,0]
    dist = dist**(0.5)
    s = s/dist
    return [s,dist] #return unit vector and distence between the two points


## Bracketing
def bracketing(X,s): ##Bounding phase method, takes a point and direction
    delta = np.float64(0.01) #start with delta = 0.01
    X = np.array(X)   #convert to numpy array for faster calculations
    s = np.array(s)
    fmin = func(X - delta*s) #f(x) value at point slightly left in the given direction
    fplus = func(X + delta*s) #f(x) at value at point slightle right
    fx = func(X) #f(x) at x
    f_eval = 3 #to count total function evaluations
    if fplus>fx and fmin>fx:
        return [X-delta*s, X+delta*s,f_eval]
    if fplus>fx>fmin : #decideing the direction of minima from X
        delta = -1 * delta
    Xless = X
    while True: ##finding bounds
        delta = delta*2
        #print(fplus,fx,'delta= ',delta, 'X= ',X)
        try:
            Xplus = X + delta*s
        except:
            print('s= ',s,'delta= ',delta)
        fx = fplus
        fplus = func(Xplus)
        f_eval += 1
        if fplus>=fx :
            return [Xless, Xplus,f_eval] #return lower and upper bounds and function evaluation count
        Xless= X
        X= Xplus

##golden section method
def golden(low,high):
    dist = new_s(low,high)[1] #get the distance between bounds
    tau = np.float64(1-(pow(5,0.5)-1)/2) #declaring tau, golden ratio
    eps = 1e-12 ## tolerance
    a = np.array(low,dtype='float64') #converting to numpy array
    b = np.array(high,dtype='float64')
    x1 = a+(b-a)*tau #getting the two points between bounds using tau
    x2 = b+(a-b)*tau
    fx1 = func(x1) #calculating function value at new points
    fx2 = func(x2)
    f_eval = 2 #function evaluation counter
    while dist>eps: #termination condition
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
    return [(x1+x2)/2, f_eval] #return minima point and function values
            

def find_minima(X,s): ##this function calls bracketing and golden section functions
    low, high, func_eval_bracket = bracketing(X,s)
    minima, func_eval_golden = golden(low,high)
    return [minima, func_eval_bracket + func_eval_golden]

def dotproduct(directions,newdir): ##gives dot product of two vectors
    for d in directions:
        if 1 - np.dot(d,newdir)<1e-3:
            return 1
    return 0


##Conjugate Direction Method
def conjugate_direction(X):
    global e
    global fx_current
    global f_calc
    f_eval = f_calc[-1]
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
        f_eval += minii[1]
        fx_current.append(func(X_new))
        f_calc.append(f_eval)
        X1 = X_new
        for s in directions[1:]:
            minii = find_minima(X_new,np.array(s))
            X_new = minii[0]
            f_eval += minii[1]
            fx_current.append(func(X_new))
            f_calc.append(f_eval)
        minii = find_minima(X_new,np.array(directions[0]))
        X_new = minii[0]
        f_eval += minii[1]
        fx_current.append(func(X_new))
        f_calc.append(f_eval)
        if abs(new_s(X_new,X1)[1])<=e:
            #print('break_new points',X_new,X1)
            break
        
        if dotproduct(directions,new_s(X_new,X1)[0]):
            return conjugate_direction(X_new)
        X = X_new
        #print(X,func(X))
        directions.pop(0)
        directions.append(new_s(X_new,X1)[0])
        #print(directions)
        iterations+=1
    fx_current.append(func(X_new))
    f_calc.append(f_eval)
    #print(X_new,'minima',penalty(X_new,0)[0],'penality',penalty(X_new,0)[1])
    return X_new


def Bracketing_phase(starting_points,e,c):##this handles solving problem using bracketing phase method.
    global R
    global fx_current
    global f_calc
    for S in starting_points[1:]:
        fx_current = [func(S)]
        f_calc = [0]
        R = 0.01
        X = np.array(S)
        #print(X,penalty(X,R)[0],penalty(X,R)[1])
        max_iter = 100
        X = conjugate_direction(X)
        #print(X,penalty(X,R)[0],penalty(X,R)[1])
        R*=c
        while max_iter:
            prev_X = X
            X = conjugate_direction(X)
            #print(X,penalty(X,R)[0],penalty(X,R)[1])
            penalty_diff = penalty(X,R)[0] - penalty(prev_X,R/c)[0]
            const_violation = penalty(X,R)[1]
            if abs(penalty_diff)<0.0001 and const_violation<0.00001: ##checking termination condition
                print('no constraint voilation')
                break
            R *= c
            max_iter -= 1
        print('Strt_point',S,'\t','minima:',X,'f(x) = ',penalty(X,0)[0],'cons_voilation = ',const_violation) #'\t'''', 'func_evaluations', f_eval,'\t','iterations:',iterations''')
        print()
        print()
##        plt.plot(f_calc[1:],fx_current[1:])
##        plt.xlabel('Func evaluations')
##        plt.ylabel('Penalty value')
##        plt.show()

Bracketing_phase(starting_points,e,c) ##the code starts from this command

	
