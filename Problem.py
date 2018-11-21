import numpy as np
import math

## info is the main function which provides objective function and constraints
## 'p' is the problem number.
def info(p):

    if p==0:
        #general function
        def objective_function(X):
            X = list(X)
            X1 = X[0]
            X2 = X[1]
            result = (X1**2+ X2 -11)**2 + (X1 + X2**2 - 7)**2
            return result

        #inequality constraints(written in from g(x)>=0)  
        def g1(X):
            X = list(X)
            return (X[0]-5)**2 + X[1]**2 -26
        def g2(X):
            X = list(X)
            return X[0]
        def g3(X):
            X = list(X)
            return X[1]
        g = [g1,g2,g3] ## list of equality const functions

        #equality constraints (written in form h(x)=0)
        h = [] ##list of equality const functions
        return objective_function,g,h


    if p==1:
        def objective_function(X): 
            n = len(X)
            X = list(X)
            prod = 1
            for i in X:
                prod*=i
            result = (n**(n/2))*prod
            return (-1)*result
        
        def g1(X): # x>0
            X = list(X)
            negative = 0
            for x in X:
                negative += (min(0,x))**2
            return -1*(negative**(0.5))

        def g2(X): #x<1
            X = list(X)
            more_than_1 = 0
            for x in X:
                more_than_1 += (max(0,x-1))**2
            return -1*(more_than_1 **(0.5))

        def h1(X):
            X = list(X)
            result = [x**2 for x in X]
            return sum(result) - 1
        g = [g1,g2]
        h = [h1]
        return objective_function, g ,h

    if p==2:
        def objective_function(X): 
            X = list(X)
            X1 = X[0]
            X2 = X[1]
            result = (X1-10)**3 + (X2-20)**3
            return result

        def g1(X):
            X = list(X)
            X1 = X[0]
            X2 = X[1]
            val = (X1-5)**2 + (X2-5)**2 - 100
            return val

        def g2(X):
            X = list(X)
            X1 = X[0]
            X2 = X[1]
            val = (X1 - 5)**2 + (X2 - 5)**2 - 82.81
            return (-1)*val

        def g3(X): #13<X1<100
            X = list(X)
            if X[0]<13:
                return (X[0]-13)
            if X[0]>100:
                return -(X[0]-100)
            else:
                return 0

        def g4(X): #0<X2<100
            X = list(X)
            if X[1]<0:
                return X[1]
            if X[1]>100:
                return -(X[1]-100)
            else:
                return 0 
            
        g = [g1,g2,g3,g4]
        h = []        
        return objective_function, g ,h

    if p==3:
        import math
        def objective_function(X): 
            X = list(X)
            X1 = X[0]
            X2 = X[1]
            result = (math.sin(2*math.pi*X1)**3)*math.sin(2*math.pi*X2)/((X1**3)*(X1+X2))
            return (-1)*result

        def g1(X): 
            X = list(X)
            X1 = X[0]
            X2 = X[1]
            val = X1**2 - X2 + 1
            return (-1)*val

        def g2(X): 
            X = list(X)
            X1 = X[0]
            X2 = X[1]
            val = 1 - X1 + (X2-4)**2
            return (-1)*val

        def g3(X): # 0<X1<10
            X = list(X)
            if X[0]<0:
                return X[0]
            if X[0]>10:
                return 10-X[0]
            else:
                return 0

        def g4(X): # 0<X2<10
            X = list(X)
            if X[1]<0:
                return X[1]
            if X[1]>10:
                return 10-X[1]
            else:
                return 0

        g = [g1,g2,g3,g4]
        h = []
        return objective_function, g ,h

    if p==4:
        def objective_function(X):
            X = list(X)
            result = sum(X[:])
            return result

        def g1(X):
            X = list(X)
            X4 = X[3]
            X6 = X[5]
            val = -1 + 0.0025*(X4+X6)
            return (-1)*val

        def g2(X): 
            X = list(X)
            X4 = X[3]
            X5 = X[4]
            X7 = X[6]
            val = -1 + 0.0025*(-X4 + X5 + X7)
            return (-1)*val

        def g3(X): 
            X = list(X)
            X6 = X[5]
            X8 = X[7]
            val = -1 + 0.01*(-X6 + X8)
            return (-1)*val

        def g4(X):
            X = list(X)
            X1 = X[0]
            X4 = X[3]
            X6 = X[5]
            val = 100*X1 - X1*X6 + 833.33252*X4 - 83333.333
            return (-1)*val

        def g5(X):
            X = list(X)
            X2 = X[1]
            X4 = X[3]
            X5 = X[4]
            X7 = X[6]
            val = X2*X4 - X2*X7 - 1250*X4 + 1250*X5
            return (-1)*val

        def g6(X): 
            X = list(X)
            X3 = X[2]
            X5 = X[4]
            X8 = X[7]
            val = X3*X5 - X3*X8 - 2500*X5 + 1250000
            return (-1)*val

        def g7(X): #variables bound
            X = list(X)
            lower = [10,100,100,1,1,1,1,1]
            upper = [10,10,10,1,1,1,1,1]
            error = np.float64(0)
            for i in range(len(X)):
                error += min(0,X[i]-10*lower[i])**2
            for i in range(len(X)):
                error += max(0,X[i]-1000*upper[i])**2
            return error**(0.5)

        g = [g1,g2,g3,g4,g5,g6,g7]
        h = []
        return objective_function, g ,h

    if p==5:
        import math
        def objective_function(X):
            X = list(X)
            prod = np.prod(X)
            return math.exp(prod)

        def g1(X): #variables bound
            X = list(X)
            limit = [2.3,2.3,3.2,3.2,3.2] ## -limit<x<limit
            error = 0
            for i in range(len(X)):
                error += min(0,X[i]+limit[i])**2
            for i in range(len(X)):
                error += max(0,X[i]-limit[i])**2
            return error**(0.5)
        
        def h1(X): 
            X = list(X)
            val = 0
            for x in X:
                val += x**2
            val -= 10
            return val

        def h2(X):
            X = list(X)
            X2 = X[1]
            X3 = X[2]
            X4 = X[3]
            X5 = X[4]
            val = X2*X3 - 5*X4*X5
            return val

        def h3(X):
            X = list(X)
            X1 = X[0]
            X2 = X[1]
            val = X1**3 + X2**3 + 1
            return val
        
        g = [g1]
        h = [h1,h2,h3]
        return objective_function, g, h





    
