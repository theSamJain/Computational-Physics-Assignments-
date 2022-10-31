#finite difference 
import numpy as np 
import matplotlib.pyplot as plt
import pandas as pd
import math
from sklearn.linear_model import LinearRegression as lr

def key(array):
	if ((array[0] == 1) and (array[1] == 0) ):
		return "D"
	elif ((array[0] == 0) and (array[1] == 1) ):
		return "N"
	elif ((array[0] != 0) and (array[1] != 0 )):
		return "R"

def TDMA(l,d,u,B):
    #l = Lower Diagonal vector, d = Main Diagonal vector
    #u = Upper Diagonal vector, B = solution vector
	n = len(B)
	w = np.zeros(n-1)
	g = np.zeros(n)
	x = np.zeros(n)
	w[0] = u[0]/d[0]
	g[0] = B[0]/d[0]
	for i in range(1,n-1):
		w[i] = u[i]/(d[i] - l[i-1]*w[i-1])
	for i in range(1,n):
		g[i] = (B[i] - l[i-1]*g[i-1])/(d[i] - l[i-1]*w[i-1])
	x[n-1] = g[n-1]
	for i in range(n-1,0,-1):
		x[i-1] = g[i-1] - w[i-1]*x[i]
	return x	
		
def finiteDiff(a, b, n, px, qx, rx, ic, bc): #define px , qx , rx using lambda x : function(x) method 
	h = (b-a)/(n-1)
	
	p = [px(a+i*h) for i in range(n)]
	q = [qx(a+i*h) for i in range(n)] 
	r = [rx(a+i*h) for i in range(n)] 	

	xea = key(ic)
	xeb = key(bc)
	
	d = [2 + h**2*q[i] for i  in range(n)]
	u = [-1+(h/2)*p[i] for i in range(n)]
	l = [-1-(h/2)*p[i] for i in range(n)] 
	B =  [-h**2*r[i] for i in range(1,n-1)]
	
	if  xea == "D":
		a11= 1
		a12 =0 
		b1 = ic[2] # ic[0] -> ic[2]
	elif xea == "N":
		a11 = d[0]
		a12 = -2
		b1 = -h**2*r[0]+2*h*l[0]*ic[2]
	elif xea == "R":
		a11 = d[0]+2*h*l[0]*(ic[0]/ic[1])
		a12 =-2
		b1 = (-h**2)*r[0]+2*h*l[0]*(ic[2]/ic[1]) # alpha3 = ic[2]
	if  xeb == "D":
		anp1np1= 1
		anp1n =0 
		bnp1 = bc[2] # bc[0]->bc[2]
	elif xeb == "N":
		anp1np1 = d[-1]
		anp1n = -2
		bnp1 = -h**2*r[-1]-2*h*u[n-1]*bc[2]
	elif xeb == "R":
		anp1np1 = d[-1]-2*h*u[n-1]*(bc[0]/bc[1])
		anp1n =-2
		bnp1 = (-h**2)*r[-1]-2*h*u[n-1]*(bc[2]/bc[1])

	B.insert(0,b1)
	B.append(bnp1)
	B=  np.array(B).reshape(-1,1)

	d[0] = a11
	d[n-1] =  anp1np1
	u[0] =  a12
	l[n-1] = anp1n
	del l[0]
	u.pop()
	w = TDMA(l,d,u,B)
	return(w)

"""
	matrix = np.zeros((int(n),int(n)))
	for j in range(n):
		for k in range(n):
			if j== 0 :
				if k == 0 :
					matrix[j][k] = a11
				elif k==1:
					 matrix[j][k] = a12
			elif j ==n :
				if k == n-1 : 
					matrix[j][k] = anp1n
				elif k == n :
					matrix[j][k] = anp1np1				
			elif j==k and j != 0 and j != n :
				matrix[j][k] =  d[j]
			elif j+1==k and j != 0 and j != n  : 
		 		matrix[j][k] = u[j]
			elif j == k +1 and j != 0 and j != n  :
				matrix[j][k]  = l[k]
	w = np.linalg.solve(matrix,B)
	print(matrix)
	return(w)
"""

def finiteDiffExt(a, b, px, qx, rx, ic, bc, k ,exactfunc ): # exact function should be in the form of lambda function 
    N = [2+2**i for i in range(k)]
    H = [(b-a)/(i-1) for i in N]
    xx = list(map(lambda x,y : [ a + i*x for i in range(y)], H , N ))
    Exact_n = list(map(lambda x : exactfunc(np.array(x)) , xx))  
    W_n = list(map(lambda z : finiteDiff(a,b,z,px , qx ,rx , ic ,bc) , N))
    wi_ui = list(map(lambda x,y : [np.abs(element1 - element2) for (element1, element2) in zip(x, y)] , W_n , Exact_n))
    max_abs_err = [max(i) for i in wi_ui]
    sqr = list(map(lambda x: np.array(x)**2 , wi_ui))
    rms = list(map(lambda x,y : np.sqrt( (1/(1+x))*np.sum(y)) , N , sqr ) )
    return xx, W_n, N, rms, max_abs_err

if __name__ == "__main__":
    #b(1)
    a , b = 0 , 1
    n = 5
    h = (b-a)/(n-1)
    px  = lambda x : 0
    qx  = lambda x : np.pi**2
    rx  = lambda x : -2*np.pi**2 * np.sin(np.pi*x)
    ic   = [1,0,0] # input from user 
    bc  = [1,0,0]

    xnum1 = [a+i*h for i in range(n)]
    ynum1 = finiteDiff(a,b,n,px,qx,rx,ic,bc)
    yexpr1 =  lambda x: np.sin(np.pi*x)
    xaxis1 = np.linspace(a,b,100)
    yexact1 = yexpr1(xaxis1)

    xx1, w_n1, N1, erreq1, maxerr1 = finiteDiffExt(a, b, px, qx, rx, ic, bc, 6, yexpr1)

    plt.plot(xaxis1, yexact1, label  = 'Analytic' )
    plt.scatter(xnum1, ynum1 , marker = '.' , color = 'red' , label = 'N = {:}'.format(n))
    plt.xlabel("x" )
    plt.ylabel("f(x)")
    plt.title("Finite Difference Method - Equation 1\nSamarth Jain and Swarnim Gupta")
    plt.legend()
    plt.show()

    #b(2)
    a, b = 0 , np.pi/2
    n = 5
    h = (b-a)/(n-1)
    px  = lambda x : 0
    qx  = lambda x : -1
    rx  = lambda x : np.sin(3*x)
    ic = [1,1,-1] # input from user 
    bc = [0,1,1]

    xnum2 = [a+i*h for i in range(n)]
    ynum2 = finiteDiff(a,b,n,px,qx,rx,ic,bc)
    yexpr2 = lambda x: (3/8)*(np.sin(x))-np.cos(x)-(1/8)*np.sin(3*x)
    xaxis2 = np.linspace(0,np.pi/2,10)
    yexact2 = yexpr2(xaxis2) 

    xx2, w_n2, N2, erreq2, maxerr2 = finiteDiffExt(a, b, px, qx, rx, ic, bc, 6, yexpr2)

    plt.plot(xaxis2, yexact2, label  = 'Analytic' )
    plt.scatter(xnum2, ynum2, marker = '.' , color = 'red' , label = 'N = {:}'.format(n))
    plt.xlabel("x" )
    plt.ylabel("f(x)")
    plt.title("Finite Difference Method - Equation 2\nSamarth Jain and Swarnim Gupta")
    plt.legend()
    plt.show()

    # c part
    ye1 = [yexpr1(i) for i in xnum1]
    print("\nN = %d\nEquation 1\n"%(n-2))
    dict1 = {"x": xnum1,
            "y calc": ynum1,
            "y exact": ye1,
            "Absolute Error": np.abs(ye1 - ynum1)
            }
    dt1 = pd.DataFrame(dict1)
    print(dt1)

    ye2 = [yexpr2(i) for i in xnum2]
    print("\n\nEquation 2")
    dict2 = {"x": xnum2,
            "y calc": ynum2,
            "y exact": ye2, 
            "Absolute Error": np.abs(ye2 - ynum2)
            }
    dt2 = pd.DataFrame(dict2)
    print(dt2)

    print("\n")

    #e part
	#function 1
    for m in range(len(xx1)):
      for n in range(len(w_n1)):
        if m==n:
            plt.plot(xx1[m], w_n1[n], marker = ".", label = " N  = {:}".format((2)**(m+1)))
    plt.plot(xaxis1, yexact1, label = "Analytic")
    plt.legend()
    plt.title("Performance of Finite Difference Method w.r.t. N \nSamarth Jain and Swarnim Gupta")
    plt.xlabel("x")
    plt.ylabel("f(x)")
    plt.show()

    #function 2
    for o in range(len(xx2)):
      for p in range(len(w_n2)):
        if o==p :
            plt.plot(xx2[o], w_n2[p], marker = ".", label = "N = {:}".format((2)**(o+1)))
    plt.plot(xaxis2, yexact2, label = "Analytic")
    plt.legend()    
    plt.title("Performance of Finite Difference Method w.r.t. N \nSamarth Jain and Swarnim Gupta")
    plt.xlabel("x")
    plt.ylabel("f(x)")
    plt.show()

    # table for rms and mas abs error
    dict_err = {"N": N1,
                "RMS Error - Eq1": erreq1, 
                "Max. Abs. Error - Eq1": maxerr1, 
                "RMS Error - Eq2": erreq2, 
                "Max. Abs. Error - Eq2": maxerr2
                }
    dt_err = pd.DataFrame(dict_err)
    print("\nError Analysis\n")
    print(dt_err)

    # f part
    log_N = [math.log(i) for i in N1]
    log_erreq1 = [math.log(i) for i in erreq1]
    log_erreq2 = [math.log(i) for i in erreq2]
    naxis = np.linspace(log_N[0], log_N[-1], 100)

    fit1 = lr().fit(np.array(log_N).reshape(-1, 1), np.array(log_erreq1).reshape(-1, 1))
    slope1 = fit1.coef_
    intc1 = fit1.intercept_
    line1 = naxis*slope1+intc1
    plt.plot(naxis, line1.flatten(), label = "Eq1: Fitted Line", zorder = -1)

    fit2 = lr().fit(np.array(log_N).reshape(-1, 1), np.array(log_erreq2).reshape(-1, 1))
    slope2 = fit2.coef_
    intc2 = fit2.intercept_
    line2 = naxis*slope2+intc2
    plt.plot(naxis, line2.flatten(), label = "Eq2: Fitted Line", zorder = -1)

    print("Slope for Equation 1: ", slope1.flatten())
    print("Slope for Equation 2: ", slope2.flatten())

    plt.scatter(log_N, log_erreq1, label = "RMS Error for Equation 1")
    plt.scatter(log_N, log_erreq2, label = "RMS Error for Equation 2")
    plt.legend()
    plt.xlabel("Log N")
    plt.ylabel("Log RMS Error")
    plt.title("Error Analysis\nSamarth Jain and Swarnim Gupta")
    plt.show()