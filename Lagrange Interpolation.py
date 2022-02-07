import numpy as np
import os
from sympy.abc import x
from sympy import *
import matplotlib.pyplot as plt
import scipy.interpolate as sci
from sympy import init_printing
init_printing(use_latex=True)
plt.style.use('dark_background')

def grf(symb,f,arx,ary,var):
    x=np.linspace(min(arx),max(arx),10000)
    y=f(x)
    plt.scatter(arx,ary,color="skyblue",marker='.')
    if(var=='y'):
        plt.title("Inverse Lagrange Interpolation\n-Samarth Jain")
    elif(var=="bes"):
        plt.title("Bessel Function: Lagrange Interpolation\n-Samarth Jain")
    elif(var=="lin"):
        plt.title("Linear Interpolation: Intensity vs Photodetector Voltage\n-Samarth Jain")
    else:
        plt.title("Lagrange Interpolation\n-Samarth Jain")
    plt.plot(x,y,color="orange",label="$P({:})={:}$".format(var,symb))
    plt.legend(loc='center',bbox_to_anchor=(0.5,-0.2))
    plt.grid(alpha=0.2)
    plt.show()

def lagf(arrx=None,arry=None,plot=True,var="x"):
    if(arrx==None and arry==None):
        arrx=list(map(float, input("Enter the Array:").strip().split()))
        arry=list(map(float, input("Enter the Array:").strip().split()))
    elif(len(arrx)!=len(arry)):
        return ValueError
    lag=0
    for i in range(len(arrx)):
        f=1
        for j in range(len(arrx)):
            if(i!=j):
                f=f*(x-arrx[j])/(arrx[i]-arrx[j])
        lag+=arry[i]*f
    lgrng=lambdify(x,lag)
    lag=simplify(lag)
    lag=N(lag,4)
    #lag=latex(lag)
    if(var=="y"):
        lag=str(lag)
        lag=lag.replace("x","y")
    if((var=="x") or (var=="y")):
        print("p({:}) = {:}".format(var,lag))
    else:
        print("p({:}) = {:}".format("x",lag))
    if(plot==True):
        grf(lag,lgrng,arrx,arry,var)
    return(lgrng)

def invlag(arrx=None,arry=None,var="y"):
    invlag=lagf(arry,arrx,var=var)
    return(invlag)

def comp(x,y):
    p0=lagf(x,y,plot=False)
    p=sci.lagrange(x,y)
    print("Scipy's Lagrange Polynomial: \n",p)
    xx=np.linspace(min(x),max(x),10000)
    plt.plot(xx,p0(xx),color="green",label="Samarth's function")
    plt.plot(xx,p(xx),color="white",label="Scipy's Interpolation")
    plt.legend()
    plt.show()

def bessel():
    b=[0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0]
    J=[1.0, 0.99, 0.96, 0.91, 0.85, 0.76, 0.67, 0.57, 0.46, 0.34, 0.22, 0.11, 0.0, -0.1, -0.18, -0.26]
    bes=lagf(b,J,var="bes")(2.3)
    besinv=invlag(b,J,var="bes")(0.5)
    print("\nJ(2.3) = {:.5}\n".format(bes))
    print("\u03B2=J\u207b\u00B9(0.5) = {:.5}".format(besinv))

def linear():
    I=[2.81, 3.24, 3.80, 4.30, 4.37, 5.29, 6.03]
    V=[0.5, 1.2, 2.1, 2.9, 3.6, 4.5, 5.7]
    invlin=invlag(I,V,var="lin")(2.4)
    print("I = V\u207b\u00B9(2.4) = {:.5}".format(invlin))

if __name__=='__main__':
    os.system('cls')
    a=[168,120,72,63]
    y=[3,7,9,10]
    print("\nQuestion 2(a): Lagrange Interpolation\n")
    lagf(a,y)
    print("\nQuestion 2(b): Inverse Lagrange Interpolation\n")
    invlag(a,y)
    print("\nQuestion 2(c): Comparison with In-built Lagrange Interpolation\n")
    comp(a,y)
    print("\nQuestion 3(a): Bessel Function\n")
    bessel()
    print("\nQuestion 3(b): Linear Interpolation: Intensity vs Photodetector Voltage\n")
    linear()