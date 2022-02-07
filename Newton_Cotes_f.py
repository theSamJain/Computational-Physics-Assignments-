import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mt
import scipy.integrate as integrate

# note: Q2c has been done in two ways. Check the function for description.

def trapz(f,a,b,n):
    h=(b-a)/n
    y=[]
    for i in range(n+1):    
        y.append(f(a+i*h))              #Value of f(x) at the nodal points
    trap=h*(f(a)+f(b))/2
    for j in range(1,len(y)-1):
        trap=trap+h*(y[j])
    return(trap)

def simps(f,a,b,n):
    h2=(b-a)/(2*n)
    simp=h2*(f(a)+f(b))/3
    for i in range(1,2*n): 
        if(i%2==0):
            simp=simp+2*h2*f(a+i*h2)/3  #y at Even Nodes 
        elif(i%2==1):
            simp=simp+4*h2*f(a+i*h2)/3  #y at Odd Nodes
    return(simp)

def grf(f,a,b):
    N=np.arange(1,100,1)    #Number of Subintervals
    H=(b-a)/N
    It=[]
    for i in N:             #Trapezoidal
        z=trapz(f,a,b,i)
        It.append(z)

    H2=(b-a)/(2*N)
    Is=[]
    for j in 2*N+1:         #Simpson
        z=simps(f,a,b,j)
        Is.append(z)

    #Plotting

    ax=plt.subplot()
    ax.get_yaxis().set_major_formatter(mt.ScalarFormatter())
    ax.get_xaxis().set_major_formatter(mt.ScalarFormatter())
    ax.get_yaxis().set_minor_formatter(mt.ScalarFormatter())
    ax.get_xaxis().set_minor_formatter(mt.ScalarFormatter())
    plt.yscale("log")
    plt.xscale("log")
    plt.grid(b=True,which='major', axis='both')
    plt.grid(b=True,which='minor', axis='both')

    plt.plot(H,It,label="Composite Trapezoidal Method ($1\leq N \leq 100$)",marker=".")
    plt.plot(H2,Is,label="Composite Simpson Method ($1\leq N \leq 100$)",marker=".")

    plt.xlabel("Step Size (h)")
    plt.ylabel("I(h) - Integral varying with h")
    plt.title("Examination of Convergence of Trapezoidal and Simpson Method \n-Samarth Jain")
    plt.legend()
    plt.show()

def Q2a():
    #Trapezoidal
    f=eval("lambda x:"+input("Enter the Function : "))
    a=float(input("a = "))
    b=float(input("b = "))
    n=int(input("N = "))
    trap=trapz(f,a,b,n)
    print("Integral is {:.8} using Trapezoidal Method".format(float(trap)))
    
def Q2b():
    #Simpson
    f=eval("lambda x:"+input("Enter the Function : "))
    a=float(input("a = "))
    b=float(input("b = "))
    n=int(input("N = "))
    simp=simps(f,a,b,n)
    print("Integral is {:.8} using Simpson Method".format(float(simp)))

#This way takes an analytical function and assumes its integral and then proceeds
def Q2c_i():
    #Example function = x**2
    #Analytical Integral = x**3/3
    f=eval("lambda x:"+"x**2")
    a=float(input("a = "))
    b=float(input("b = "))
    n=int(input("N = "))
    anlytc=(b**3-a**3)/3
    trap=trapz(f,a,b,n)
    err=abs(anlytc-trap)
    simp=simps(f,a,b,n)
    err2=abs(anlytc-simp)
    print("For f(x) = x**2 : \nTruncation Error(Trapezoidal Method) = {:e}\nTruncation Error(Simpson Method) = {:e}".format(float(err),float(err2)))

#This method takes arbitrary func from user and uses scipy.quad to calc integral
#So it is numerical vs numerical and uses an in-built library!!!
def Q2c_ii():
    f=eval("lambda x:"+input("Enter the Function : "))
    a=float(input("a = "))
    b=float(input("b = "))
    n=int(input("N = "))
    qd,errq=integrate.quad(f,a,b)
    trap=trapz(f,a,b,n)
    err=abs(qd-trap)
    simp=simps(f,a,b,n)
    err2=abs(qd-simp)
    print("\nTruncation Error in Quadrature = {:e}".format(errq))
    print("Truncation Error between Quadrature and Trapezoidal Method = {:e}\nTruncation Error between Quadrature and Simpson Method = {:e}".format(float(err),float(err2)))

def Q2d():
    f=eval("lambda x:"+input("Enter the Function : "))
    a=float(input("a = "))
    b=float(input("b = "))
    grf(f,a,b)

def Q3a():
    v=np.array([0.0, 0.5, 2.0, 4.05, 8.0, 12.5, 18.0, 24.5, 32.0, 40.5, 50.0])
    c=np.array([0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
    h=(c[-1]-c[0])/(len(c)-1)
    pwrt=h*(v[0]+v[-1])/2 
    for j in range(1,len(v)-1):
        pwrt=pwrt+h*(v[j])          #Power using Trapezoidal Method

    pwrs=h*(v[0]+v[-1])/3
    for i in range(1,len(v)-1): 
        if(i%2==0):
            pwrs=pwrs+2*h*v[i]/3
        elif(i%2==1):
            pwrs=pwrs+4*h*v[i]/3    #Power using Simpson Method

    print("\nPower is {:.8} Joules (Using Trapezoidal Method)\nPower is {:.8} Joules (Using Simpson Method)".format(float(pwrt), float(pwrs)))

    plt.scatter(c,v,marker=".")
    plt.ylabel("Volage (V)")
    plt.xlabel("Current (mA)")
    plt.show()

def Q3b():
    f=eval("lambda x:"+input("Enter the Function : "))
    a=float(input("a = "))
    b=float(input("b = "))
    n=int(input("N = "))
    h=(b-a)/n
    #Trapezoidal
    trap=trapz(f,a,b,n)
    print("Integral is {:.8} using Trapezoidal Method".format(float(trap)))
    #Simpson
    simp=simps(f,a,b,n)
    print("Integral is {:.8} using Simpson Method".format(float(simp)))
    
    print("f(h) = {:.8}".format(float(f(h))))
    
    print("\nPLEASE NOTE: For graph, the step size 'h' is being recalculated for a given range of 'n'\n")
    #Plotting
    grf(f,a,b)

if __name__ == '__main__':
    print("\nQuestion 2(a)")
    Q2a()
    print("\n\nQuestion 2(b)")
    Q2b()
    print("\n\nQuestion 2(c) - Analytical-Numerical (Example Function = x**2) ")
    Q2c_i()
    print("\n\nQuestion 2(c) - Numerical-Numerical")
    Q2c_ii()
    print("\n\nQuestion 2(d)")
    Q2d()
    print("\n\nQuestion 3(a)")
    Q3a()
    print("\n\nQuestion 3(b)")
    Q3b()