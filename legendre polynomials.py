import numpy as np
import scipy.integrate as intg
import math as mt
import os
import sympy
from sympy.abc import x
import matplotlib.pyplot as plt
from scipy.special import legendre as lgdr
import random
plt.style.use('dark_background')

def krondelta(n,m):
    if(n==m):
        return(1)
    else:
        return(0) 

def grf(f,fdv,n):
    x=np.linspace(-1,1,10**4)
    y=[]; ydv=[]
    for i in range(len(x)):
        y.append(f(x[i]))
    plt.style.use('dark_background')
    plt.plot(x,y,color="orange",label="$P_{:}$".format(n))
    if(fdv!=None):
        for i in range(len(x)):
            ydv.append(fdv(x[i]))
        plt.plot(x,ydv,color="pink",label="$P_{:}'$".format(n))
    plt.grid(alpha=0.2)
    #plt.ylim(-1,1)         #When we need to limit y
    plt.legend(loc="upper left")
    plt.title("Legendre Polynomial\n-Samarth Jain")
    plt.show()

def table(data,head):
    line='_'*len(head)*12+'____'
    for i in head:
        print("{0:^12}".format(i),end=" ")
    print("\n",line)
    for row in data:
        for val in row:
            print("{0:^12.2}".format(val), end=" ")
        print("\n") 

def randgen():
    arr=[]
    for i in range (1,4):
        lf=lp(i)[1]       #Q2a : call correct function according to the question
        #lf=ld(i)[1]        #Q2b : call correct function according to the question
        arr.append(lf)
    x=[]; y1=[]; y2=[]; y3=[]
    for i in range(1001):
        #x.append(random.uniform(-1,1))
        #z=-1+0.002*i       #Only for leg00_
        #x.append(z)
        y1.append(arr[0](x[i]))
        y2.append(arr[1](x[i]))
        y3.append(arr[2](x[i]))
    dout=np.column_stack((x,y1,y2,y3))
    np.savetxt(r'C:\Users\jains\Desktop\College Physics\Python MPhy\leg00.dat',dout,fmt=("%.5f")) #leg00 leg01

def gamma():
    np1=float(input("Enter the value of n for \u0393(n): "))
    if (np1<=0):
        print("Gamma function is divergent for negative real numbers :)")
        exit()
    f=lambda x: (x**(np1-1))*np.exp(-x)
    gamma_np1,err=intg.quad(f,0,np.inf)
    print("\u0393({:}) = {:.5}".format(float(np1), float(gamma_np1)))

def lp(nval=None,out=None,plot=False):
    if(nval==None):
        n=int(input("Enter the value of n: "))
    else:
        n=nval
    pn=0                    #Legendre Polynomial
    if(n%2==0):
        m=int(n/2)
    elif(n%2!=0):
        m=int((n-1)/2)
    for i in range (0,m+1):
        pn+=(((-1)**i)*mt.factorial((2*n-2*i))*(x**(n-2*i)))/((2**n)*mt.factorial(i)*mt.factorial(n-i)*mt.factorial(n-2*i))
    lf=sympy.lambdify(x,pn)
    if(out==None):          #When we want to print Pn
        print("\nP{:} = {:}".format(int(n),pn))
    if(plot==True):         #When we want to plot only Pn
        grf(lf,None,n)         
    return(pn,lf,n)         

def ld(out=None,plot=False):
    f,lf,n=lp()
    sfdv=f.diff(x)
    lfdv=sympy.lambdify(x,sfdv)
    if(f.is_constant is True):  #Because Matplotlib has hard time plotting constant functions!
        lf=np.ones(101)
    if(sfdv.is_constant is True):
        lfdv=np.ones(101)
    if(out==None):
        print("\nP{:}' = {:}".format(int(n),sfdv))
    if(plot==True):
        grf(lf,lfdv,n)          #When we want to plot both Pn and Pn'
    return(sfdv,lfdv)

def lcomp():
    f=lp()
    pn=lgdr(f[2])
    print("\nP{:} (Scipy's) = \n{:}".format(int(f[2]),pn))
    grf(pn,f[1],f[2])

def Q2a():
    #randgen()
    x,y1,y2,y3=np.loadtxt(r'C:\Users\jains\Desktop\College Physics\Python MPhy\leg00.dat',unpack=True)
    data=np.loadtxt(r'C:\Users\jains\Desktop\College Physics\Python MPhy\leg00.dat')
    head=["x","P1(x)","P2(x)","P3(x)"]
    table(data,head)
    plt.scatter(x,y1,label="$P_1$",color='white',marker=".",s=1.5)
    plt.scatter(x,y2,label="$P_2$",color='pink',marker=".",s=1.5)
    plt.scatter(x,y3,label="$P_3$",color='skyblue',marker=".",s=1.5)
    plt.title("$P_1(x)$, $P_2(x)$, $P_3(x)$ - 100 Random Points\n-Samarth Jain")
    plt.legend()
    plt.show()

def Q2b():
    #randgen()
    x,y1,y2,y3=np.loadtxt(r'C:\Users\jains\Desktop\College Physics\Python MPhy\leg01.dat',unpack=True)
    yp1=[]
    p1=lp(1)[1]
    for i in range(len(x)):
        yp1.append(p1(x[i]))
    data=np.loadtxt(r'C:\Users\jains\Desktop\College Physics\Python MPhy\leg01.dat')
    head=["x","P1'(x)","P2'(x)","P3'(x)"]
    table(data,head)
    plt.scatter(x,y1,label="$P_1'$",color='white',marker=".",s=1.5)
    plt.scatter(x,y2,label="$P_2'$",color='pink',marker=".",s=1.5)
    plt.scatter(x,y3,label="$P_3'$",color='skyblue',marker=".",s=1.5)
    plt.scatter(x,yp1,label="$P_1$",color='yellow',marker=".",s=1.5)
    plt.title("$P_1'(x)$, $P_2'(x)$, $P_3'(x), P_1(x)$ - 100 Random Points\n-Samarth Jain")
    plt.legend()
    plt.show()

def Q2c1():
    n=2
    p2=lp(n)[1]
    x,pdv1,pdv2=np.loadtxt(r'C:\Users\jains\Desktop\College Physics\Python MPhy\leg01.dat',unpack=True, usecols=[0,1,2])
    rhs=[]; lhs=[]; z=[]
    for i in range(len(x)):
        z1=p2(x[i])
        z.append(z1)
        lhs.append(2*z1)
    for i in range(len(x)):
        z2=x[i]*pdv2[i]-pdv1[i]
        rhs.append(z2)
    n_0=[n]*len(x)
    n_1=[n-1]*len(x)
    data=np.column_stack((n_0,n_1,x,z,pdv1,pdv2,lhs,rhs))
    np.savetxt(r'C:\Users\jains\Desktop\College Physics\Python MPhy\leg02.dat',data,header="n_0,n_1,x,z,pdv1,pdv2,lhs,rhs",fmt='%.5f')
    head=["n","n-1","x","P2(x)","P1'(x)","P2'(x)","LHS","RHS"]
    table(data,head)
    print("\nYay! n*Pn(x)=x*Pn-1'(x)-Pn-1'(x) is True! (Atleast for n=2)\n")

def Q2c2():
    n=2
    x,p1,p2,p3=np.loadtxt(r'C:\Users\jains\Desktop\College Physics\Python MPhy\leg00.dat',unpack=True)
    rhs=[]; lhs=[]; z=[]
    for i in range(len(x)):
        z1=5*x[i]*p2[i]
        lhs.append(z1)
    for i in range(len(x)):
        z2=3*p3[i]+2*p1[i]
        rhs.append(z2)
    n_0=[n]*len(x)
    n_1=[n-1]*len(x)
    np1=[n+1]*len(x)
    data=np.column_stack((n_0,n_1,np1,x,p1,p2,p3,lhs,rhs))
    np.savetxt(r'C:\Users\jains\Desktop\College Physics\Python MPhy\leg03.dat',data,header='n_0,n_1,np1,x,p1,p2,p3,lhs,rhs',fmt="%.5f")
    head=["n","n-1","n+1","x","P1(x)","P2(x)","P3(x)","LHS","RHS"]
    table(data,head)
    print("\nYay! (2*n+1)*x*Pn(x)=(n+1)*Pn+1(x)-n*Pn-1'(x) is True! (Atleast for n=2)\n")

def Q2c3():
    n=3
    x,p1,p2,p3=np.loadtxt(r'C:\Users\jains\Desktop\College Physics\Python MPhy\leg00.dat',unpack=True)
    rhs=[]; lhs=[]; z=[]
    for i in range(len(x)):
        z1=3*p3[i]
        lhs.append(z1)
    for i in range(len(x)):
        z2=5*x[i]*p2[i]-2*p1[i]
        rhs.append(z2)
    n_0=[n]*len(x)
    n_1=[n-1]*len(x)
    n_2=[n-2]*len(x)
    data=np.column_stack((n_0,n_1,n_2,x,p1,p2,p3,lhs,rhs))
    np.savetxt(r'C:\Users\jains\Desktop\College Physics\Python MPhy\leg04.dat',data,header='n_0,n_1,n_2,x,p1,p2,p3,lhs,rhs')
    head=["n","n-1","n-2","x","P1(x)","P2(x)","P3(x)","LHS","RHS"]
    table(data,head)
    print("\nYay! n*Pn(x)=(2*n-1)*x*Pn-1(x)-(n-1)*Pn-2(x) is True! (Atleast for n=3)\n")

def Q2d_i():            #File version of Q2d with 1000 points
    x,p1,p2,p3=np.loadtxt(r'C:\Users\jains\Desktop\College Physics\Python MPhy\leg00_.dat',unpack=True,usecols=[0,1,2,3])
    data=[p1,p2,p3]
    rhsCol=[]; intglCol=np.array([]); n=[]; m=[]
    for i in range(3):
        pn=data[i]
        for j in range(3):
            y=[]
            pm=data[j]
            for k in range(len(x)):
                y.append(pn[k]*pm[k])
            intgl=intg.simpson(y,x)
            intglCol=np.append(intglCol,intgl)
            rhs=2/(2*(i+1)+1)*krondelta(i+1,j+1)
            rhsCol.append(rhs)
        n.append(i+1)
        m=n
    print("For Integral:\n{:^12} {:^10} {:^10} {:^10}\n".format("N/M",n[0],n[1],n[2]))
    for i in range(3):
        print("{:^12} {:^10.4} {:^10.4} {:^10.4}".format(m[i], intglCol[i],intglCol[i+1],intglCol[i+2]))
        print("\n")
    print("\nFor Kronecker Delta:\n{:^12} {:^10} {:^10} {:^10}\n".format("N/M",n[0],n[1],n[2]))
    for i in range(3):
        print("{:^12} {:^10.4} {:^10.4} {:^10.4}".format(m[i], rhsCol[i],rhsCol[i+1],rhsCol[i+2]))
        print("\n")
    
def Q2d_ii():           #Another version of Q2d
    rhsCol=[]; intglCol=[]; n=[]; m=[]
    for i in range(1,4):
        pn=lp(int(i),out=0)[0]
        for j in range(1,4):
            n.append(i)
            m.append(j)
            pm=lp(int(j),out=0)[0]
            intgrnd=pn*pm
            intgrnd=sympy.lambdify(x,intgrnd)
            intgl,err=intg.quad(intgrnd,-1,1)
            intglCol.append(intgl)
            rhs=2/(2*i+1)*krondelta(i,j)
            rhsCol.append(rhs)
    data=np.column_stack((n,m,intglCol,rhsCol))
    head=["n","m","Integral","RHS"]
    table(data,head)

if __name__=='__main__':
    os.system('cls')
    print("\nQuestion 1(a): Gamma Function\n")
    gamma()
    print("\n\nQuestion 1(b): Legendre Polynomials for n\n")
    lp(plot=True)
    print("\n\nQuestion 1(c): Derivatives of Legendre Polynomials for n\n")
    ld(plot=True)
    print("\n\nQuestion 1(d): Comparison of Self-made and In-built Functions\n")
    lcomp()
    print("\n\nQuestion 2(a): Storing Legendre Polynomials for random x in a file\n")
    Q2a()
    print("\n\nQuestion 2(b): Storing Derivatives of Legendre Polynomials for random x in a file\n")
    Q2b()
    print("\n\nQuestion 2(c)(1): Proving: n*Pn(x)=x*Pn-1'(x)-Pn-1'(x)")
    Q2c1()
    print("\n\nQuestion 2(c)(2): Proving: (2*n+1)*x*Pn(x)=(n+1)*Pn+1(x)-n*Pn-1'(x)")
    Q2c2()
    print("\n\nQuestion 2(c)(3): Proving: (2*n-1)*x*Pn-1(x)-(n-1)*Pn-2(x)")
    Q2c3()
    print("\n\nQuestion 2(d): Proving Orthonormality Relation using values from 'leg00.dat'\n")
    Q2d_i()
    print("\n\nAdditional Version of Question 2(d): Proving Orthonormality Relation for any x\n")
    Q2d_ii()