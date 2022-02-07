import numpy as np
import sympy as sp
from sympy.abc import x,y
import matplotlib.pyplot as plt
import os

def table(head,data,splfr=None,title=None):
    t=0
    print("{:^}".format(title))
    line='_'*len(head)*15+'_'*10
    for i in head:
        print("{:^13}".format(i),end=" ")
    print("\n",line)
    for row in data:
        for val in row:
            print("{:^13.5e}".format(val), end=" ")
        if((splfr!=None) and (t==0)):
            for i in range(len(splfr)):
                print("{:^13}".format(splfr[i]),end=" ")
                t=1
        print("\n") 

def grf(xarr,yarr,title,xt=None,yt=None,yarr2=None,labels=None,onept=None,app=None):
    if(labels==None):
        labels=[]
        labels.append("Integrated Points")
    plt.scatter(xarr,yarr,marker=".",label=labels[0])
    if(yarr2!=None):
        plt.scatter(xarr,yarr2,marker=".",label=labels[1])
    plt.title(title)
    if(onept!=None):
        plt.scatter(onept[0],onept[1],marker="o",label="Average Life")
    if(xt!=None):
        plt.xlabel(xt)
        plt.ylabel(yt)
    else:
        plt.xlabel("X Axis")
        plt.ylabel("Y Axis")
    if(app!=None):
        ticks=np.arange(0,app+onept[0]/2,onept[0])
        lbls=[i for i in range(len(ticks))]
        plt.xticks(ticks,lbls)
    plt.legend()
    plt.grid(alpha=0.5)
    plt.show()

def analyt(yo,x,tau):
    y=[]
    for i in range(len(x)):
        y.append(yo*np.exp(-x[i]/tau))
    return(y)

def elr(yo,a,b,n=None,h=None,f=None,prnt=None,plot=None):
    if(n!=None):
        h=(b-a)/n       #step size
    elif(h!=None):
        n=int((b-a)/h)
    elif((n!=None) and (h!=None)):
        return os.error("Can not give both n and h")
    xo=a
    if (f==None):
        f=eval("lambda y,x:"+input("Enter the slope function y': "))
        slopefunc=f
    else:
        slopefunc=sp.lambdify([y,x],f)
    yarr=[]; xarr=[]
    yin=yo; xin=xo
    for i in range(n+1):
        yarr.append(yin)
        xarr.append(xin)
        ynext=yin+h*slopefunc(yin,xin)
        yin=ynext
        xin=xo+h*(i+1)
    if(prnt==None):
        head=["x","y(x)","a","b","N","h"]
        title="Euler's Method (Forward)\n"
        data=np.column_stack((xarr,yarr))
        fr=[a,b,n,h]
        table(head,data,fr,title)
    if(plot==None):
        f=sp.sympify(f)
        # print(f)
        title="Euler's Method (Forward) - Samarth Jain\ny'=$"+sp.latex(f)+"$"
        grf(xarr,yarr,title)
    return(xarr,yarr)

def rk2(yo,a,b,n=None,h=None,f=None,prnt=None,plot=None):
    if(n!=None):
        h=(b-a)/n       #When we have n
    elif(h!=None):      #When we have h
        n=int((b-a)/h)
    elif((n!=None) and (h!=None)):
        return os.error("Can not give both n and h")
    xo=a
    if (f==None):
        f=eval("lambda y,x:"+input("Enter the slope function y': "))
        slopefunc=f
    else:
        slopefunc=sp.lambdify([y,x],f)
    yarr=[]; xarr=[]
    yin=yo; xin=xo
    for i in range(n+1):
        yarr.append(yin)
        xarr.append(xin)
        k1=h*slopefunc(yin,xin)
        # ynext=yin+k1; xnext=xin+h         #Answer not matching!
        ynext=yin+k1/2; xnext=xin+h/2
        k2=h*slopefunc(ynext,xnext)
        yf=yin+k2
        yin=yf
        # xin=xnext                         #For the other method
        xin=xnext+h/2
    if(prnt==None):
        head=["x","y(x)","a","b","N","h"]
        data=np.column_stack((xarr,yarr))
        fr=[a,b,n,h]
        title="Runge-Kutta Method\n"
        table(head,data,fr,title)
    if(plot==None):
        f=sp.sympify(f)
        # print(f)
        title="Runge-Kutta Method - Samarth Jain\ny'=$"+sp.latex(f)+"$"
        grf(xarr,yarr,title)
    return(xarr,yarr)

def rk_elr_comp(f,a,b,yo,n=None,h=None):
    rk=rk2(a=a,b=b,n=n,yo=yo,f=f,prnt=False,plot=False)
    eul=elr(a=a,b=b,n=n,yo=yo,f=f,prnt=False,plot=False)
    plt.scatter(rk[0],rk[1],label="From Runge-Kutta",marker=".")
    plt.scatter(eul[0],eul[1],label="From Euler",marker=".")
    plt.legend()
    plt.grid(alpha=0.5)
    err=[abs((i-j)) for (i,j) in zip(eul[1],rk[1])]
    data=np.column_stack((eul[0],eul[1],rk[1],err))
    head=["x","y(x)-Euler","y(x)-RK","Abs. Error","h"]
    tlt="Euler v. Runge-Kutta\n-Samarth Jain"
    h=[(b-a)/n]
    table(head=head,data=data,title=tlt,splfr=h)
    plt.xlabel("X axis")
    plt.ylabel("Y axis")
    plt.title(tlt)
    plt.show()

def hvsi(a,b,f):
    fig,axs=plt.subplots(1,2)
    c=5
    n=[10**i for i in range(1,c+1)]
    h=[(b-a)/n[i] for i in range(c)]
    for i in range(c):
        xx,yy=elr(a=a,b=b,n=n[i],yo=1,f=f,prnt=False,plot=False)
        axs[0].plot(xx,yy,label="h=%.1e"%h[i])
    for i in range(c):
        xx,yy=rk2(a=a,b=b,n=n[i],yo=1,f=f,prnt=False,plot=False)
        axs[1].plot(xx,yy,label="h=%.1e"%h[i])
    axs[0].set(title="Euler's Method",xlabel='x', ylabel='f(x)')
    axs[1].set(title="Runge-Kutta's Method",xlabel='x', ylabel='f(x)')
    axs[0].grid()
    axs[1].grid()
    axs[0].legend()
    axs[1].legend()
    fig.suptitle("Variation with Step Size(h)")
    plt.tight_layout()
    plt.show()

def decay(b=None,n=None):
    t1="Radio Active Decay - Samarth Jain"
    t2=t1+"\n$\dot N=-\lambda N=\dfrac{-N}{\\tau}$"
    No=20000
    hl=4
    tau=hl/np.log(2)
    if((b==None) and (n==None)):
        b=5*tau
        n=20
    f=-y/tau
    h=tau/n
    dec=rk2(a=0,b=b,h=h,yo=No,f=f,plot=False,prnt=False)        #Passing h instead of n BUT can pass n instead of h as well!
    dece=elr(a=0,b=b,h=h,yo=No,f=f,plot=False,prnt=False)[1]
    al=analyt(yo=No,x=dec[0],tau=tau)
    errrk=[abs(i-j) for (i,j) in zip(al,dec[1])]
    erre=[abs(i-j) for (i,j) in zip(al,dece)]
    pt=analyt(yo=No,x=[tau],tau=tau)
    onept=[tau,pt]
    data=np.column_stack((dec[0],dec[1],dece,errrk,erre))
    head=["Time (years)","N-RK","N-Euler","RK Abs. Err","Elr Abs. Err","No","Half Life","Step Size"]
    fr=[No,hl,h]
    table(head,data,fr,t1)
    labels=["Runge-Kutta","Euler"]
    grf(dec[0],dec[1],t2,"Tau (%.4f years)"%tau,"Population (N)",dece,labels,onept,app=b)
    grf(dec[0],errrk,title="Radioactivity - Errors",xt="Time (years)",yt="Absolute Error",labels=labels,yarr2=erre) #Error Plot

def RC(b=None,n=None):
    Vo=10
    r=1000
    c=10**-6
    tau=r*c
    f=-y/tau
    t1="RC Circuit - Samarth Jain"
    t2=t1+"\n$\dot V=\dfrac{-V}{RC}=\dfrac{-V}{\\tau}$"
    if((b==None) and (n==None)):
        b=5*tau
        n=10
    h=tau/n
    rc=rk2(a=0,b=b,h=h,yo=Vo,f=f,plot=False,prnt=False)     #Passing h
    rce=elr(a=0,b=b,h=h,yo=Vo,f=f,plot=False,prnt=False)[1] #Passing h
    al=analyt(yo=Vo,x=rc[0],tau=tau)
    errrk=[abs(i-j) for (i,j) in zip(al,rc[1])]
    erre=[abs(i-j) for (i,j) in zip(al,rce)]
    pt=analyt(yo=Vo,x=[tau],tau=tau)
    onept=[tau,pt]
    data=np.column_stack((rc[0],rc[1],rce,errrk,erre))
    head=["Time (Secs.)","V.-RK","V.-Euler","RK Abs. Err","Elr Abs. Err","Vo (V.)","R. (Ohm)","C. (Farads)","Step Size"]
    fr=[Vo,r,c,h]
    table(head,data,fr,t1)
    labels=["Runge-Kutta","Euler"]
    grf(rc[0],rc[1],t2,"Tau (%.4f s.)"%tau,"Volatge (V.)",rce,labels,onept,app=b)
    grf(rc[0],errrk,title="RC - Errors",xt="Time (s)",yt="Absolute Error",labels=labels,yarr2=erre) #Error Plot

def stokes(eta,rad,mass,vo,b=None,n=None):
    tau=mass/(6*np.pi*eta*rad)
    f=-y/tau
    t1="Stoke's Law - Samarth Jain"
    t2=t1+"\n$m\dot v=-6\pi \eta av=\dfrac{-v}{\\tau}$"
    if((b==None) and (n==None)):
        b=5*tau
        n=10
    h=tau/n
    stkrk=rk2(a=0,b=b,h=h,yo=vo,f=f,plot=False,prnt=False)  #Passing h
    stke=elr(a=0,b=b,h=h,yo=vo,f=f,plot=False,prnt=False)[1]#Passing h
    al=analyt(yo=vo,x=stkrk[0],tau=tau)
    errrk=[abs(i-j) for (i,j) in zip(al,stkrk[1])]
    erre=[abs(i-j) for (i,j) in zip(al,stke)]
    pt=analyt(yo=vo,x=[tau],tau=tau)
    onept=[tau,pt]
    data=np.column_stack((stkrk[0],stkrk[1],stke,errrk,erre))
    head=["Time (s.)","Vel.-RK","Vel.-Euler","RK Abs. Err","Elr Abs. Err","Mass (gm.)","vo (m/s)","Radius (m.)","Viscosity (Ns/m)","Step Size"]
    fr=[mass,vo,rad,eta,h]
    table(head,data,fr,t1)
    labels=["Runge-Kutta","Euler"]
    grf(stkrk[0],stkrk[1],t2,"Tau (%.4f s.)"%tau,r'$\mathrm{Velocity\;}\left(\dfrac{m}{s}\right)$',stke,labels,onept,app=b)
    grf(stkrk[0],errrk,title="Stokes - Errors",xt="Time (s)",yt="Absolute Error",labels=labels,yarr2=erre) #Error Plot

if __name__=='__main__':
    os.system('cls')
    f=x+y+x*y           #Example Function
    elr(a=0,b=0.5,h=0.05,yo=1,f=f)
    rk2(a=0,b=0.5,n=10,yo=1,f=f)
    print("Comparison of Euler's and Runge-Kutta's Method:")
    rk_elr_comp(a=0,b=0.5,n=10,yo=1,f=f)
    print("Variation with h")
    hvsi(a=0,b=1,f=f)
    print("Question 3(a):")
    decay()
    print("Question 3(b):")
    RC()
    print("Question 3(b):")
    stokes(eta=20,rad=0.5,mass=100,vo=10)