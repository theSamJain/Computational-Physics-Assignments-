import matplotlib.pyplot as pt
import matplotlib.ticker as ticker
import numpy as np
from sklearn.linear_model import LinearRegression as lr
import pandas as pd
from matplotlib.ticker import (AutoMinorLocator, MultipleLocator)

# Arrays
x = np.array([233266189.9,299682091.6,335317019.6,409652196.1,593570126]).reshape(-1,1)
y = np.array([1.515699465,1.520245036,1.522510879,1.533210179,1.5360083]).reshape(-1,1)

# Fitting Params 
    # Nomal Fitting ie fit xy
model = lr().fit(x,y)
slope1 = float(model.coef_)
intc1 = model.intercept_

    #Fit YX for Cauchy 
model2 = lr().fit(y,x)
slope2 = 1/float(model2.coef_)
intc2 = model2.intercept_*(-1)/float(model2.coef_)

    #Best line ie average of fit xy and fit yx
slopef = (slope1+slope2)/2
intcf = (intc1+intc2)/2


#GRAPHING EQUATIONS

xx=np.linspace(min(x)-1,max(x)+1,100)
ax=pt.axes()
ax.set_facecolor('black')

    #Line of fit xy
y2=slope1*xx+intc1

    #Line of fit yx
y3=slope2*xx+intc2 #for cauchy ie fit yx

    #Average line of Fit xy and Fit yx
yav=slopef*xx+intcf #average line of fit xy and fit yx for cauchy

    #Uncertainties for Cauchy
yc=[x*slope1+intc1 for x in x]
unc=[]
for i in range(len(y)):
    unc.append((float(yc[i]-y[i])))
#print("UNC",unc)
weights=[1/unct**2 for unct in unc]

    #Weighted lsf
wls=lr().fit(x,y,sample_weight=weights)
wlsslope=wls.coef_
wlsintc=wls.intercept_
ywls=wlsslope*xx+wlsintc


#FUNCTIONS FOR CHECKING CORRECTNESS
def sr(x,y,m,c):                                #sum of residuals
    Y=np.array([])
    i=0
    sd=0
    while(i<len(x)):
        z=m*x[i]+c
        d=(y[i]-z)
        sd+=d
        Y=np.append(Y,z).reshape(-1,1)
        i+=1
    return(sd)
    
def ssr(x,y,m,c):                               #sum of residual squared
    i=0
    s=0
    while(i<len(x)):
        z=m*x[i]+c
        ssq=(y[i]-z)**2
        s+=ssq
        i+=1
    return(s)

def ess(x,y):                               #Explained sum squares
    i=0
    t=0
    sum2=0
    se=0
    while(i<len(x)):
        sum2+=x[i]
        i+=1
    m=sum2/len(x)
    while(t<len(x)):
        se+=(x[t]-m)**2
        t+=1
    ess=se**0.5
    return(ess)

def sumxsq(x):
    i=0
    sum=0
    while(i<len(x)):
        sum+=(x[i])**2
        i+=1
    summ=sum**0.5
    return(summ)

#CALCULATIONS OF VARIOUS QUANTITIES FOR FIT XY
a=(ssr(x,y,wlsslope,wlsintc))/(len(x)-2)
A=a**0.5
B=ess(x,y)
err=A/B
r_sq=model.score(x,y)
errintc=A*(sumxsq(x))/((len(x)**0.5)*B)


#CALCULATING CHI sq
yerr=[]
yf=[x*wlsslope+wlsintc for x in x]
for i in range(len(y)):
    yerr.append((float(yf[i]-y[i])))
sumy=0
for j in range(len(yerr)):
    sumy+=float(yerr[j]**2)*weights[j]
chi2=sumy


#CALCULATING >1sig POINTS
ind=[]
sig=np.std(y)
errbars=[]
for i in range(len(x)):
    e=(float(errintc)**2+float(err/x[i])**2)**0.5
    errbars.append(float(e))
#print(errbars)
for i in range(len(errbars)):
    if (errbars[i]>sig):
        ind.append(i)


#DISPERSIVE POWER AND ERROR
dsp=(yf[4]-yf[0])/((yf[4]+yf[0])/2-1) #num/denom
dsperr=2*((((yf[4]+yf[0])/2-1) - (yf[4]-yf[0]))*(yerr[4]) - (((yf[4]+yf[0])/2-1) + (yf[4]-yf[0]))*(yerr[0]))/(((yf[4]+yf[0])/2-1)**2)
#formula for error = 2*((denom - num)(del wavel 2) - (denom+num)(del wavel 1))/denom**2


#PRINTING QUANTITIES DEMONSTRATING CORRECTNESS OF THE FIT XY

print("\nSlope of WLS Fitted Line\t\t\t=\t{:.1e}\nIntercept of WLS Fitted Line\t\t\t=\t{:.4}\nCoefficient of Determination of WLS Fitted Line\t=\t{:.4}\nCorrelation Coefficient of WLS Fitted Line\t=\t{:.4}\nStandard Error of Slope in WLS Fitted Line\t=\t{:.1e}\nStandard Error of Intercept in WLS Fitted Line\t=\t{:.0e}\nChi Square\t\t\t\t\t=\t{:.5}\n".format(float(wlsslope), float(wlsintc), float(wls.score(x,y)),  float(wls.score(x,y)**0.5), float(err), float(errintc), float(chi2)))
'''
print("\n{:} points have deviation larger than 1 Sigma with respect to the fitted line, where Sigma = {:.5}".format(int(len(ind)), float(sig)))
for i in range(len(ind)):
    print("{:}. ({:.e},{:})\n".format(int(i), float(x[ind[i]]), float(y[ind[i]])))

print("\n\nCauchy's Constants are as follows:")
print("B = {:.2} +/- {:.0} cm^2 and A = {:.4} +/- {:.0} ".format(float(wlsslope), float(err), float(wlsintc), float(errintc)))
print("\n\nDispersive Power = {:.2e}\nError in Dispersive Power = {:.0e}\n".format(float(dsp), float(dsperr)))


#PLOT
    #Plots and Scatters
pt.scatter(x,y,marker='o', color="xkcd:sky blue")
pt.plot(xx,y2,c="gold",linestyle='--', linewidth=1, label="OLS Fitted Line")
#for cauchy
#pt.plot(xx,y3,c='xkcd:mint green',linestyle='--',linewidth=1, label="Fitted wrt to $\dfrac{1}{\lambda^2}$")
#pt.plot(xx,yav,c="orange",linewidth=1, label="Best Fit")
pt.plot(xx,ywls,c="pink",linewidth=1,label="WLS Fitted Line")

    #Errorbars 
pt.errorbar(x,y,yerr=errbars,c="xkcd:sky blue",fmt="o")

    #Grids
pt.grid(True)
#ax.grid(which='both')
ax.grid(which='both',alpha=0.2)
ax.grid(which='major',alpha=0.5)
#pt.xticks(np.arange(min(x), max(x)+1, 2.0))

    #Title
pt.title("Team Number 20\nCauchy's Constant\n $\mu$ v. $\dfrac{1}{\lambda^2}$")       #NAME OF EXPERIMENT
pt.ylabel("Refractive Index. ($\mu$)\n")                    #Y LABEL
pt.xlabel("Wavelength ($\dfrac{1}{\lambda^2}$) (cm.)\n")    #X LABEL

pt.legend()
pt.show()
'''