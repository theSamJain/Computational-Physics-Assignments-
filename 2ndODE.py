import numpy as np
import matplotlib.pyplot as plt
import math
import table

def func(y, x, Q, consts=None):
    if (Q == 1):
        k, m = consts
        xx, v = y
        dxdt = v
        dydt = -k*xx/m
        eqn = [dxdt, dydt]
    elif (Q == 2):
        k, b, m = consts
        xx, v = y
        dxdt = v
        dydt = (-b*v - k*xx)/m
        eqn = [dxdt, dydt]
    elif (Q == 3):
        g, l = consts
        xx, w = y
        dxdt = w
        dydt = -g*xx/l 
        eqn = [dxdt, dydt]
    elif (Q == 4):
        w0, k, m = consts
        xa, va, xb, vb = y
        dxadt = va
        dxbdt = vb
        dyadt = -(w0**2+k/m)*xa
        dybdt = -(w0**2+k/m)*xb 
        eqn = [dxadt, dyadt, dxbdt, dybdt]
    return eqn

def rk2(y0, x, Q, consts=None, h=None):
    if(h == None):
        h = (x[-1] - x[0]) / len(x)
    yarr = []
    yin = y0
    xin = [x[0]]
    for i in range(len(x)):
        yarr.append(yin)
        k1 = [h * ele for ele in func(yin, xin, Q, consts)]
        yn = [e1 + e2/2 for (e1, e2) in zip(yin, k1)]
        xn = [e1 + h/2 for e1 in xin]
        k2 = [h * ele for ele in func(yn, xn, Q, consts)]
        yn = [e1 + e2/2 for (e1, e2) in zip(yin, k2)]
        yf = [iny + e1 for (iny, e1) in zip(yin, k2)]
        yin = yf
        xin = [e1 + h/2 for e1 in xn]
    yarr=np.array(yarr).reshape(-1, len(yin))
    return(yarr)

def rk4(y0, x, Q, consts=None, h=None):
    if(h == None):
        h = (x[-1] - x[0]) / len(x)
    yarr = []
    yin = y0
    xin = [x[0]]
    for i in range(len(x)):
        yarr.append(yin)
        k1 = [h * ele for ele in func(yin, xin, Q, consts)]
        yn = [e1 + e2/2 for (e1, e2) in zip(yin, k1)]
        xn = [e1 + h/2 for e1 in xin]
        k2 = [h * ele for ele in func(yn, xn, Q, consts)]
        yn = [e1 + e2/2 for (e1, e2) in zip(yin, k2)]
        k3 = [h * ele for ele in func(yn, xn, Q, consts)]
        yn = [e1 + e2 for (e1, e2) in zip(yin, k3)]
        xn = [e1 + h for e1 in xin]
        k4= [h * ele for ele in func(yn, xn, Q, consts)]
        yf = [ini_y + (e1 + 2 * (e2 + e3) + e4) / 6 for (ini_y,e1,e2,e3,e4) in zip(yin, k1, k2, k3, k4)]
        yin = yf
        xin = [e1 + h/2 for e1 in xn]
    yarr=np.array(yarr).reshape(-1, len(yin))
    return(yarr)

def plot(t, solrk2, solrk4, label, title, yl, tp, tl):
    tks = np.arange(0, t[-1], tp)
    fig, ax = plt.subplots(1, 2)
    ax[0].plot(t, solrk2[:, 0], marker = "*", label = label[0])
    ax[0].plot(t, solrk4[:, 0], marker = ".", label = label[1])
    ax[0].set(xticks = tks, xlabel = "Time Periods", ylabel = yl[0], title = title[0])
    ax[0].set_xticklabels(tl)
    ax[0].legend(loc = 'upper left')
    ax[0].grid()
    ax[1].plot(t, solrk2[:, 1], marker = "*", label = label[2])
    ax[1].plot(t, solrk4[:, 1], marker = ".", label = label[3])
    ax[1].set(xticks = tks, xlabel = "Time Periods", ylabel = yl[1], title = title[0])
    ax[1].set_xticklabels(tl)
    ax[1].set(title = title[1])
    ax[1].legend(loc = 'upper left')
    ax[1].grid()
    fig.suptitle(title[2])
    fig.tight_layout()
    plt.show()

t = np.linspace(0, 50, 200)

consts1 = [2, 1]
y01 = [1, 0]
sol02 = rk2(y0 = y01, x = t, consts = consts1, Q = 1)
sol04 = rk4(y0 = y01, x = t, Q = 1, consts = consts1)
tp = 2*np.pi*(consts1[1]/consts1[0])**0.5
tl = [str(i) for i in range(math.ceil(t[-1]/tp))]
head = ["x-RK2", "x-RK4", "% Error","v-RK2", "v-RK4", "% Error", "k", "m"]
err1 = [(i-j)/j*100 for (i, j) in zip(sol02[:, 0], sol04[:, 0])]
err2 = [(i-j)/j*100 for (i, j) in zip(sol02[:, 1], sol04[:, 1])]
data = np.column_stack((sol02[:, 0], sol04[:, 0], err1, sol02[:, 1], sol04[:, 1], err2))
title1 = ["Position - SHM", "Velocity - SHM", "Simple Harmonic Motion"]
yl = ["Position $x$", "Velocity $v$"]
label1 = ["x -> RK2", "x -> RK4", "v -> RK2", "v -> RK4"]
table.table(head, data, splfr = consts1, title = title1[2])
plot(t, sol02, sol04, label = label1, title = title1, yl = yl, tp = tp, tl = tl)
plt.plot(t, err1, marker = ".", label = "x %Error")
plt.plot(t, err2, marker = ".", label = "v %Error")
plt.title("Error Plot")
plt.xlabel("Time")

consts2 = [0.5, 0.1, 1]
y2 = [1, 0]
sol12 = rk2(y0 = y2, x = t, consts = consts2, Q = 2)
sol14 = rk4(y0 = y2, x = t, Q = 2, consts = consts2)
tp = 2*np.pi* 1/( consts2[0]/consts2[2]-(consts2[1]/(2*consts2[0]))**2 )**0.5
tl = [str(i) for i in range(math.ceil(t[-1]/tp))]
head = ["x-RK2", "x-RK4", "% Error","v-RK2", "v-RK4", "% Error", "k", "b", "m"]
err1 = [(i-j)/j*100 for (i, j) in zip(sol12[:, 0], sol14[:, 0])]
err2 = [(i-j)/j*100 for (i, j) in zip(sol12[:, 1], sol14[:, 1])]
data = np.column_stack((sol02[:, 0], sol04[:, 0], err1, sol02[:, 1], sol04[:, 1], err2))
title2 = ["Position - Damped SHM", "Velocity - Damped SHM", "Damped Simple Harmonic Motion"]
yl = ["Position $x$", "Velocity $v$"]
label2 = ["x -> RK2", "x -> RK4", "v -> RK2", "v -> RK4"]
table.table(head, data, splfr = consts2, title = title2[2])
plot(t, sol12, sol14, label = label2, title = title2, yl = yl, tp = tp, tl = tl)
plt.plot(t, err1, marker = ".", label = "x %Error")
plt.plot(t, err2, marker = ".", label = "v %Error")
plt.title("Error Plot")
plt.xlabel("Time")

consts3 = [9.81/6, 1]
y2 = [np.pi/6, 0]
sol22 = rk2(y0 = y2, x = t, consts = consts3, Q = 3)
sol24 = rk4(y0 = y2, x = t, Q = 3, consts = consts3)
tp = 2*np.pi*(consts3[1]/consts3[0])**0.5
tl = [str(i) for i in range(math.ceil(t[-1]/tp))]
head = ["\u0398-RK2", "\u0398-RK4", "% Error","\u03c9-RK2", "\u03c9-RK4", "% Error", "g", "L"]
err1 = [(i-j)/j*100 for (i, j) in zip(sol22[:, 0], sol24[:, 0])]
err2 = [(i-j)/j*100 for (i, j) in zip(sol22[:, 1], sol24[:, 1])]
data = np.column_stack((sol22[:, 0], sol24[:, 0], err1, sol22[:, 1], sol24[:, 1], err2))
title3 = ["Ang-Position - Simple Pendulum", "Ang-Velocity - Simple Pendulum", "Simple Pendulum on Moon"]
yl = ["Angle $\\theta$", "Angular Velocity $\omega$"]
label3 = ["$\omega$ -> RK2", "$\omega$ -> RK4", "$\\theta$ -> RK2", "$\\theta$ -> RK4"]
table.table(head, data, splfr = consts3, title = title3[2])
plot(t, sol22, sol24, label = label3, title = title3, yl = yl, tp = tp, tl = tl)
plt.plot(t, err1, marker = ".", label = "x %Error")
plt.plot(t, err2, marker = ".", label = "v %Error")
plt.title("Error Plot")
plt.xlabel("Time")

consts4 = [2*np.pi/5, 1, 1]
y3 = [1, 0, 2, 0]
sol32 = rk2(y0 = y3, x = t, consts = consts4, Q = 4)
sol34 = rk4(y0 = y3, x = t, Q = 4, consts = consts4)
tp = 2*np.pi*1/(consts4[0]**2 + consts4[1]/consts4[2])**0.5
tl = [str(i) for i in range(math.ceil(t[-1]/tp))]
head = ["xA-RK2", "xA-RK4", "% Error", "vA-RK2", "vA-RK4", "% Error", "xB-RK2", "xB-RK4", "% Error", "vB-RK2", "vB-RK4", "% Error" "w0", "k", "m"]
err1 = [(i-j)/j*100 for (i, j) in zip(sol32[:, 0], sol34[:, 0])]
err2 = [(i-j)/j*100 for (i, j) in zip(sol32[:, 1], sol34[:, 1])]
err3 = [(i-j)/j*100 for (i, j) in zip(sol32[:, 2], sol34[:, 2])]
err4 = [(i-j)/j*100 for (i, j) in zip(sol32[:, 3], sol34[:, 3])]
data = np.column_stack((sol32[:, 0], sol34[:, 0], err1, sol32[:, 1], sol34[:, 1], err2, sol32[:, 2], sol34[:, 2], err3, sol32[:, 3], sol34[:, 3], err4))
title41 = ["$x_A$", "$v_A$", "Coupled Oscillations"]
title42 = ["$x_B$", "$v_B$", "Coupled Oscillations"]
yl1 = ["$x_A$", "$v_A$"]
yl2 = ["$x_B$", "$v_B$"]
label31 = ["$x_A$ -> RK2", "$x_A$ -> RK4", "$v_A$ -> RK2", "$v_A$ -> RK4"]
label32 = ["$x_B$ -> RK2", "$x_B$ -> RK4", "$v_B$ -> RK2", "$v_B$ -> RK4"]

arr1 = np.array([ [i, j] for (i, j) in zip(sol32[:, 2], sol32[:, 3])]).reshape(-1,2)
arr2 = np.array([ [i, j] for (i, j) in zip(sol34[:, 2], sol34[:, 3])]).reshape(-1,2)
table.table(head, data, splfr = consts4, title = title41[2])
plot(t, sol32, sol34, label = label31, title = title41, yl = yl1, tp = tp, tl = tl)
plot(t, solrk2 = arr1, solrk4 = arr2, label = label32, title = title42, yl = yl2, tp = tp, tl = tl)
plt.plot(t, err1, marker = ".", label = "xA %Error")
plt.plot(t, err2, marker = "*", label = "vA %Error")
plt.plot(t, err3, marker = "p", label = "xB %Error")
plt.plot(t, err4, marker = "x", label = "vB %Error")
plt.title("Error Plot")
plt.xlabel("Time")

plt.legend()
plt.show()