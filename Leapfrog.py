import numpy as np
import matplotlib.pyplot as plt
import os
import time

begin = time.time()

def leapfrog(incd, t, consts = None, CO = None):
    G, M = consts
    dt = (t[1] - t[0])

    dirname = "{:} - ({:.3e}), ({:3e}), ({:5e})".format(str(CO), (incd[0][0]**2 + incd[0][1]**2)**0.5, (incd[1][0]**2+incd[1][1]**2)**0.5, dt)
    os.makedirs(r'C:\Users\jains\Desktop\Voyager-Trajectory\Samarth Jain\%s '%dirname, exist_ok=False)

    # Note: Do Not Change the Order!
    # Initialise
    x = np.zeros(len(t))    ; x[0] = incd[0][0]
    vx = np.zeros(len(t))   ; vx[0] = incd[1][0]
    y = np.zeros(len(t))    ; y[0] = incd[0][1]
    vy = np.zeros(len(t))   ; vy[0] = incd[1][1]
    r = np.zeros(len(t))    ; r[0] = (x[0]**2+y[0]**2)**0.5
    ax = np.zeros(len(t))   ; ax[0] = G*M*x[0]/(r[0]**3)
    ay = np.zeros(len(t))   ; ay[0] = G*M*y[0]/(r[0]**3)

    for i in range(len(t)-1):
        # Velocity Vectors
        vx[i + 1] = vx[i] -  G*M*x[i]*dt/(r[i]**3)
        vy[i + 1] = vy[i] -  G*M*y[i]/(r[i]**3)*dt
        # Position Vectors
        x[i + 1] = x[i] + dt*vx[i + 1]
        y[i + 1] = y[i] + dt*vy[i + 1]
        # Radius
        r[i+1] = (x[i+1]**2+y[i+1]**2)**0.5
        # Acceleration
        ax[i+1] = G*M*x[i+1]/(r[i+1]**3)
        ay[i+1] = G*M*y[i+1]/(r[i+1]**3)

    data = np.column_stack((t, x, y, vx, vy, ax, ay))
    np.savetxt(r'C:\Users\jains\Desktop\Voyager-Trajectory\Samarth Jain\%s\%s.txt'%(dirname, CO), data, header= "time xpos ypos velx vely accx accy")

    return(data)

if __name__ == '__main__':
    year = 365*24*60*60     # seconds
    G = 6.6743e-11
    M = 1.9891e30
    consts = [G,M]

    t = np.linspace(0,year,int(year/60))
    in_posE = [1.47098e11, 0]
    in_velE = [0, 30288]

    # in_posJ = [740.52e9, 0]
    # in_posJ = []

    incdE = [in_posE, in_velE]
    leapfrog(incdE, t, consts, 'Earth')
    time.sleep(0)
    end = time.time()
    print(end-begin)