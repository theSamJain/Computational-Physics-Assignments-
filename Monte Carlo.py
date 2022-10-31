import random
import numpy as np
import matplotlib.pyplot as plt

counter = 0 
# rad = float(input("Enter Radius = "))
rad = 0.5
n = int(input("n = "))
x_sqr, y_sqr = np.linspace(-rad, rad, n), np.linspace(-rad, rad, n)

xarr, yarr = [], []
for i in range(n):
    x, y = random.uniform(-rad, rad), random.uniform(-rad, rad)
    if( (x**2 + y**2) <= rad**2 ):
        xarr.append(x)
        yarr.append(y)
        counter += 1

print("counter = ", counter) 
pi = 4*counter/n
print("pi = ", pi)

sqr = plt.Rectangle((-rad, -rad), 2*rad, 2*rad, fc = "white", ec = "red")
circ = plt.Circle((0, 0), rad, fc = "white", ec = "orange")
plt.gca().add_patch(sqr)
plt.gca().add_patch(circ)
plt.axis('scaled')
plt.scatter(xarr, yarr, marker = '.', zorder = 10)
plt.show()