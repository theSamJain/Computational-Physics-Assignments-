# To write functions that generate Taylor's Series for e^x, sin(x) and cos(x)

import math as m
from sympy.abc import x

def Exp(n):
    expo = 1
    for i in range(1, n+1):
        expo += (x**i)/m.factorial(i)
    print("exp(x, %s) = "%n, expo)

def Sin(n):
    sine = x
    for i in range(3, n+2, 2):
        sine += ((-1)**((i-1)/2))*(x**i)/(m.factorial(i))
    print("sin(x, %s) = "%n, sine)

def Cosine(n):
    cosine = 1
    for i in range(2, n+1, 2):
        cosine += ((-1)**(i/2))*(x**i)/(m.factorial(i))
    print("cos(x, %s) = "%n,cosine)

if __name__ == '__main__':
    n = int(input("Enter the number of terms (n): "))
    Exp(n)
    Sin(n)
    Cosine(n)
