import numpy as np

def rk4(func, y0, a, b, n):
    h = (b - a) / n
    yarr = []
    yin = y0
    x_axis = np.linspace(a, b, n+1)
    for i in range(n):
        yarr.append(yin)
        k1 = [h * ele for ele in func(yin, x_axis[i])]
        yn = [e1 + e2/2 for (e1, e2) in zip(yin, k1)]
        k2 = [h * ele for ele in func(yn, x_axis[i]+h/2)]
        yn = [e1 + e2/2 for (e1, e2) in zip(yin, k2)]
        k3 = [h * ele for ele in func(yn, x_axis[i]+h/2)]
        yn = [e1 + e2 for (e1, e2) in zip(yin, k3)]
        k4 = [h * ele for ele in func(yn, x_axis[i]+h)]
        yf = [ini_y + (e1 + 2 * (e2 + e3) + e4) / 6 for (ini_y,e1,e2,e3,e4) in zip(yin, k1, k2, k3, k4)]
        yin = yf
    yarr = np.array(yarr).reshape(-1, len(yin))

    res_matrix = []
    for i in range(yarr.shape[1]):
        y_i = list(yarr[:, i:i+1].flatten())
        res_matrix.append(y_i)
    return(res_matrix)