import math
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import fsolve
import scipy
import pylab

a = -2
b = 3
cnts_n = [5, 10, 15, 20, 30, 40, 50, 60, 70]
cnts_n = [i + 1 for i in cnts_n]
T_CNT = 137
CUR_N = 5 + 1
print(math.e)
def func_f(x):
  #  return (x**2 * math.sin(2*x-3))
    return np.e**(-x) + 5*(x-1)**2

def func_t(j):
    return (a + (b-a)/T_CNT * j)

def func_node(i, N):
    return a + ((b-a) * i) / (N-1)

def getNodes(N):
    arr_N = []
    for i in range(N):
        arr_N.append(func_node(i, N))
    return np.array(arr_N)

def createArray(f, start, end, step):
    x = np.arange(start, end, step)
    y = []
    for x_i in x:
        y_i = f(x_i)
        y.append(y_i)
    return x, np.array(y)

def func_L(t, x, N):
    curL = 0
    for k in range(N):
        P = 1
        for i in range(N):
            if(k != i):
                P *= (t - x[i])/(x[k] - x[i])
        curL += func_f(x[k]) * P
    return curL

def getListL(x, arr_t, N):
    arr_L = []
    arr_f_to_t = []
    for t in arr_t:
        arr_f_to_t.append(func_f(t))
        arr_L.append(func_L(t, x, N))
    return np.array(arr_L), np.array(arr_f_to_t)

def getMaxDist(arr1, arr2):
    res = 0
    for i in range(len(arr1)):
        res = max(abs(arr1[i] - arr2[i]), res)
    return res

def getListR(arr_t):
    list_r = {}
    for N in cnts_n:
        arr_nodes = getNodes(N)
        arr_L, arr_f_to_t = getListL(arr_nodes, arr_t, N)
        list_r[N-1] = getMaxDist(arr_f_to_t, arr_L)
    return (list_r)

def main():
    x, y = createArray(func_f, a, b, 0.01)

    arr_t = [func_t(i) for i in range(0, T_CNT)]
    arr_nodes = getNodes(CUR_N)
    arr_L, arr_f_to_t = getListL(arr_nodes, arr_t, CUR_N)
    list_r_1 = getListR(arr_t)

   
    print("Вычисление погрешности интерполяции с равноотстоящими узлами")
    for r in list_r_1:
        print(f"N={r} -> {list_r_1[r]}")
    print()
    

    list_r_2 = {}
    for N in cnts_n:
        list_v = [math.cos( math.pi  * ( (1 + 2*i)/(2*(N)) ) ) for i in range(N)]
        list_z = [(a+b)/2 + ((b-a)/2)*list_v[i] for i in range(N)]
        
        arr_L_2, arr_f_to_t_2 = getListL(list_z, arr_t, N)
        list_r_2[N-1] = getMaxDist(arr_f_to_t_2, arr_L_2)

    print("Вычисление погрешности интерполяции с Чебышевскими узлами")
    for r in list_r_2:
        print(f"N={r} -> {list_r_2[r]}")

    print()


    l_arr_x = [a + ((b-a)*i)/(CUR_N-1) for i in range(CUR_N)]
    l_y = [func_f(i) for i in l_arr_x]
    l_f = scipy.interpolate.interp1d(l_arr_x, l_y)

    r_3 = []
    for i in arr_t:
        r_3.append(func_f(i) - l_f(i))


    list_r_3 = {}
    for N in cnts_n:
        tmp_arr_x = [a + ((b-a)*i)/(N-1) for i in range(N)]
        tmp_l_y = [func_f(i) for i in tmp_arr_x]
        tmp_l_f = scipy.interpolate.interp1d(tmp_arr_x, tmp_l_y)
        arr_L_2 = [tmp_l_f(x) for x in arr_t]
        arr_f_to_t_2 = [func_f(x) for x in arr_t]
        list_r_3[N-1] = getMaxDist(arr_f_to_t_2, arr_L_2)

    print("Вычисление погрешности при линейной интерполяции")
    for r in list_r_3:
        print(f"N={r} -> {list_r_3[r]}")
    print()





    fig, axes = plt.subplots(ncols=3, nrows=2, sharex=False, sharey=False)
    fig.subplots_adjust(hspace=0.5, wspace=0.2)

    pylab.subplot (2, 3, 1)
    plt.plot(x, y, label="f(x)")
    plt.xlabel("x")
    plt.title('f(x)=y')
    plt.ylabel("y")
    plt.grid(True)
    plt.legend()
        

    pylab.subplot (2, 3, 2)
    plt.plot(arr_t, arr_L, label="L(t[j])")
    plt.plot(arr_t, arr_f_to_t, label="f(t[j])")
    plt.title("f(с равноотст. узлами")
    plt.legend()
    plt.grid(True)

    pylab.subplot (2, 3, 3)
    plt.plot(arr_t, [arr_f_to_t[i] - arr_L[i] for i in range(len(arr_L))])
    plt.title("f(t[j-t[j])")
    plt.grid(True)


    pylab.subplot (2, 3, 4)
    plt.plot(list_r_1.keys(), list_r_1.values(), label="Равноотст. узлами")
    plt.plot(list_r_2.keys(), list_r_2.values(), label="Чебышевскими узлами")
    #plt.plot(list_r_2.keys(), list_r_2.values(), label="Линейная интерпаляция")
    plt.legend()
    plt.yscale('log')
    plt.title('f(t[j-t[j])')
    plt.grid(True)
    

    pylab.subplot (2, 3, 5)
    plt.plot(x, [l_f(i) for i in x])
    plt.title('Линейная интерполяция')
    plt.grid(True)

    pylab.subplot (2, 3, 6)
    plt.plot(arr_t, [abs(i) for i in r_3 ])
    plt.title('f(t[j])-LS(t[j])')
    plt.grid(True)
    
    plt.show()


if __name__ == '__main__':
    main()