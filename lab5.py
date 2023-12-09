import math
import matplotlib.pyplot as plt
import numpy as np
import pylab
from numpy.polynomial import Polynomial


# slava
# list_x = [0.25, 0.31, 0.36, 0.39, 0.43, 0.47, 0.52, 0.56, 0.64, 0.66, 0.71]
# list_y = [1.031, 1.048, 1.066, 1.107, 1.194, 1.233, 1.218, 1.161, 1.121, 1.102, 1.056] 

# simple
list_x = [0.2, 0.31, 0.36, 0.39, 0.43, 0.47, 0.52, 0.56, 0.62, 0.66, 0.70]
list_y = [0.2474, 0.2951, 0.3623, 0.3902, 0.4069, 0.3529, 0.3226, 0.3112, 0.3772, 0.4931, 0.6018] 
# moe
list_x = [1, 1.3, 1.4, 1.5, 1.7, 2.0, 2.3, 2.5, 2.8, 3.1, 3.3]
list_y = [0.179, 0.114, 0.106, 0.124, 0.130, 0.178, 0.222, 0.242, 0.295, 0.310, 0.331]

n = len(list_x)

if(len(list_x) != len(list_y)):
    print(f"Error, len(x):{len(list_x)} != len(y):{len(list_y)}")
    exit(0)

def create_matrix_k(vec, k):
    result_arr = [[] for _ in range(k)]

    for i in range(k*k):
        result_arr[i // k].append(vec[i // k + i % k])

    return np.array(result_arr)

def func_k(k, s, b):
    matrix_A = create_matrix_k(s, k)
    matrix_B = np.array([b[i] for i in range(k)])

    solve = np.linalg.solve(matrix_A, matrix_B)

    p = []
    
    for x in list_x:
        tmp = 0
        for i in range(k):
            tmp += solve[i] * x**i
        p.append(tmp)    
    
    return np.array(p)

def find_max_diff(p):
    max_diff = 0   

    for i in range(n):
        max_diff = max(max_diff, abs(list_y[i] - p[i]))
    
    return max_diff


def find_diff(p):
    list_diff = []   

    for i in range(n):
        list_diff.append(abs(list_y[i] - p[i]))
    
    standart_diff = (1/(n)) * math.sqrt(sum([i**2 for i in list_diff]))

    return list_diff, standart_diff


        
def main():
    print("Start working..")

    list_s = []
    for j in range(n):
        tmp = 0
        for x in list_x:
            tmp += x**j
        list_s.append(tmp)
    
    list_b = []
    for j in range(n):
        tmp = 0
        for i in range(n):
            tmp += (list_y[i] * list_x[i]**j)
        list_b.append(tmp)


    p2 = func_k(2 + 1, list_s, list_b)
    p3 = func_k(3 + 1, list_s, list_b)

    approximation = Polynomial.fit(list_x, list_y, 8)
    p8 = np.array([approximation(x) for x in list_x])
    
    tmp_dict = {"p2": p2, "p3": p3, "p8": p8}
    dict_diff = {}


    for i in tmp_dict.keys():
        list_diff, standart_diff = find_diff(tmp_dict[i])
        print(f"Абсолют погреш. у {i} = {round(max(list_diff), 6)}, а среднеквадрат. погреш. = {round(standart_diff, 6)}")

        dict_diff[i] = list_diff

    mina = [0, 1000]


    for i in range(n):
        approx = Polynomial.fit(list_x, list_y, i)
        p = np.array([approx(x) for x in list_x])
        list_diff, standart_diff = find_diff(p)
        if(mina[1] > standart_diff):
            mina = [i, standart_diff]
        print(f"При m={i}, среднеквадрат. погрешность будет: {standart_diff}")

    print(f"Минимальная погрешность была при: {mina[0]}")

    # graphics
    fig, axes = plt.subplots(ncols=4, nrows=1, sharex=False, sharey=False)
    fig.subplots_adjust(hspace=0.5, wspace=0.2)

    pylab.subplot (1, 4, 1)
    plt.plot(list_x, list_y, label="f(x)")
    plt.xlabel("x")
    plt.title('f(x)=y')
    plt.ylabel("y")
    plt.grid(True)
    plt.legend()

    pylab.subplot (1, 4, 2)
    for i in tmp_dict.keys():
        plt.plot(list_x, tmp_dict[i], label=f"{i}(t)")
    plt.xlabel("x, t")
    plt.ylabel("y")
    plt.grid(True)
    plt.legend()

    pylab.subplot (1, 4, 3)
    for i in tmp_dict.keys():
        plt.plot(list_x, dict_diff[i], label=f"{i}(t)")
    plt.xlabel("x, t")
    plt.ylabel("y")
    plt.grid(True)
    plt.legend()

    points = 1000
    print(f"Точек: {points}")

    pylab.subplot (1, 4, 4)
    for i in [2,3,6,8,10]:
        new_list_x = [ i/1000 for i in range(1000, 3000)]
        approx = Polynomial.fit(list_x, list_y, i)
        new_list_y = [approx(x) for x in new_list_x]
        plt.plot(new_list_x, new_list_y, label=f"{i}(t)")



    plt.xlabel("x, t")
    plt.ylabel("y")
    plt.grid(True)
    plt.legend()


    plt.show()

if __name__ == "__main__":
    main()