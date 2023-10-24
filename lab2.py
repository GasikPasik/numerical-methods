import numpy as np


import matplotlib.pyplot as plt
from scipy.optimize import fsolve
import pylab

list_r = 0


def solve_system_by_inverse(A, b):
    x = np.linalg.solve(A, b)
    return x


def solve_system_by_iterative_method(A, b, tol=1e-10):
    x = np.zeros(A.shape[1])
    r = b - A @ x

    while np.linalg.norm(r) > tol:
        J = np.dot(r, A @ r) / np.dot(A @ r, A @ r)
        x = x + r / J
        r = r - J * r

    return x

def MatrixToVec(A, vec):
    res = np.array([0.0]*len(vec))
    tmp = A*vec
    
    for i in range(len(res)):
        res[i] = sum(tmp[i])

  
    return (res)

def MMR(A, b, max_iter = 1000, tol=1e-10):
    global list_r
    list_x = []
    list_x.append(np.array([0]*len(b)))

    list_r = []
    list_r.append(MatrixToVec(A, list_x[0]) - b)
 
    tmp = (MatrixToVec(A, list_r[0]))
    t = np.dot( list_r[0], tmp) /  np.dot( tmp, tmp )
    
    list_x.append(list_x[0] - t*(MatrixToVec(A, list_x[0]) - b))
    E = (np.eye(A.shape[0]))
    for i in range(1, max_iter):
        tmp = (MatrixToVec(A, list_r[-1]))
        t = np.dot( list_r[-1], tmp) /  np.dot( tmp, tmp )
        list_x.append(list_x[-1] - t*(MatrixToVec(A, list_x[-1]) - b))
        list_r.append(MatrixToVec(A, list_x[-1]) - b)
    
    norm = np.linalg.norm(E - t*A, np.inf)
    print("\nНорма матрицы, s: " + str(norm))
    if(norm <= 1):
        print("Условие сходимости соблюдены!")
    else:
        print("Условие сходимости не соблюдены!")
    print()

    return list_x

def printMatrixTwo(A, B, nameA, nameB):
    maxLenA = 0
    for i in A:
        for j in i:
            maxLenA = max(maxLenA, len(str(round(j, 3))))
    maxLenB = 0
    for i in B:
        for j in i:
            maxLenB = max(maxLenB, len(str(round(j, 3))))
  
    print(f" {nameA}:" + " "*(maxLenA*len(A[0]) - len(nameA)  + (len(A[0]) + 1)) + f"  {nameB}:")
    for i in range(len(A)):
        print("|", end="")
        for j in A[i]:
            print(str(j) + " "*(maxLenA-len(str(j))) + "|", end = "")
        print("\t|", end="")
        for j in B[i]:
            print(str(round(j, 2)) + " "*(maxLenB-len(str(round(j, 2)))) + "|", end = "")
        print()
    print()


def printMatrixToSystem(A, b):
    print("System: ")
    for i in range(len(A)):
        print(f"\t{i}: x({1})*{A[i][0]}", end="")
        for j in range(1, len(A[i])):
            print(f" + x({j+1})*{A[i][j]}", end="")
        print(f" = {b[i]}")

def solve(A, b):

    inverse_A = np.linalg.inv(A)
    printMatrixTwo(A, inverse_A, "A", "A^-1")
    print(" b:")
    for i in b:
        print(f"|{i}" + " "*(7-len(str(i))) + "|")
    print("")

    printMatrixToSystem(A, b)

    list_x = MMR(A, b, 10)
    x = list_x[-1]

    print("Solve: ")
    for i in range(len(A)):
        print(f"\t{i}: {round(x[0], 5)}*{A[i][0]}", end="")
        for j in range(1, len(A[i])):
            print(f" + {round(x[j], 5)}*{A[i][j]}", end="")
        print(f" = {b[i]}")
    print()

    for i in range(len(A)):
        print(f"\t{i}: {round(x[0]*A[i][0], 5)}", end="")
        sum = x[0]*A[i][0] 
        for j in range(1, len(A[i])):
            print(f" + {round(x[j]*A[i][j], 5)}", end="")
            sum += x[j]*A[i][j]
        print(f" = {sum}({round(sum, 3)})")
    print()

    for i in range(len(x)):
        print(f"x{i}={round(x[i],3)}, ", end = "")
    print()

    root_x = solve_system_by_inverse(A, b)
    print("\nSolve by standart func: \n\t", end="")
    for i in range(len(x)):
        print(f"x{i}={round(x[i],3)}, ", end = "")
    print("\n\n\n")

    pogreshnost = [abs(np.linalg.norm(i - list_x[-1])) for i in list_x]
    print(pogreshnost)
    plt.plot(range(0, len(pogreshnost)), pogreshnost)
    plt.xlabel('iteration')
    plt.ylabel('Погрешности')
    plt.yscale('log')
    plt.title('Погрешность овтетов')
    plt.grid(True)
    plt.show()




A = np.array([[1, 2], [3, 4]])

def main(k):
    
    D = ([[1.342, 0.432, -0.599, 0.202],
                 [0.202, 1.342, 0.432, -0.599], 
                 [-0.599, 0.202, 1.342, 0.432],
                 [0.432, -0.599, 0.202, 1.342]])
    C = np.array([[0.02, 0, 0, 0],
                 [0, 0.02, 0, 0],
                 [0, 0, 0.02, 0],
                 [0, 0, 0, 0.02]])

    A = D + k*C
    b = np.array([1.941, -0.23, -1.941, 0.23])
    print("k: " + str(k) + "\n")
    solve(A,b)

    


main(21)