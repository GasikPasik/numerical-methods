
from math import sin, pi, e
from scipy import integrate


def f(x):
    # return sin(x**2 -1) ** 2
    return e**(-x) + 5*(x - 1)**2


def main():
    a = -1
    a = -2
    
    b = pi/2
    b = 3 

    CONST_K = 10

    resI, err = integrate.quad(f, a, b)

    print(f"Решение интеграла встроенными методами: {resI}")

    n = [ 2**i for i in range(1, CONST_K+1) ]

    h = [(b-a)/nk for nk in n]

    j = [] 
    for k in range(CONST_K):
        j.append(0)
        for i in range(n[k]): # - 1
            j[k] += (h[k]/2) * (f(a + i*h[k]) + f(a + h[k] + i*h[k]))
    
    z = [abs(resI - jk) for jk in j]

    print("\nПогрешности: ")
    for k in range(CONST_K):
        print(f"{k}: При nk = {n[k]}: {z[k]}")


    r = [ (abs(j[k] - j[k-1]))/3 for k in range(1, CONST_K)]

    print("\nПогрешность по методу Рунге:")
    for k in range(CONST_K-1):
        print(f"{k}: При nk = {n[k]}: {r[k]}")

    print()
    for k in range(1, CONST_K-1):
        print(f"Сотношение r{k-1} и r{k}: {r[k-1]/r[k]}")

    

    x_list = {
        4: [0.06943184,
         0.33000948,
         0.66999052,
         0.93056815],

        5:[0.04691008,
           0.23076534,
           0.5,
           0.76923466,
           0.95308992],

        6: [0.03376524,
        0.16939531,
        0.38069041,
        0.61930959,
        0.83060469,
        0.96623475]}
    
    c_list = {
        4: [0.17392742,
         0.32607258,
         0.32607258,
         0.17392742],

        5:[0.11846344,
           0.23931433,
           0.28444444,
           0.23931433,
           0.11846344],

        6: [0.08566225,
         0.18038079,
         0.23395697,
         0.23395697,
         0.18038079,
         0.08566225]}


    for cnt_n in x_list.keys(): 
        x = x_list[cnt_n]
        C = c_list[cnt_n]
        m = len(x)

        z = [ a + (b-a)*x[j] for j in range(m)]

        G = (b-a) * sum([ C[j]*f(z[j]) for j in range(1, m) ])
        print(f"\nПри len(x) = {cnt_n}")
        print(f"{G=}")        
        print(f"|I - G| = {abs(resI-G)}")
        breakpoint()


if __name__ == "__main__":
    main()