
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
import pylab

def f(x):
    return np.exp(-x) - 5*(x-1)**2

def derivative_of_f(x):
    return -np.exp(-x) - 10*(x-1)

def simple_iteration_method(g, x_start, epsilon, max_iterations):
    x = x_start
    x_arr = [x]

    for i in range(max_iterations):
        next_x = g(x)
        x_arr.append(next_x)

        if abs(next_x - x) < epsilon:
            return (next_x, x_arr)

        x = next_x

    print("Метод не сошелся за заданное количество итераций.")
    return (None, x_arr)

def g(x, m = 1):
    return x - f(x) / derivative_of_f(x)


x = np.linspace(-2, 2, 100)
y = f(x)

epsilon = 1e-10
root_x0 = fsolve(f, [0, 5])
print(f"root(f(x) {root_x0}")

x_start_1 = 1.5
root_1, x_arr_1 = simple_iteration_method(g, x_start_1, epsilon, 50)
print(f"Приближенное решение уравнения (старт: {x_start_1}) f(x) = 0: {root_1}")
print(x_arr_1)

x_start_2 = 0.2
root_2, x_arr_2 = simple_iteration_method(g, x_start_2, epsilon, 100)
print(f"Приближенное решение уравнения (старт: {x_start_2}): {root_2}")
print(x_arr_2)




if(root_1 != None and root_2 != None):
        difference = abs(root_1 - root_2)
        print("Модуль разности приближенных решений:", difference)

print(f"Погрешность root_1: {abs(g(x_arr_1[-1]) - x_arr_1[-1])}")
print(f"Погрешность root_2: {abs(g(x_arr_2[-1]) - x_arr_2[-1])}")



fig, axes = plt.subplots(nrows=4, ncols=4, figsize=(4, 4))
fig.tight_layout()

pylab.subplot (2, 2, 1)
plt.plot(x, y)
if(root_1 != None):
    plt.scatter([root_1], [f(root_1)], color='red')
if(root_2 != None):
    plt.scatter([root_2], [f(root_2)], color='red')
plt.xlabel('x')
plt.ylabel('f(x)')
plt.title('f(x) = e^(-x)-5*(x-1)^2')
plt.grid(True)


pylab.subplot (2, 2, 3)
plt.plot(x, derivative_of_f(x))
plt.xlabel('x')
plt.ylabel('f\'(x)')
plt.title('f\'(x) = e^(-x)-10*(x-1)')
plt.grid(True)

pylab.subplot (2, 2, 2)
plt.plot(range(0, len(x_arr_1)), x_arr_1)
plt.plot(range(0, len(x_arr_1)), [abs(root_x0[1] - i) for i in x_arr_1])
plt.xlabel('iteration')
plt.ylabel('x')
plt.title('Первое приближение')
plt.grid(True)

pylab.subplot (2, 2, 4)
plt.plot(range(0, len(x_arr_2)), x_arr_2)
plt.plot(range(0, len(x_arr_2)), [abs(root_x0[0] - i) for i in x_arr_2])
plt.xlabel('iteration')
plt.ylabel('x')
plt.title('Второе приближение')
plt.grid(True)

pylab.subplot (2, 3, 5)

difference_array = [abs(x_arr_2[i] - x_arr_1[i]) for i in range(min(len(x_arr_2), len(x_arr_1)))]
plt.plot(range(0, len(difference_array)), difference_array)

plt.plot(range(0, len(difference_array)), [abs(g(x)-x) for x in difference_array])
plt.xlabel('iteration')
plt.ylabel('Значение')
plt.title('Сравнение 5 и 6 пунктов графически')
plt.grid(True)


plt.show()