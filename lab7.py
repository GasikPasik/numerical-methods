
import matplotlib.pyplot as plt
import pylab

a = 0
b = 0

def f(t, u):
    return u**2 / (t+1)

def solve(n):

    h = (b-a)/n
    t = [a + k*h for k in range(0, n)]

    return h, t

def main():
    global a, b
    a = 0.75
    b = 3

    n = 5000

    x = [1]
    h, t = solve(n)

    for k in range(1, n):
        x.append(x[k-1] + h * f(t[k-1], x[k-1]))


    fig, axes = plt.subplots(ncols=2, nrows=2, sharex=False, sharey=False)
    fig.subplots_adjust(hspace=0.5, wspace=0.2)

    pylab.subplot (2, 2, 1)
    plt.plot(t, x)
    plt.xlabel("t[k]")
    plt.ylabel("x [k]")
    plt.grid(True)

    n *= 2
    y = [1]
    h, t = solve(n)

    for k in range(1, n):
        y.append(y[k-1] + h * f(t[k-1], y[k-1]))

    z1 = [abs(x[m] - y[2*m])  for m in range(0, n//2)]
    print(f"max(z1): {max(z1)}")


    print("Разностная схема Эйлера-Коши")

    n = 200
    h, t = solve(n)
    x1 = [1]
    for k in range(1, n):
        k -= 1
        x1.append( x1[k] + h/2 * (f(t[k], x1[k]) + f(t[k+1], x1[k] + h*f(t[k], x1[k]))) )


    pylab.subplot (2, 2, 2)
    plt.plot(t, x1)
    plt.xlabel("t[k]")
    plt.ylabel("x [k]")
    plt.grid(True)
    plt.legend()

    

    plt.show()

    

if __name__ == "__main__":
    main()