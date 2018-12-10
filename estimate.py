import sympy as sy
import numpy as np
from sympy.functions import sin, cos, ln, exp
import matplotlib.pyplot as plt
plt.style.use("ggplot")

# Factorial function
def factorial(n):
    if n <= 0:
        return 1
    else:
        return n * factorial(n - 1)

# Taylor approximation at x0 of the function 'function'
def taylor(function, x0, n, x = sy.Symbol('x')):
    i = 0
    p = 0
    while i <= n:
        p = p + (function.diff(x, i).subs(x, x0))/(factorial(i))*(x - x0)**i
        i += 1
    return p

def lagrange(function, x0, n, x = sy.Symbol('x')):
    i = n + 1
    p = (function.diff(x, i).subs(x, x0))/(factorial(i))*(x - x0)**i
    return p

x = sy.Symbol('x')
f = exp(x) * ln(1 + x)

func = taylor(f, 0, 5)
taylor_lambda = sy.lambdify(x, func, "numpy")
print('Taylor expansion at n=' + str(5), func)

lagr = lagrange(f, 0, 5)
print('Lagrange remainder at n=' + str(5), lagr)

#plt.plot(x1, taylor_lambda(x1), label = 'Order '+ str(j))