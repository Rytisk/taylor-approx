import sympy as sy
import numpy as np
from sympy.functions import sin, cos, ln, exp
import matplotlib.pyplot as plt
import sys
import math

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

def lagrange(function, x0, n, x = sy.Symbol('x'), E = sy.Symbol('E')):
	i = n + 1
	p = (function.diff(x, i).subs(x, E))/(factorial(i))*(x - x0)**i
	return p

def run(f, n):
	func = taylor(f, 0, n)
	print('Taylor expansion at n=' + str(n), func)

	lagr = lagrange(f, 0, n)
	lagr_lambda = sy.lambdify((x,E), lagr, "numpy")
	print('Lagrange remainder', lagr)
	print('\n')
	print("Accuracy: {0:.10f}".format(lagr_lambda(1, 1)))

def run_till_accuracy(f, acc):
	l_acc = 0
	n = 0
	while True:
		lagr = lagrange(f, 0, n)
		lagr_lambda = sy.lambdify((x,E), lagr, "numpy")
		l_acc = abs(lagr_lambda(1, 1))
		if l_acc <= acc:
			break
		n = n + 1
	print("To reach accuracy of {}, {} terms are needed".format(acc, n))
	print("Accuracy: {0:.10f}".format(l_acc))

def usage():
	print("-n [order]")
	print("-a [accuracy]")
	exit(0)

if __name__ == "__main__":
	x = sy.Symbol('x')
	E = sy.Symbol('E')
	f = exp(x) * ln(1 + x)

	if(len(sys.argv) == 3):
		arg = sys.argv[1]
		if arg == "-n":
			n = int(sys.argv[2])
			run(f, n)
		elif arg == "-a":
			acc = float(sys.argv[2])
			run_till_accuracy(f, acc)
		else:
			usage()
	else:
		usage()	