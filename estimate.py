import sys
import sympy as sy
import numpy as np
from sympy.functions import ln, exp
import matplotlib.pyplot as plt

plt.style.use("ggplot")

def factorial(n):
	if n <= 0:
		return 1
	else:
		return n * factorial(n - 1)

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

def run(f, n, interval):
	func = taylor(f, 0, n)
	print('Taylor expansion at n=' + str(n), func)
	lagr = lagrange(f, 0, n)
	lagr_lambda = sy.lambdify((x,E), lagr, "numpy")
	print('Lagrange remainder', lagr)
	print('\n')
	low_x = abs(lagr_lambda(interval[1], 0))
	high_x = abs(lagr_lambda(interval[1], interval[1]))
	acc = low_x if low_x > high_x else high_x
	print("Accuracy: {0:.20f}".format(acc))


def run_till_accuracy(f, acc, interval):
	l_acc = 0
	n = 0
	while True:
		lagr = lagrange(f, 0, n)
		lagr_lambda = sy.lambdify((x,E), lagr, "numpy")
		l_acc = abs(lagr_lambda(interval[1], interval[1]))
		if l_acc <= acc:
			break
		n = n + 1
	print("To reach accuracy of {}, {} terms are needed".format(acc, n))
	print("Accuracy: {0:.10f}".format(l_acc))

def usage():
	print("-n [order] [from_x] [to_x]")
	print("-a [accuracy] [from_x] [to_x]")
	exit(0)

if __name__ == "__main__":
	x = sy.Symbol('x')
	E = sy.Symbol('E')
	f = exp(x) * ln(1 + x)

	if(len(sys.argv) == 5):
		arg = sys.argv[1]
		if arg == "-n":
			n = int(sys.argv[2])
			from_x = float(sys.argv[3])
			to_x = float(sys.argv[4])
			run(f, n, [from_x, to_x])
		elif arg == "-a":
			acc = float(sys.argv[2])
			from_x = float(sys.argv[3])
			to_x = float(sys.argv[4])
			run_till_accuracy(f, acc, [from_x, to_x])
		else:
			usage()
	else:
		usage()
