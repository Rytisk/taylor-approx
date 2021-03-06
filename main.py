import sympy as sy
import numpy as np
from sympy.functions import sin,cos, exp, log
import matplotlib.pyplot as plt
import math

plt.style.use("ggplot")

# Define the variable and the function to approximate
x = sy.Symbol('x')
f = exp(x) * log(1+x)

# Factorial function
def factorial(n):
	if n <= 0:
		return 1
	else:
		return n*factorial(n-1)

# Taylor approximation at x0 of the function 'function'
def taylor(function,x0,n):
	i = 0
	p = 0
	while i <= n:
		p = p + (function.diff(x,i).subs(x,x0))/(factorial(i))*(x-x0)**i
		i += 1
	return p

# Plot results
def plot():
	x_lims = [-1/2,2]
	x1 = np.linspace(x_lims[0],x_lims[1],800)
	y1 = []
	# Approximate up until 10 starting from 1 and using steps of 2
	for j in range(1,6,1):
		func = taylor(f,0,j)
		print('Taylor expansion at n='+str(j),func)
		for k in x1:
			y1.append(func.subs(x,k))
		plt.plot(x1,y1,label='order '+str(j))
		y1 = []
	# Plot the function to approximate
	fs = list(map(lambda x0: math.exp(x0) * math.log(1+x0), x1))
	plt.plot(x1, fs,label='f of x')
	plt.xlim(x_lims)
	plt.ylim([-1/2,2])
	plt.xlabel('x')
	plt.ylabel('y')
	plt.legend()
	plt.grid(True)
	plt.title('Taylor series approximation')
	plt.show()

plot()