import numpy as np

def fun(x, y):

    return np.vstack((y[1], -np.exp(y[0])))

#Implement evaluation of the boundary condition residuals:

def bc(ya, yb):

    return np.array([ya[0], yb[0]])

#Define the initial mesh with 5 nodes:

x = np.linspace(0, 1, 5)

#This problem is known to have two solutions. To obtain both of them, we use two different initial guesses for y. We denote them by subscripts a and b.

y_a = np.zeros((2, x.size))

y_b = np.zeros((2, x.size))
y_b[0] = 3

#Now we are ready to run the solver.

from scipy.integrate import solve_bvp

res_a = solve_bvp(fun, bc, x, y_a)
res_b = solve_bvp(fun, bc, x, y_b)

#Let’s plot the two found solutions. We take an advantage of having the solution in a spline form to produce a smooth plot.
x_plot = np.linspace(0, 1, 100)
y_plot_a = res_a.sol(x_plot)[0]
y_plot_b = res_b.sol(x_plot)[0]

import matplotlib.pyplot as plt

plt.plot(x_plot, y_plot_a, label='y_a')
plt.plot(x_plot, y_plot_b, label='y_b')
plt.legend()
plt.xlabel("x")
plt.ylabel("y")
plt.show()