import numpy as np
from matplotlib import rc
import matplotlib.pyplot as plt

# Figure dpi
dpi = 72
save_name = ["unstableneg.png", "stableneg.png", "unstablepos.png", "stablepos.png"]

def plot_cobweb(f, b, k_max, gamma, k_d, n, x0, interval, nmax=50, n_arrows=5, name_fig=None):

    x = np.linspace(0, 1, 500)
    fig = plt.figure(figsize=(600/dpi, 450/dpi), dpi=dpi)
    ax = fig.add_subplot(111)

    # Plot y = f(x) and y = x
    ax.plot(x, f(x, b, k_max, gamma, k_d, n), c='#444444', lw=2)
    ax.plot(x, x, c='#444444', lw=2, linestyle='dashed')

    # Iterate x = f(x) for nmax steps, starting at (x0, 0).
    px, py = np.empty((2, nmax+1, 1))
    px[0], py[0] = x0, 0
    for i in range(1, nmax, 2):
        px[i] = px[i-1]
        py[i] = f(px[i], b, k_max, gamma, k_d, n)
        px[i+1] = py[i]
        py[i+1] = py[i]

    # Plot the path traced out by the iteration.
    ax.plot(px, py, c='b')
    ax.set_ylim((0,1))
    ax.set_xlim((0,1))

    # Add arrows
    arrow_interval = 1  # Determines how frequently arrows are drawn
    for n in range(1, n_arrows, arrow_interval):
        # Calculate the start and end points of the arrow
        start_x = px[n - 1][0]
        start_y = py[n - 1][0]
        if n < n_arrows - 1:  # To avoid index out of range
            end_x = (px[n][0] - start_x)/2
            end_y = (py[n][0] - start_y)/2

            # Check for zero-length arrow to avoid plotting issues
            if end_x != 0 or end_y != 0:
                # Draw the arrow
                ax.arrow(start_x, start_y, end_x, end_y, head_width=0.02, head_length=0.03, fc='red', ec='red')

    # Annotate and tidy the plot.
    ax.minorticks_on()
    ax.grid(which='minor', alpha=0.5)
    ax.grid(which='major', alpha=0.5)
    ax.set_aspect('equal')
    ax.set_xlabel('$q_{i-1}$')
    ax.set_ylabel(f.latex_label)
    ax.set_title(interval)

    plt.savefig(name_fig, dpi=dpi)
    plt.show()


class AnnotatedFunction:
    # A class representing a mathematical function.

    def __init__(self, func, latex_label):
        self.func = func
        self.latex_label = latex_label

    def __call__(self, *args, **kwargs):
        return self.func(*args, **kwargs)


diff_eq = AnnotatedFunction(lambda x,b,k_max,gamma,k_d,n: b*k_max/(gamma*(1 + ((x/k_d)**n))), r'$q_i$')
diff_eq_rep = AnnotatedFunction(lambda x,b,k_max,gamma,k_d,n: b*k_max/(gamma*(1 + ((k_d/x)**n))), r'$q_i$')
h_diff = AnnotatedFunction(lambda x,b,k_max,gamma,k_d,n: -(n/(k_max*k_d))*((x/k_d)**(n-1))*((k_max/(1 + ((x/k_d)**n)))**2), r'$q_i$')
h_diff_rep = AnnotatedFunction(lambda x,b,k_max,gamma,k_d,n: (n/(k_max*k_d))*((k_d/x)**(n+1))*((k_max/(1 + ((k_d/x)**n)))**2), r'$q_i$')

condition = ["$bh'/\Gamma<-1$","$-1<bh'/\Gamma<0$","$1<bh'/\Gamma$","$0<bh'/\Gamma<1$"]
b = [0.5, 0.8, 0.8, 0.5]
k_max = [1.8, 1.9, 2, 2.5]
gamma = [1, 1.8, 0.7, 0.8]
k_d = [0.2, 0.7, 0.7, 0.8]
n = [2, 2, 2, 1]
x0 = 0.4
x = np.linspace(0, 1, 500)

from scipy.optimize import fsolve


def equations(var, *params):
    x = var
    bb, kk_max, ggamma, kk_d, nn, j = params
    if j <= 1:
        return diff_eq(x, bb, kk_max, ggamma, kk_d, nn) - x
    else:
        return diff_eq_rep(x, bb, kk_max, ggamma, kk_d, nn) - x

initial_guess = 0.5
solution = []
data = [[0]*3 for _ in range(4)]

for i in range(4):
    data[i][0] = i+1
    data[i][1] = fsolve(equations, initial_guess, args=(b[i], k_max[i], gamma[i], k_d[i], n[i], i))
    if i <= 1:
        data[i][2] = b[i]*h_diff(data[i][1], b[i], k_max[i], gamma[i], k_d[i], n[i])/gamma[i]
        plot_cobweb(diff_eq, b[i], k_max[i], gamma[i], k_d[i], n[i], x0, interval=condition[i], n_arrows=12, name_fig=save_name[i])
    else:
        data[i][2] = b[i]*h_diff_rep(data[i][1], b[i], k_max[i], gamma[i], k_d[i], n[i])/gamma[i]
        plot_cobweb(diff_eq_rep, b[i], k_max[i], gamma[i], k_d[i], n[i], x0, interval=condition[i], n_arrows=12, name_fig=save_name[i])


from tabulate import tabulate
print (tabulate(data, headers=["Scenario", "Fixed point", "Stability condition"]))



