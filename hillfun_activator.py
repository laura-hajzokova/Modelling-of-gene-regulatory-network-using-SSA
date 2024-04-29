import matplotlib
matplotlib.use('tkagg')
import matplotlib.pyplot as plt
from scipy.integrate import odeint
import numpy as np

n1 = 1
n2 = 2
n3 = 3
n4 = 4
n_infty = 1000

k_max = 1
K = 0.5
x = np.linspace(0, 1, num=100)+0.001
params1 = [k_max, K, n1]
params2 = [k_max, K, n2]
params3 = [k_max, K, n3]
params4 = [k_max, K, n4]
params5 = [k_max, K, n_infty]

def hill(x, params):

    k_max = params[0]
    K = params[1]
    n = params[2]

    X = x

    k = k_max / (1 + (K / X) ** n)

    return(k)


y1 = hill(x,params1)
y2 = hill(x,params2)
y3 = hill(x,params3)
y4 = hill(x,params4)
y5 = hill(x,params5)

f, ax = plt.subplots(1)

line1, = ax.plot(x, y1, color="k", label="n=1")
line2, = ax.plot(x, y2, color="b", label="n=2")
line3, = ax.plot(x, y3, color="r", label="n=3")
line4, = ax.plot(x, y4, color="g", label="n=4")
line5, = ax.plot(x, y5, color="m", label="n->∞")

ax.set_ylabel('f(X)')
ax.set_xlabel('X')
plt.title("Hillova funkcia pre aktivator")
ax.legend(handles=[line1, line2, line3, line4, line5], loc="upper left")

plt.show()
