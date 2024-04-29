import matplotlib
matplotlib.use('tkagg')
import matplotlib.pyplot as plt
from scipy.integrate import odeint
import numpy as np
from tabulate import tabulate

y0 = [0, 0, 0, 0, 0, 0]

t = np.linspace(0,300, num=100)

gamma = [0.08, 0.05, 0.05]
Gamma = [0.008, 0.05, 0.05]
k1 = 1
k_max = 2
k_d = 2
K = [0.01, 0.05, 0.05]
n = 2


def h(x, kmax, kd, ni):
    return kmax / (1 + ((x/kd)**ni))


def sim(variables, t, gamma, Gamma, k1, k_max, k_d, K, n):

    m1 = variables[0]
    p1 = variables[1]
    m2 = variables[2]
    p2 = variables[3]
    m3 = variables[4]
    p3 = variables[5]

    dm1dt = k1 - gamma[0] * m1
    dp1dt = K[0]*m1 - Gamma[0]*p1
    dm2dt = h(p1, k_max, k_d, n) - gamma[1]*m2
    dp2dt = K[1]*m2 - Gamma[1]*p2
    dm3dt = h(p2, k_max, k_d, n) - gamma[2] * m3
    dp3dt = K[2] * m3 - Gamma[2] * p3

    return [dm1dt, dp1dt, dm2dt, dp2dt, dm3dt, dp3dt]


y = odeint(sim, y0, t, args=(gamma, Gamma, k1, k_max, k_d, K, n,))

dpi = 90
fig = plt.figure(figsize=(500/dpi, 450/dpi), dpi=dpi)
ax = fig.add_subplot(111)
#f, ax = plt.subplots(1)

m, p = np.empty((2, 3, 1))
m[0], p[0] = k1/gamma[0], K[0]*k1/(gamma[0]*Gamma[0])
plt.axhline(y = m[0], color="lightsteelblue",linestyle="--", label = "stable state")
plt.axhline(y = p[0], color="lightsteelblue",linestyle="--")
for i in range(1,3):
    m[i] = h(p[i-1],k_max,k_d,n)/ gamma[i]
    plt.axhline(y = m[i], color="lightsteelblue",linestyle="--")
    p[i] = K[i]/Gamma[i] * m[i]
    plt.axhline(y = p[i], color="lightsteelblue",linestyle="--")

labels = ["m1", "p1", "m2", "p2", "m3", "p3"]
col = ["b", "b", "g", "g", "r", "r"]
line_i = []

for i in range(6):
    if i%2 == 1:
        line_i.append(ax.plot(t, y[:, i], color=col[i], label=labels[i])[0])
    else:
        line_i.append(ax.plot(t, y[:, i], color=col[i], label=labels[i], linestyle="--")[0])



steady_state = []
for j in range(3):
    steady_state.append([m[j], p[j]])

# create header
head = ["mRNA", "Protein"]

# display table
print(tabulate(steady_state, headers=head, tablefmt="grid"))

ax.set_ylabel('Koncentrácia')
ax.set_xlabel('Čas')
plt.title("Trojstupňová genetická regulačná kaskáda")
ax.legend(handles=line_i, loc="upper right")

plt.grid()
plt.show()
