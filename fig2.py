import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import degree
from RattleBack import RattleBack

abch = (20,3,2,1)
ABCD = (2,16,17,-0.2)
sigma = 1
init = [0.5*degree, 0.5*degree, 0,0,0,5]
t = np.linspace(0,10,600)

y = RattleBack(abch, ABCD, init, t, sigma)
alpha,beta,gamma = y.T[:3]
delta = np.arccos(np.cos(alpha)*np.cos(beta))

plt.figure(figsize=(5,5.5))

plt.subplot(211)
plt.plot(t, gamma/degree)
plt.ylabel(r'$\gamma$  / deg')

plt.subplot(212)
plt.plot(t, delta/degree)
plt.ylabel(r'$\delta$  / deg')
plt.xlabel(r'$t$ = time  / sec')
plt.tight_layout()
plt.savefig('fig2.eps')
plt.show()
