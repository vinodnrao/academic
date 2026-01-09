# Python Code for Figure 2 in `Protocols for counterfactual and twin-field quantum digital signature'

import matplotlib.pyplot as plt
import numpy as np

f = np.linspace(0.001,1,100)
#n = np.linspace(1,100,100)

x = -((f * np.log2(f+0.0001)) + ((1-f) * np.log2(1-f+0.0001)))

e1 = (4/(1 - f+0.0001))
e2 = (4/((1 + f) * (1 - x)))

plt.figure(dpi=1200)
plt.plot (f,e1,label="Singleton",linestyle="solid")
plt.plot (f,e2,label="Hamming",linestyle="dashed")

plt.xlabel(r'ratio r of injected bits', fontsize=12, fontweight='bold')
plt.ylabel(r'total number N of bits', fontsize=12, fontweight='bold')
plt.grid()

plt.xlim([0, 0.4])
plt.ylim([0, 100])
plt.legend()
#plt.show()

plt.savefig("error_correction.png", dpi = 2000)
plt.savefig("Figure2.pdf", bbox_inches='tight')





