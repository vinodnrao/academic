# Python Code for Figure 3 in `Protocols for counterfactual and twin-field quantum digital signature'

import numpy as np
import matplotlib.pyplot as plt

# Define the range of w
w = np.linspace(0.01, 1, 100)

# Define the function for error rate e(w, r)
def error_rate(w, r):
    return 0.5 * (((w * (1 + r)) + r) / (1 + r + 2 * w * (1 + r)))

# Define the Shannon entropy for binary variable
def entropy(e):
    return - (e * np.log2(e) + (1 - e) * np.log2(1 - e))

# Define mutual information difference I(AB) - I(AE)
def mutual_info(e, w):
    x = entropy(e)
    return 1 - x - (w / 2)

# List of r values to plot
r_values = [0.1, 0.01]
styles = ['solid', 'dashed']
labels = [r"$r = 0.1$", r"$r = 0.01$"]

plt.figure()

for r, style, label in zip(r_values, styles, labels):
    e = error_rate(w, r)
    I = mutual_info(e, w)
    plt.plot(e, I, label=label, linestyle=style)

plt.xlabel(r'Error rate', fontsize=12, fontweight='bold')
plt.ylabel(r'$I(AB) - I(AE)$', fontsize=12, fontweight='bold')
plt.grid(True)
plt.xlim([0, 0.25])
plt.ylim([0, 0.8])
plt.legend()

plt.savefig("error_rate.png", dpi=2000)
plt.savefig("Figure3.pdf", bbox_inches='tight')
plt.show()

