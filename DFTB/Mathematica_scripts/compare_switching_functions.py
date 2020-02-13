import numpy as np
from matplotlib import pyplot as plt

approx_0p1 = np.loadtxt("approx_0p1.dat")
approx_0p2 = np.loadtxt("approx_0p2.dat")
approx_0p333 = np.loadtxt("approx_0p333.dat")

exact_0p1 = np.loadtxt("exact_0p1.dat")
exact_0p2 = np.loadtxt("exact_0p2.dat")
exact_0p333 = np.loadtxt("exact_0p333.dat")

plt.xlabel("$R_{AB}$ / Bohr", fontsize=17)
plt.ylabel("$\gamma^{lr}_{AB}$ / Hartree", fontsize=17)

lw=3
plt.plot(exact_0p1[:,0], exact_0p1[:,1], color="blue", lw=lw, label="$\\omega = \\frac{1}{10}$ (exact)")
plt.plot(approx_0p1[:,0], approx_0p1[:,1], color="blue", lw=lw, ls="-.", label="$\\omega = \\frac{1}{10}$ (approx.)")
plt.plot(exact_0p2[:,0], exact_0p2[:,1], color="red", lw=lw, label="$\\omega = \\frac{1}{5}$ (exact)")
plt.plot(approx_0p2[:,0], approx_0p2[:,1], color="red", ls="-.", lw=lw, label="$\\omega = \\frac{1}{5}$ (approx.)")
plt.plot(exact_0p333[:,0], exact_0p333[:,1], color="green", lw=lw, label="$\\omega = \\frac{1}{3}$ (exact)")
plt.plot(approx_0p333[:,0], approx_0p333[:,1], color="green", ls="-.", lw=lw, label="$\\omega = \\frac{1}{3}$ (approx.)")

plt.legend()
plt.savefig("comparison_gamma_lr_exact_approx.png")
plt.show()

