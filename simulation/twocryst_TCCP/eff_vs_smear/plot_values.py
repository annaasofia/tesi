import numpy as np
import matplotlib.pyplot as plt

s_x, s_y, eff_gaus, eff_gaus_err, eff_3, eff_3_err = np.loadtxt("values_xcoll.csv", delimiter=",", skiprows=1, unpack=True)
o_x, o_y = np.argsort(s_x), np.argsort(s_y)
idx = np.where(s_x == 4.5)[0]

plt.figure(1, figsize=(8,6))
# plt.errorbar(s_x[o_x], eff_gaus[o_x], yerr=eff_gaus_err[o_x], fmt="o--", color="orange", capsize=5, label="xcoll gaus")
plt.errorbar(s_x[o_x], eff_3[o_x], yerr=eff_3_err[o_x], fmt="o--", color="royalblue", label="xcoll 3σ")
plt.errorbar(s_x[idx], eff_3[idx], yerr=eff_3_err[idx], fmt="o", color="orange", zorder=5, label="measured TCCP")
# plt.title(r"TCCP: $\varepsilon$ vs $\sigma_{detector}$", fontsize=16)
plt.title(r"TCCP: efficiency vs detector uncertainty", fontsize=16)
plt.xlabel(r"$\sigma_{det}$ [$\mu$rad]", fontsize=14)
plt.ylabel("Channeling efficiency [%]", fontsize=14)
plt.xticks(fontsize=11)
plt.yticks(fontsize=11)
plt.grid(alpha=0.3)
plt.legend(fontsize=12)
plt.savefig("eff_vs_sigma.png", dpi=200)

x, y, _, _, z, _ = np.loadtxt("values_diff_xcoll.csv", delimiter=",", unpack=True)

plt.figure(2, figsize=(7,6))
eff_map = plt.scatter(x, y, c=z, cmap='plasma', marker='s', s=350)
plt.colorbar(eff_map, label="Channeling efficiency [%]")
plt.title(r"Efficiency map $\sigma_{in}$ - $\sigma_{out}$")
plt.xlabel(r"$\sigma_{in}$ [$\mu$rad]", fontsize=14)
plt.ylabel(r"$\sigma_{out}$ [$\mu$rad]", fontsize=14)
plt.xticks(range(16))
plt.yticks(range(16))
plt.savefig("eff_vs_sigma_diff.png", dpi=200)

# plt.show()




