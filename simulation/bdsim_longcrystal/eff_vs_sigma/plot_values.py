import sys
import numpy as np
import matplotlib.pyplot as plt

s_x, s_y, eff_gaus, eff_gaus_err, eff_3, eff_3_err = np.loadtxt("values_bdsim.csv", delimiter=",", skiprows=1, unpack=True)
s_x2, s_y2, eff_gaus2, eff_gaus_err2, eff_32, eff_3_err2 = np.loadtxt("values_xcoll.csv", delimiter=",", skiprows=1, unpack=True)
o_x, o_y = np.argsort(s_x), np.argsort(s_y)
o_x2, o_y2 = np.argsort(s_x2), np.argsort(s_y2)

plt.figure(1, figsize=(8,6))
# plt.errorbar(s_x[o_x], eff_gaus[o_x], yerr=eff_gaus_err[o_x], fmt="o--", color="orange", capsize=5, label="xcoll gaus")
plt.errorbar(s_x[o_x], eff_3[o_x], yerr=eff_3_err[o_x], fmt="o--", color="#68238B", label="Xcoll")
plt.errorbar(s_x2[o_x2], eff_32[o_x2], yerr=eff_3_err2[o_x2], fmt="o--", color="#982916", label="BDSIM")
# plt.title(r"TCCP: $\varepsilon$ vs $\sigma_{detector}$", fontsize=16)
plt.title(r"Long crystal: efficiency vs detector uncertainty", fontsize=16)
plt.xlabel(r"$\sigma_{detector}$ [$\mu$rad]", fontsize=14)
plt.ylabel("Channeling efficiency [%]", fontsize=14)
plt.xticks(fontsize=11)
plt.yticks(fontsize=11)
plt.grid(alpha=0.3)
plt.legend(fontsize=12)
plt.savefig("eff_vs_sigma.png", dpi=200)

# x, y, _, _, z, _ = np.loadtxt("values_diff_xcoll.csv", delimiter=",", unpack=True)

# plt.figure(2, figsize=(7,6))
# eff_map = plt.scatter(x, y, c=z, cmap='plasma', marker='s', s=350)
# plt.colorbar(eff_map, label="Channeling efficiency [%]")
# plt.title(r"Efficiency map $\sigma_{in}$ & $\sigma_{out}$", fontsize=16)
# plt.xlabel(r"$\sigma_{in}$ [$\mu$rad]", fontsize=16)
# plt.ylabel(r"$\sigma_{out}$ [$\mu$rad]", fontsize=16)
# plt.xticks(range(16), fontsize=12)
# plt.yticks(range(16), fontsize=12)

# z_mid = (z.min() + z.max()) / 2
# for xi, yi, zi in zip(x, y, z):
#     if xi == yi:
#         color = "white" if zi < z_mid else "black"
#         plt.text(
#             xi, yi, f"{zi:.1f}",
#             color=color,
#             ha="center", va="center",
#             fontsize=8, fontweight="bold"
#         )
# plt.savefig("eff_vs_sigma_diff.png", dpi=200)

# plt.show()




