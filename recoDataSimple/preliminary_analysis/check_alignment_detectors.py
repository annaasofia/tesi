import ROOT
import sys
sys.path.append("..")
import plotting_utils as pu
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import os

ROOT.ROOT.EnableImplicitMT()
ROOT.gStyle.SetOptStat(0)

file_num = 8430
filename = f"../data/recoDataSimple_{file_num}_xtalMerging.root"

print(f"Analyzing {filename} for alignment offset...")

os.makedirs(f"../plots_{file_num}_py", exist_ok=True)


CRYSTAL_BOX = {
    # file_num: (x_min, x_max, y_min, y_max)
    8430: (-1.1044, 0.8965, -5.0250, 7.7750),
    8431: (-1.1044, 0.8965, -5.1000, 7.7000),
    8650: (-0.2744, 1.7256, -5.7450, 7.0550),
    8655: (0.5536, 2.5536, -6.9150, 5.8850),
    8656: (0.5546, 2.5546, -6.8750, 5.9250),
}

x_min, x_max, y_min, y_max = CRYSTAL_BOX[file_num]
 
if x_min == x_max or y_min == y_max:
    print("WARNING: crystal box not set (x_min/x_max/y_min/y_max). "
          "Paste the values printed by the main macro before running this.")
    sys.exit(1)

margin_x = 1.5  # mm
margin_y = 1.5  # mm

df = ROOT.RDataFrame("simpleEvent", filename)
df = df.Filter("SingleTrack == 1", "Single Track")

# Spatial Cut: Select particles that definitely missed the crystal.
air_cut = (
    f"(Tracks.d0_x < {x_min - margin_x} || Tracks.d0_x > {x_max + margin_x} || "
    f"Tracks.d0_y < {y_min - margin_y} || Tracks.d0_y > {y_max + margin_y})"
)
df_air = df.Filter(air_cut, "Air only (missed crystal)")

N_phys = df.Count().GetValue()
N_air = df_air.Count().GetValue()
print(f"Air-only tracks: {N_air} / {N_phys} "
      f"({N_air / N_phys * 100:.2f}% of quality-filtered tracks)")

# Calculate the difference between downstream and upstream projections
# If perfectly aligned, d0Out - d0 should be centered exactly at 0.
df_diff = df_air.Define("diff_x", "Tracks.d0Out_x - Tracks.d0_x") \
                .Define("diff_y", "Tracks.d0Out_y - Tracks.d0_y")

# Create Histograms
h_diff_x = df_diff.Histo1D(("h_diff_x", "Alignment Offset X (Air Particles); x(out) - x(in) [mm]; Counts", 4000, -0.4, 0.4), "diff_x")
h_diff_y = df_diff.Histo1D(("h_diff_y", "Alignment Offset Y (Air Particles); y(out) - y(in) [mm]; Counts", 4000, -0.4, 0.4), "diff_y")

hx = h_diff_x.GetValue()
# Fit X
fit_x = ROOT.TF1("fit_x", "gaus", -0.2, 0.2)
fit_x.SetParameters(h_diff_x.GetMaximum(), h_diff_x.GetMean(), h_diff_x.GetRMS())
hx.Fit(fit_x, "RQ")
offset_x = fit_x.GetParameter(1)
offset_x_err = fit_x.GetParError(1)
sigma_x = fit_x.GetParameter(2)

hy = h_diff_y.GetValue()
# Fit Y
fit_y = ROOT.TF1("fit_y", "gaus", -0.2, 0.2)
fit_y.SetParameters(h_diff_y.GetMaximum(), h_diff_y.GetMean(), h_diff_y.GetRMS())
hy.Fit(fit_y, "RQ")
offset_y = fit_y.GetParameter(1)
offset_y_err = fit_y.GetParError(1)
sigma_y = fit_y.GetParameter(2)

# c1 = ROOT.TCanvas("c1", "Detector Alignment Check", 1200, 600)
# c1.Divide(2, 1)
# c1.cd(1)
# hx.SetFillColor(ROOT.kAzure - 3)
# hx.Draw()
# line0x = ROOT.TLine(0, 0, 0, h_diff_x.GetMaximum())
# line0x.SetLineColor(ROOT.kBlack); line0x.SetLineStyle(2); line0x.SetLineWidth(2); line0x.Draw("SAME")
# line1x = ROOT.TLine(offset_x, 0, offset_x, h_diff_x.GetMaximum())
# line1x.SetLineColor(ROOT.kRed); line1x.SetLineStyle(2); line1x.SetLineWidth(2); line1x.Draw("SAME")
# c1.cd(2)
# hy.SetFillColor(ROOT.kOrange + 1)
# hy.Draw()
# line0y = ROOT.TLine(0, 0, 0, h_diff_y.GetMaximum())
# line0y.SetLineColor(ROOT.kBlack); line0y.SetLineStyle(2); line0y.SetLineWidth(2); line0y.Draw("SAME")
# line1y = ROOT.TLine(offset_y, 0, offset_y, h_diff_y.GetMaximum())
# line1y.SetLineColor(ROOT.kRed); line1y.SetLineStyle(2); line1y.SetLineWidth(2); line1y.Draw("SAME")
# c1.Update()

print("=" * 50)
print("GAUSSIAN FIT (offset between IN and OUT detector systems):")
print(f"\tOffset x = {offset_x:.4f} +/- {offset_x_err:.4f} mm  (sigma = {sigma_x:.4f} mm)")
print(f"\t\tChi^2/NDF = {(fit_x.GetChisquare() / fit_x.GetNDF()):.3f}")
print(f"\t\tSkewness = {h_diff_x.GetSkewness():.4f}")
print(f"\tOffset y = {offset_y:.4f} +/- {offset_y_err:.4f} mm  (sigma = {sigma_y:.4f} mm)")
print(f"\t\tChi^2/NDF = {(fit_y.GetChisquare() / fit_y.GetNDF()):.3f}")
print(f"\t\tSkewness = {h_diff_y.GetSkewness():.4f}")


h_prof_x = df_air.Profile1D(
    ("h_prof_x", "x(out) vs x(in) (air only); x(in) [mm]; x(out) [mm]",
     100, -5, 15), "Tracks.d0_x", "Tracks.d0Out_x")
h_prof_y = df_air.Profile1D(
    ("h_prof_y", "y(out) vs y(in) (air only); y(in) [mm]; y(out) [mm]",
     100, -15, 15), "Tracks.d0_y", "Tracks.d0Out_y")
 
h_prof_x_val = h_prof_x.GetValue()
h_prof_y_val = h_prof_y.GetValue()
 
lin_x = ROOT.TF1("lin_x", "pol1", -5, 15)
h_prof_x_val.Fit(lin_x, "RQ0")
slope_x = lin_x.GetParameter(1)
intercept_x = lin_x.GetParameter(0)
 
lin_y = ROOT.TF1("lin_y", "pol1", -15, 15)
h_prof_y_val.Fit(lin_y, "RQ0")
slope_y = lin_y.GetParameter(1)
intercept_y = lin_y.GetParameter(0)
 
print("=" * 50)
print("LINEAR FIT CROSS-CHECK (d0Out = slope * d0 + intercept):")
print(f"\tx: slope = {slope_x:.4f} intercept = {intercept_x:.4f} mm")
print(f"\ty: slope = {slope_y:.4f} intercept = {intercept_y:.4f} mm")
# print("\t(slope far from 1 would suggest rotation/scale, not just a shift)")


# c2 = ROOT.TCanvas("c2", "Detector alignment - profile fits", 1400, 700)
# c2.Divide(2, 1)
# c2.cd(1)
# lin_x.SetLineColor(ROOT.kRed)
# h_prof_x_val.Draw()
# lin_x.Draw("SAME")
# h_prof_x_val.Draw("SAME")
# c2.cd(2)
# lin_y.SetLineColor(ROOT.kRed)
# h_prof_y_val.Draw()
# lin_y.Draw("SAME")
# h_prof_y_val.Draw("SAME")
# c2.Update()

fig1, axs1 = plt.subplots(1, 2, figsize=(14, 6))

# X Offset
pu.plot_histo1d(h_diff_x, fit_func=fit_x, ax=axs1[0], style="fill", color="cornflowerblue", alpha=1.0,
                xlabel="x(out) - x(in) [mm]", ylabel="Counts", data_label=None, fit_label=None, title="Alignment Offset x")
axs1[0].axvline(0, color='black', ls='--', lw=2, label='Zero')
axs1[0].axvline(offset_x, color='red', ls='--', lw=2, label='Fit Offset')
# axs1[0].get_legend().remove()
axs1[0].legend(loc='upper right', frameon=False)
# pu.add_info_box(axs1[0], [f"Offset x = {offset_x:.4f} mm", f"Sigma = {sigma_x:.4f} mm"], loc="upper left")

# Y Offset
pu.plot_histo1d(h_diff_y, fit_func=fit_y, ax=axs1[1], style="fill", color="orange", alpha=1.0,
                xlabel="y(out) - y(in) [mm]", ylabel="Counts", data_label=None, fit_label=None, title="Alignment Offset y")
axs1[1].axvline(0, color='black', ls='--', lw=2, label='Zero')
axs1[1].axvline(offset_y, color='red', ls='--', lw=2, label='Fit Offset')
# axs1[1].get_legend().remove()
axs1[1].legend(loc='upper right', frameon=False)
# pu.add_info_box(axs1[1], [f"Offset y = {offset_y:.4f} mm", f"Sigma = {sigma_y:.4f} mm"], loc="upper left")

fig1.tight_layout()
fig1.savefig(f"../plots_{file_num}_py/alignment_check.pdf")
fig1.savefig(f"../plots_{file_num}_py/alignment_check.png")
plt.close(fig1)



fig2, axs2 = plt.subplots(1, 2, figsize=(14, 6))

# Profile X (usando style="errorbar" è perfetto per i TProfile)
pu.plot_histo1d(h_prof_x_val, fit_func=lin_x, ax=axs2[0], style="errorbar", color="black",
                xlabel="x(in) [mm]", ylabel="x(out) [mm]", title="x(out) vs x(in) (air only)")
pu.add_info_box(axs2[0], [f"Slope = {slope_x:.4f}", f"Intercept = {intercept_x:.4f} mm"], loc="upper left", fontsize=14)

# Profile Y
pu.plot_histo1d(h_prof_y_val, fit_func=lin_y, ax=axs2[1], style="errorbar", color="black",
                xlabel="y(in) [mm]", ylabel="y(out) [mm]", title="y(out) vs y(in) (air only)")
pu.add_info_box(axs2[1], [f"Slope = {slope_y:.4f}", f"Intercept = {intercept_y:.4f} mm"], loc="upper left", fontsize=14)

fig2.tight_layout()
fig2.savefig(f"../plots_{file_num}_py/linear_check.pdf")
fig2.savefig(f"../plots_{file_num}_py/linear_check.png")
plt.close(fig2)

# c1.SaveAs(f"../plots_{file_num}/alignment_check.png")
# c1.SaveAs(f"../plots_{file_num}/alignment_check.pdf")
# c2.SaveAs(f"../plots_{file_num}/linear_check.png")
# c2.SaveAs(f"../plots_{file_num}/linear_check.pdf")

 

