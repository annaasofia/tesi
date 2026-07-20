import ROOT
import sys

ROOT.ROOT.EnableImplicitMT()
ROOT.gStyle.SetOptStat(0)

file_num = 8430
filename = f"recoDataSimple_{file_num}_xtalMerging.root"

print(f"Analyzing {filename} for alignment offset...")


CRYSTAL_BOX = {
    # file_num: (x_min, x_max, y_min, y_max)
    8430: (-1.1000, 0.9000, -2.8288, 5.1712),
    8431: (-1.1000, 0.9000, -2.8090, 5.1910),
    8650: (-0.2700, 1.7300, -3.2639, 4.7361),
    8655: (0.5700, 2.5700, -3.4964, 4.5036),
    8656: (0.5700, 2.5700, -3.5103, 4.4897),
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
h_diff_x = df_diff.Histo1D(("h_diff_x", "Alignment Offset X (Air Particles); x_{out} - x_{in} [mm]; Counts", 4000, -0.4, 0.4), "diff_x")
h_diff_y = df_diff.Histo1D(("h_diff_y", "Alignment Offset Y (Air Particles); y_{out} - y_{in} [mm]; Counts", 4000, -0.4, 0.4), "diff_y")

# Fit with Gaussian to extract the exact offset
c1 = ROOT.TCanvas("c1", "Detector Alignment Check", 1200, 600)
c1.Divide(2, 1)

# X Offset Fit
c1.cd(1)
hx = h_diff_x.GetValue()
hx.SetFillColor(ROOT.kAzure - 3)
hx.Draw()

# Fit X
fit_x = ROOT.TF1("fit_x", "gaus", -0.2, 0.2)
fit_x.SetParameters(h_diff_x.GetMaximum(), h_diff_x.GetMean(), h_diff_x.GetRMS())
hx.Fit(fit_x, "RQ")
offset_x = fit_x.GetParameter(1)
offset_x_err = fit_x.GetParError(1)
sigma_x = fit_x.GetParameter(2)

line0x = ROOT.TLine(0, 0, 0, h_diff_x.GetMaximum())
line0x.SetLineColor(ROOT.kBlack); line0x.SetLineStyle(2); line0x.SetLineWidth(2); line0x.Draw("SAME")
line1x = ROOT.TLine(offset_x, 0, offset_x, h_diff_x.GetMaximum())
line1x.SetLineColor(ROOT.kRed); line1x.SetLineStyle(2); line1x.SetLineWidth(2); line1x.Draw("SAME")

# Y Offset Fit
c1.cd(2)
hy = h_diff_y.GetValue()
hy.SetFillColor(ROOT.kOrange + 1)
hy.Draw()

# Fit Y
fit_y = ROOT.TF1("fit_y", "gaus", -0.2, 0.2)
fit_y.SetParameters(h_diff_y.GetMaximum(), h_diff_y.GetMean(), h_diff_y.GetRMS())
hy.Fit(fit_y, "RQ")
offset_y = fit_y.GetParameter(1)
offset_y_err = fit_y.GetParError(1)
sigma_y = fit_y.GetParameter(2)

line0y = ROOT.TLine(0, 0, 0, h_diff_y.GetMaximum())
line0y.SetLineColor(ROOT.kBlack); line0y.SetLineStyle(2); line0y.SetLineWidth(2); line0y.Draw("SAME")
line1y = ROOT.TLine(offset_y, 0, offset_y, h_diff_y.GetMaximum())
line1y.SetLineColor(ROOT.kRed); line1y.SetLineStyle(2); line1y.SetLineWidth(2); line1y.Draw("SAME")

c1.Update()

 
print("=" * 50)
print("GAUSSIAN FIT (offset between IN and OUT detector systems):")
print(f"\tOffset x = {offset_x:.4f} +/- {offset_x_err:.4f} mm  (sigma = {sigma_x:.4f} mm)")
print(f"\t\tChi^2/NDF = {(fit_x.GetChisquare() / fit_x.GetNDF()):.3f}")
print(f"\t\tSkewness = {h_diff_x.GetSkewness():.4f}")
print(f"\tOffset y = {offset_y:.4f} +/- {offset_y_err:.4f} mm  (sigma = {sigma_y:.4f} mm)")
print(f"\t\tChi^2/NDF = {(fit_y.GetChisquare() / fit_y.GetNDF()):.3f}")
print(f"\t\tSkewness = {h_diff_y.GetSkewness():.4f}")

h_prof_x = df_air.Profile1D(
    ("h_prof_x", "x_{out} vs x_{in} (air only); x_{in} [mm]; x_{out} [mm]",
     100, -5, 15), "Tracks.d0_x", "Tracks.d0Out_x")
h_prof_y = df_air.Profile1D(
    ("h_prof_y", "y_{out} vs y_{in} (air only); y_{in} [mm]; y_{out} [mm]",
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


c2 = ROOT.TCanvas("c2", "Detector alignment - profile fits", 1400, 700)
c2.Divide(2, 1)
c2.cd(1)
lin_x.SetLineColor(ROOT.kRed)
h_prof_x_val.Draw()
lin_x.Draw("SAME")
h_prof_x_val.Draw("SAME")
c2.cd(2)
lin_y.SetLineColor(ROOT.kRed)
h_prof_y_val.Draw()
lin_y.Draw("SAME")
h_prof_y_val.Draw("SAME")
c2.Update()

# c1.SaveAs(f"./plots_{file_num}/alignment_check.png")
# c1.SaveAs(f"./plots_{file_num}/alignment_check.pdf")
c2.SaveAs(f"./plots_{file_num}/linear_check.png")
c2.SaveAs(f"./plots_{file_num}/linear_check.pdf")
 
# # 2D scatter of the air-only sample, to visually confirm these
# # tracks really are outside the crystal footprint
# h_d0_air = df_air.Histo2D(
#     ("h_d0_air", "Air-only tracks (incoming); d0_x [mm]; d0_y [mm]",
#      200, -20, 20, 200, -20, 20), "Tracks.d0_x", "Tracks.d0_y")
# c3 = ROOT.TCanvas("c3", "Air-only sample sanity check", 900, 800)
# h_d0_air.Draw("COLZ")
# box = ROOT.TBox(x_min, y_min, x_max, y_max)
# box.SetLineColor(ROOT.kRed); box.SetLineWidth(2); box.SetFillStyle(0)
# box.Draw("SAME")
# c3.Update()


 

