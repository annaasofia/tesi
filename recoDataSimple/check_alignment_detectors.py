import ROOT
import sys

ROOT.ROOT.EnableImplicitMT()
ROOT.gStyle.SetOptStat(1111)

file_num = 8430
filename = f"recoDataSimple_{file_num}_xtalMerging.root"

print(f"Analyzing {filename} for alignment offset...")

x_min, x_max = 0.0, 0.0
y_min, y_max = 0.0, 0.0
 
if x_min == x_max or y_min == y_max:
    print("WARNING: crystal box not set (x_min/x_max/y_min/y_max). "
          "Paste the values printed by the main macro before running this.")
    sys.exit(1)

margin_x = 1.5  # mm
margin_y = 1.5  # mm


df = ROOT.RDataFrame("simpleEvent", filename)
df = df.Filter("SingleTrack == 1", "Single Track")

# Spatial Cut: Select particles that definitely missed the crystal.
air_cut = "Tracks.d0_x > 5.0 || Tracks.d0_x < -4.0"
df_air = df.Filter(air_cut, "Air only (missed crystal)")

# 4. Calculate the difference between downstream and upstream projections
# If perfectly aligned, d0Out - d0 should be centered exactly at 0.
df_diff = df_air.Define("diff_x", "Tracks.d0Out_x - Tracks.d0_x") \
                .Define("diff_y", "Tracks.d0Out_y - Tracks.d0_y")

# 5. Create Histograms
# A range of +/- 2 mm should be plenty to catch small sub-millimeter offsets
h_diff_x = df_diff.Histo1D(("h_diff_x", "Alignment Offset X (Air Particles); d0Out_x - d0_x [mm]; Counts", 200, -2.0, 2.0), "diff_x")
h_diff_y = df_diff.Histo1D(("h_diff_y", "Alignment Offset Y (Air Particles); d0Out_y - d0_y [mm]; Counts", 200, -2.0, 2.0), "diff_y")

# 6. Fit with Gaussian to extract the exact offset
c1 = ROOT.TCanvas("c1", "Detector Alignment Check", 1200, 600)
c1.Divide(2, 1)

# X Offset Fit
c1.cd(1)
hx = h_diff_x.GetValue()
hx.SetFillColor(ROOT.kAzure - 3)
hx.Draw()

# Fit X
fit_x = ROOT.TF1("fit_x", "gaus", -1.0, 1.0)
hx.Fit(fit_x, "RQ")
offset_x = fit_x.GetParameter(1)
print(f"--> Measured X Offset (delta_x): {offset_x:.4f} mm")

# Y Offset Fit
c1.cd(2)
hy = h_diff_y.GetValue()
hy.SetFillColor(ROOT.kOrange + 1)
hy.Draw()

# Fit Y
fit_y = ROOT.TF1("fit_y", "gaus", -1.0, 1.0)
hy.Fit(fit_y, "RQ")
offset_y = fit_y.GetParameter(1)
print(f"--> Measured Y Offset (delta_y): {offset_y:.4f} mm")

c1.Update()
c1.SaveAs(f"alignment_check_{file_num}.png")

# Keep the program open to view the canvas
input("Press Enter to exit...")