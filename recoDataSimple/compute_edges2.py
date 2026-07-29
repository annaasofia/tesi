import ROOT
import math
import sys
import numpy as np
from array import array



ROOT.ROOT.EnableImplicitMT() 

# file = input("File number: ")
file = 8430
filename = "recoDataSimple_" + str(file) + "_xtalMerging.root"
files = ["recoDataSimple_8430_xtalMerging.root", "recoDataSimple_8431_xtalMerging.root"]

ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetPalette(ROOT.kBird)

# VARIABLES
theta_L = 0 # Lindhard angle for TCCP and TCCPA in urad
width = 0 # width of the crystal in mm (for spatial cut)
height = 0 # height of the crystal in mm (for spatial cut)
delta_x = 0 # shift in x (for spatial cut)
theta_0 = 0 # centro della tua isola di channeling - looking at 2D "h scan" find on which x (thetaIn) the channeling isle is centered
max_value = 0 # histograms range
deflection_peak = 0 # urad

if file in [8430, 8431, 8655, 8656]:
    theta_L = 12.992
    deflection_peak = 6010.0
    width = 12.8
    height = 2
    delta_x = 0.5*74*6.01*pow(10,-3) # 0.222 mm
    max_value = 8000
elif file in [8650]:
    theta_L = 14.294
    deflection_peak = 6010.0
    width = 12.8
    height = 2
    delta_x = 0.5*74*6.01*pow(10,-3) # 0.222 mm
    max_value = 8000
else:
    print("Run not found...")
    sys.exit(1)

# ROOT DATA FRAME
df = ROOT.RDataFrame("simpleEvent", filename)
# df = ROOT.RDataFrame("simpleEvent", files)
print("="*50)
print(f"Analyzing {filename} ...")

# FILTERING the data: single tracks and conversion from rad to urad
df_phys = df.Filter("SingleTrack == 1")
df_phys = df_phys.Define("thetaIn_x", "Tracks.thetaIn_x * 1e6").Define("Deltatheta_x", "(Tracks.thetaOut_x - Tracks.thetaIn_x) * 1e6")
print("="*50)

# HISTOGRAMS PRE SPATIAL CUT
h_d0_xy = df_phys.Histo2D(("h_d0_xy", "Incoming d0_x vs d0_y of all particles; d0_x [mm]; d0_y [mm]", 1700, -10, 20, 1700, -20, 20), "Tracks.d0_x", "Tracks.d0_y")
h_d0_Out_xy = df_phys.Histo2D(("h_d0_Out_xy", "Outgoing d0_x vs d0_y of all particles; d0_x [mm]; d0_y [mm]", 1700, -10, 20, 1700, -20, 20), "Tracks.d0Out_x", "Tracks.d0Out_y")
c0 = ROOT.TCanvas("c0", "d0_x and d0_y of all particles", 1700, 900)
c0.Divide(2, 2)
c0.cd(1); h_d0_xy.Draw("COLZ")
x1, y1, x2, y2 = -2.0, -8.0, 2.0, 8.0
box1 = ROOT.TBox(x1, y1, x2, y2); box1.SetLineColor(ROOT.kMagenta); box1.SetLineWidth(2); box1.SetFillStyle(0); box1.Draw("SAME")
c0.cd(2); h_d0_Out_xy.Draw("COLZ")
box2 = ROOT.TBox(x1, y1, x2, y2); box2.SetLineColor(ROOT.kMagenta); box2.SetLineWidth(2); box2.SetFillStyle(0); box2.Draw("SAME")

# FILTER 3: spatial cut (d0_x and d0_y within the crystal area)
# as first thing i choose a cut on Deltatheta_x to select only channeled particles
preliminary_cut = 0
h_defl_before = df_phys.Histo1D(("h_defl_before", "Angular Deflection; #Delta#theta_{x} [#murad]; Counts", 500, -2000, max_value), "Deltatheta_x") 
h_defl_val = h_defl_before.GetValue()

pre_fit = ROOT.TF1("pre_fit", "gaus", deflection_peak - 500, deflection_peak + 500)
h_defl_val.Fit(pre_fit, "RQ0") # R=Range, Q=Quiet, 0=NoDraw
pre_mean = pre_fit.GetParameter(1); pre_sigma = pre_fit.GetParameter(2)
preliminary_cut = pre_mean - 3.0 * pre_sigma
df_cut = df_phys.Filter(f"Deltatheta_x > {preliminary_cut}", "Preliminary cut on Deltatheta_x to select channeled particles")

# i look at the beam profile of d0 of those who channeled, and find the region (width x length)
h_d0_xy_ch = df_cut.Histo2D(("h_d0_xy_ch", "Incoming beam - channeled particles; d0_x [mm]; d0_y [mm]", 5000, -2, 3, 1700, -8, 9), "Tracks.d0_x", "Tracks.d0_y")
h_d0_Out_xy_ch = df_cut.Histo2D(("h_d0_Out_xy_ch", "Outgoing beam - channeled particles; d0Out_x [mm]; d0Out_y [mm]", 5000, -2, 3, 1700, -8, 9), "Tracks.d0Out_x", "Tracks.d0Out_y")

h_d0_x_ch = df_cut.Histo1D(("h_d0_x_ch", "d0_x of channeled particles; d0_x [mm]; Counts", 5000, -2, 3), "Tracks.d0_x")
h_d0_y_ch = df_cut.Histo1D(("h_d0_y_ch", "d0_y of channeled particles; d0_y [mm]; Counts", 1700, -8, 9), "Tracks.d0_y")
h_d0_Out_x_ch = df_cut.Histo1D(("h_d0_Out_x_ch", "d0Out_x of channeled particles; d0Out_x [mm]; Counts", 5000, -2, 3), "Tracks.d0Out_x")
h_d0_Out_y_ch = df_cut.Histo1D(("h_d0_Out_y_ch", "d0Out_y of channeled particles; d0Out_y [mm]; Counts", 1700, -8, 9), "Tracks.d0Out_y")


f_gaus_in = ROOT.TF1("f_gaus_in", "gaus", -8, 9)
h_d0_y_ch.Fit(f_gaus_in, "RQ0")
y_c_in, y_c_in_err = f_gaus_in.GetParameter(1), f_gaus_in.GetParError(1)
f_gaus_out = ROOT.TF1("f_gaus_out", "gaus", -8, 9)
h_d0_Out_y_ch.Fit(f_gaus_out, "RQ0")
y_c_out, y_c_out_err = f_gaus_out.GetParameter(1), f_gaus_out.GetParError(1)

w_in, w_out = 1.0/y_c_in_err**2, 1.0/y_c_out_err**2
y_c = (y_c_in*w_in + y_c_out*w_out) / (w_in + w_out)
y_c_err = (1.0/(w_in + w_out))**0.5

print(f"mean(in) = {y_c_in:.3f} +/- {y_c_in_err:.3f}")
print(f"mean(out) = {y_c_out:.3f} +/- {y_c_out_err:.3f}")

if file in [0]:
    print(f"mean = {y_c:.3f} +/- {y_c_err:.3f}")
    y_min, y_max = y_c - width / 2.0, y_c + width / 2.0

if file in [8430, 8431, 8650, 8655, 8656]:
    # sliding window method along y
    h_y_in = h_d0_y_ch.GetValue()
    max_particles = -1
    best_y_min_in = 0
    for i in range(1, h_y_in.GetNbinsX() + 1):

        current_y_start = h_y_in.GetBinLowEdge(i)
        current_y_end = current_y_start + width
        
        bin_start = i
        bin_end = h_y_in.FindBin(current_y_end)
        
        particles_in_window = h_y_in.Integral(bin_start, bin_end)
        
        if particles_in_window > max_particles:
            max_particles = particles_in_window
            best_y_min_in = current_y_start

    y_min_in = best_y_min_in
    y_max_in = best_y_min_in + width
    print(f"upstream window: [{y_min_in:.3f}:{y_max_in:.3f}]")

    h_y_out = h_d0_Out_y_ch.GetValue()
    max_particles = -1
    best_y_min_out = 0
    for i in range(1, h_y_out.GetNbinsX() + 1):

        current_y_start = h_y_out.GetBinLowEdge(i)
        current_y_end = current_y_start + width
        
        bin_start = i
        bin_end = h_y_out.FindBin(current_y_end)
        
        particles_in_window = h_y_out.Integral(bin_start, bin_end)
        
        if particles_in_window > max_particles:
            max_particles = particles_in_window
            best_y_min_out = current_y_start

    y_min_out = best_y_min_out
    y_max_out = best_y_min_out + width
    print(f"downstream window: [{y_min_out:.3f}:{y_max_out:.3f}]")

    y_min, y_max = (y_min_in + y_min_out)/2.0, (y_max_in + y_max_out)/2.0

# sliding window method along x since the distribution is not symmetric neither gaussian
h_xy = h_d0_Out_xy_ch.GetValue() if hasattr(h_d0_Out_xy_ch, 'GetValue') else h_d0_Out_xy_ch
y_bin_min = h_xy.GetYaxis().FindBin(y_min)
y_bin_max = h_xy.GetYaxis().FindBin(y_max)
h_x_restricted = h_xy.ProjectionX("h_x_restricted", y_bin_min, y_bin_max)

max_particles = -1
best_x_min = 0
for i in range(1, h_x_restricted.GetNbinsX() + 1):

    current_x_start = h_x_restricted.GetBinLowEdge(i)
    current_x_end = current_x_start + height
    
    bin_start = i
    bin_end = h_x_restricted.FindBin(current_x_end)
    
    particles_in_window = h_x_restricted.Integral(bin_start, bin_end)
    
    if particles_in_window > max_particles:
        max_particles = particles_in_window
        best_x_min = current_x_start

best_x_max = best_x_min + height
print(f"\tTaglio ottimale in x trovato: [{best_x_min:.3f} mm, {best_x_max:.3f} mm]")

x_min, x_max = best_x_min - delta_x, best_x_max - delta_x
print(f"\tTaglio finale in x (dopo shift di {delta_x:.3f} mm): x = [{x_min:.4f} mm, {x_max:.4f} mm], y = [{y_min:.4f} mm, {y_max:.4f} mm]")

c0.cd(3); h_d0_xy_ch.Draw("COLZ")
box3 = ROOT.TBox(x_min, y_min, x_max, y_max); box3.SetLineColor(ROOT.kRed); box3.SetLineWidth(2); box3.SetFillStyle(0); box3.Draw("SAME")
c0.cd(4); h_d0_Out_xy_ch.Draw("COLZ")
box4 = ROOT.TBox(x_min+delta_x, y_min, x_max+delta_x, y_max); box4.SetLineColor(ROOT.kRed); box4.SetLineWidth(2); box4.SetFillStyle(0); box4.Draw("SAME")
box5 = ROOT.TBox(x_min, y_min, x_max, y_max); box5.SetLineColor(ROOT.kRed); box5.SetLineWidth(1); box5.SetLineStyle(2); box5.SetFillStyle(0); box5.Draw("SAME")
c0.Update()

c0b = ROOT.TCanvas("c0b", "d0_x and d0_y of channeled particles", 1200, 900)
h_d0_xy.Draw("COLZ"); h_d0_xy.SetTitle("Incoming beam of all particles; d0_x [mm]; d0_y [mm]")
box_selected = ROOT.TBox(x_min, y_min, x_max, y_max); box_selected.SetLineColor(ROOT.kRed); box_selected.SetLineWidth(2); box_selected.SetFillStyle(0); box_selected.Draw("SAME")
c0b.Update()

c1 = ROOT.TCanvas("c1", "d0_x and d0_y of channeled particles", 1000, 900)
pad_center = ROOT.TPad("pad_center", "pad_center", 0, 0, 0.65, 0.65)
pad_center.Draw()
pad_top = ROOT.TPad("pad_top", "pad_top", 0.0, 0.60, 0.65, 1.0)
pad_top.Draw()
pad_right = ROOT.TPad("pad_right", "pad_right", 0.60, 0.0, 1.0, 0.65)
pad_right.Draw()
h_d0_Out_xy_ch.SetTitle(""); pad_center.cd(); h_d0_Out_xy_ch.Draw("COL"); box4.Draw("SAME")
h_d0_Out_x_ch.SetTitle(""); h_d0_Out_x_ch.GetXaxis().SetTitle(""); pad_top.cd(); h_d0_Out_x_ch.SetFillColor(ROOT.kAzure-3); h_d0_Out_x_ch.Draw("BAR X+")
l1 = ROOT.TLine(x_min+delta_x, 0, x_min+delta_x, h_d0_Out_x_ch.GetMaximum()); l1.SetLineColor(ROOT.kRed); l1.SetLineStyle(2); l1.SetLineWidth(2); l1.Draw("SAME")
l2 = ROOT.TLine(x_max+delta_x, 0, x_max+delta_x, h_d0_Out_x_ch.GetMaximum()); l2.SetLineColor(ROOT.kRed); l2.SetLineStyle(2); l2.SetLineWidth(2); l2.Draw("SAME")
h_d0_Out_y_ch.SetTitle(""); h_d0_Out_y_ch.GetXaxis().SetTitle(""); pad_right.cd(); h_d0_Out_y_ch.SetFillColor(ROOT.kAzure-3); h_d0_Out_y_ch.Draw("HBAR Y+")
l3 = ROOT.TLine(0, y_min, h_d0_Out_y_ch.GetMaximum(), y_min); l3.SetLineColor(ROOT.kRed); l3.SetLineStyle(2); l3.SetLineWidth(2); l3.Draw("SAME")
l4 = ROOT.TLine(0, y_max, h_d0_Out_y_ch.GetMaximum(), y_max); l4.SetLineColor(ROOT.kRed); l4.SetLineStyle(2); l4.SetLineWidth(2); l4.Draw("SAME")
g_fit_out = ROOT.TGraph(); step = (y_max - y_min) / 1700
for i in range(1700):
    y_val = y_min + i * step
    x_val = f_gaus_out.Eval(y_val) 
    g_fit_out.SetPoint(i, x_val, y_val) 
g_fit_out.SetLineColor(ROOT.kBlue+2); g_fit_out.SetLineWidth(2); g_fit_out.Draw("L SAME")
c1.Update()

c2 = ROOT.TCanvas("c2", "d0_x and d0_y of channeled particles", 1000, 900)
pad_center = ROOT.TPad("pad_center", "pad_center", 0, 0, 0.65, 0.65)
pad_center.Draw()
pad_top = ROOT.TPad("pad_top", "pad_top", 0.0, 0.60, 0.65, 1.0)
pad_top.Draw()
pad_right = ROOT.TPad("pad_right", "pad_right", 0.60, 0.0, 1.0, 0.65)
pad_right.Draw()
h_d0_xy_ch.SetTitle(""); pad_center.cd(); h_d0_xy_ch.Draw("COL"); box3.Draw("SAME")
h_d0_x_ch.SetTitle(""); h_d0_x_ch.GetXaxis().SetTitle(""); pad_top.cd(); h_d0_x_ch.SetFillColor(ROOT.kAzure-3); h_d0_x_ch.Draw("BAR X+"); 
l5 = ROOT.TLine(x_min, 0, x_min, h_d0_x_ch.GetMaximum()); l5.SetLineColor(ROOT.kRed); l5.SetLineStyle(2); l5.SetLineWidth(2); l5.Draw("SAME")
l6 = ROOT.TLine(x_max, 0, x_max, h_d0_x_ch.GetMaximum()); l6.SetLineColor(ROOT.kRed); l6.SetLineStyle(2); l6.SetLineWidth(2); l6.Draw("SAME")
h_d0_y_ch.SetTitle(""); h_d0_y_ch.GetXaxis().SetTitle(""); pad_right.cd(); h_d0_y_ch.SetFillColor(ROOT.kAzure-3); h_d0_y_ch.Draw("HBAR Y+");
g_fit_in = ROOT.TGraph(); step = (y_max - y_min) / 1700
for i in range(1700):
    y_val = y_min + i * step
    x_val = f_gaus_in.Eval(y_val) 
    g_fit_in.SetPoint(i, x_val, y_val) 
g_fit_in.SetLineColor(ROOT.kBlue+2); g_fit_in.SetLineWidth(2); g_fit_in.Draw("L SAME")
l3.Draw("SAME"); l4.Draw("SAME")
c2.Update()

# now that we have the crystal area/position, we can apply the spatial cut to the entire dataframe
spatial_cut = f"Tracks.d0_x > {x_min} && Tracks.d0_x < {x_max} && Tracks.d0_y > {y_min} && Tracks.d0_y < {y_max}"
df_phys = df_phys.Filter(spatial_cut, "Spatial Cut (Crystal Area)")