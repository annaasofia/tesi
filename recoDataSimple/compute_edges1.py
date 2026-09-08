import ROOT
import math
import sys
import os
import numpy as np
from array import array
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import Rectangle
import plotting_utils as pu

ROOT.ROOT.EnableImplicitMT() 
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetPalette(ROOT.kBird)

file = 8430
filename = "data/recoDataSimple_" + str(file) + "_xtalMerging.root"
files = ["data/recoDataSimple_8430_xtalMerging.root", "data/recoDataSimple_8431_xtalMerging.root"]

os.makedirs(f"plots_{file}_edges1", exist_ok=True)

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

# FILTERING the data (single tracks) and conversion from rad to urad
df_phys = df.Filter("SingleTrack == 1")
df_phys = df_phys.Define("thetaIn_x", "Tracks.thetaIn_x * 1e6").Define("Deltatheta_x", "(Tracks.thetaOut_x - Tracks.thetaIn_x) * 1e6")
print("="*50)

# HISTOGRAMS PRE SPATIAL CUT
h_d0_xy = df_phys.Histo2D(("h_d0_xy", "Incoming x vs y of all particles; x [mm]; y [mm]", 1700, -10, 20, 1700, -20, 20), "Tracks.d0_x", "Tracks.d0_y")
h_d0_Out_xy = df_phys.Histo2D(("h_d0_Out_xy", "Outgoing x vs y of all particles; x [mm]; y [mm]", 1700, -10, 20, 1700, -20, 20), "Tracks.d0Out_x", "Tracks.d0Out_y")

# FILTER: spatial cut (d0_x and d0_y within the crystal area)
# as first thing i choose a cut on Deltatheta_x to select only channeled particles
preliminary_cut = 0
h_defl_before = df_phys.Histo1D(("h_defl_before", "Angular Deflection; #Delta#theta_{x} [#murad]; Counts", 500, -2000, max_value), "Deltatheta_x") 
h_defl_val = h_defl_before.GetValue()

pre_fit = ROOT.TF1("pre_fit", "gaus", deflection_peak - 500, deflection_peak + 500)
h_defl_val.Fit(pre_fit, "RQ0") # R=Range, Q=Quiet, 0=NoDraw
pre_mean = pre_fit.GetParameter(1); pre_sigma = pre_fit.GetParameter(2)
preliminary_cut = pre_mean - 3.0 * pre_sigma
# print(f"\tPreliminary cut on Deltatheta_x to select channeled particles: {preliminary_cut:.2f} urad (mean = {pre_mean:.2f}, sigma = {pre_sigma:.2f})")
df_cut = df_phys.Filter(f"Deltatheta_x > {preliminary_cut}", "Preliminary cut on Deltatheta_x to select channeled particles")

# i look at the beam profile of d0 of those who channeled, and find the region (width x length)
h_d0_xy_ch = df_cut.Histo2D(("h_d0_xy_ch", "Incoming beam - channeled particles; x [mm]; y [mm]", 5000, -2, 3, 1700, -8, 9), "Tracks.d0_x", "Tracks.d0_y")
h_d0_Out_xy_ch = df_cut.Histo2D(("h_d0_Out_xy_ch", "Outgoing beam - channeled particles; x [mm]; y [mm]", 5000, -2, 3, 1700, -8, 9), "Tracks.d0Out_x", "Tracks.d0Out_y")

h_d0_x_ch = df_cut.Histo1D(("h_d0_x_ch", "x of channeled particles; x [mm]; Counts", 5000, -2, 3), "Tracks.d0_x")
h_d0_y_ch = df_cut.Histo1D(("h_d0_y_ch", "y of channeled particles; y [mm]; Counts", 1700, -8, 9), "Tracks.d0_y")
h_d0_Out_x_ch = df_cut.Histo1D(("h_d0_Out_x_ch", "x of channeled particles; x [mm]; Counts", 5000, -2, 3), "Tracks.d0Out_x")
h_d0_Out_y_ch = df_cut.Histo1D(("h_d0_Out_y_ch", "y of channeled particles; y [mm]; Counts", 1700, -8, 9), "Tracks.d0Out_y")


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

# Definisci tutte le operazioni lazy PRIMA di estrarre i valori
mean_node_in_x  = df_cut.Mean("Tracks.d0Err_x")
mean_node_in_y  = df_cut.Mean("Tracks.d0Err_y")
mean_node_out_x = df_cut.Mean("Tracks.d0OutErr_x")
mean_node_out_y = df_cut.Mean("Tracks.d0OutErr_y")
# La prima chiamata a GetValue() fa partire l'unico loop sui dati necessario
mean_d0Err_in_x  = mean_node_in_x.GetValue()
mean_d0Err_in_y  = mean_node_in_y.GetValue()
mean_d0Err_out_x = mean_node_out_x.GetValue()
mean_d0Err_out_y = mean_node_out_y.GetValue()

# uncertainty due to binning choice in the sliding window (uniform distribution over the bin width -> variance = bin_width^2/12)
bin_w_y = h_y_in.GetBinWidth(1)
bin_w_x = h_x_restricted.GetBinWidth(1)
disc_err_y = bin_w_y / math.sqrt(12)
disc_err_x = bin_w_x / math.sqrt(12)

# quadrature combination: detector resolution + discretization
y_min_err = math.sqrt(mean_d0Err_in_y**2 + disc_err_y**2)
y_max_err = math.sqrt(mean_d0Err_out_y**2 + disc_err_y**2)
x_min_err = math.sqrt(mean_d0Err_in_x**2 + disc_err_x**2)
x_max_err = math.sqrt(mean_d0Err_out_x**2 + disc_err_x**2)

print(f"x = [{x_min:.4f} ± {x_min_err:.4f}, {x_max:.4f} ± {x_max_err:.4f}] mm")
print(f"y = [{y_min:.4f} ± {y_min_err:.4f}, {y_max:.4f} ± {y_max_err:.4f}] mm")

# x1, y1, x2, y2 = -2.0, -8.0, 2.0, 8.0
# box1 = ROOT.TBox(x1, y1, x2, y2); box1.SetLineColor(ROOT.kMagenta); box1.SetLineWidth(2); box1.SetFillStyle(0); box1.Draw("SAME")
# box2 = ROOT.TBox(x1, y1, x2, y2); box2.SetLineColor(ROOT.kMagenta); box2.SetLineWidth(2); box2.SetFillStyle(0); box2.Draw("SAME")
# box3 = ROOT.TBox(x_min, y_min, x_max, y_max); box3.SetLineColor(ROOT.kRed); box3.SetLineWidth(2); box3.SetFillStyle(0); box3.Draw("SAME")
# box4 = ROOT.TBox(x_min+delta_x, y_min, x_max+delta_x, y_max); box4.SetLineColor(ROOT.kRed); box4.SetLineWidth(2); box4.SetFillStyle(0); box4.Draw("SAME")
# box5 = ROOT.TBox(x_min, y_min, x_max, y_max); box5.SetLineColor(ROOT.kRed); box5.SetLineWidth(1); box5.SetLineStyle(2); box5.SetFillStyle(0); box5.Draw("SAME")

# Creazione di copie a bassa risoluzione ESCLUSIVAMENTE per Matplotlib
# Riduciamo i bin di un fattore 5x5 e 10x10 per non bloccare il rendering
plot_h_d0_xy = h_d0_xy.GetValue().Clone("plot_h_d0_xy")
plot_h_d0_xy.Rebin2D(5, 5)

plot_h_d0_Out_xy = h_d0_Out_xy.GetValue().Clone("plot_h_d0_Out_xy")
plot_h_d0_Out_xy.Rebin2D(5, 5)

plot_h_d0_xy_ch = h_d0_xy_ch.GetValue().Clone("plot_h_d0_xy_ch")
plot_h_d0_xy_ch.Rebin2D(10, 10)

plot_h_d0_Out_xy_ch = h_d0_Out_xy_ch.GetValue().Clone("plot_h_d0_Out_xy_ch")
plot_h_d0_Out_xy_ch.Rebin2D(10, 10)

fig0, axs0 = plt.subplots(2, 2, figsize=(14, 10))
# 1) Incoming all
pu.plot_histo2d(plot_h_d0_xy, ax=axs0[0,0], xlabel="x [mm]", ylabel="y [mm]", title="Incoming d0_x vs d0_y (All)")
axs0[0,0].add_patch(Rectangle((-2.0, -8.0), 4.0, 16.0, fill=False, edgecolor='magenta', lw=2))
# 2) Outgoing all
pu.plot_histo2d(plot_h_d0_Out_xy, ax=axs0[0,1], xlabel="x [mm]", ylabel="y [mm]", title="Outgoing d0_x vs d0_y (All)")
axs0[0,1].add_patch(Rectangle((-2.0, -8.0), 4.0, 16.0, fill=False, edgecolor='magenta', lw=2))
# 3) Incoming channeled
pu.plot_histo2d(plot_h_d0_xy_ch, ax=axs0[1,0], xlabel="x [mm]", ylabel="y [mm]", title="Incoming d0_x vs d0_y (Channeled)")
axs0[1,0].add_patch(Rectangle((x_min, y_min), x_max-x_min, y_max-y_min, fill=False, edgecolor='red', lw=2))
# 4) Outgoing channeled
pu.plot_histo2d(plot_h_d0_Out_xy_ch, ax=axs0[1,1], xlabel="x [mm]", ylabel="y [mm]", title="Outgoing d0_x vs d0_y (Channeled)")
axs0[1,1].add_patch(Rectangle((x_min+delta_x, y_min), x_max-x_min, y_max-y_min, fill=False, edgecolor='red', lw=2))
axs0[1,1].add_patch(Rectangle((x_min, y_min), x_max-x_min, y_max-y_min, fill=False, edgecolor='red', lw=1, ls='--'))
fig0.tight_layout()
fig0.savefig(f"plots_{file}_edges1/spatial_cut_box.pdf")
fig0.savefig(f"plots_{file}_edges1/spatial_cut_box.png")
plt.close(fig0)

fig0b, ax0b = plt.subplots(figsize=(8, 7))
pu.plot_histo2d(plot_h_d0_xy, ax=ax0b, xlabel="x [mm]", ylabel="y [mm]", title="Incoming beam of all particles")
ax0b.add_patch(Rectangle((x_min, y_min), x_max-x_min, y_max-y_min, fill=False, edgecolor='red', lw=2))
fig0b.savefig(f"plots_{file}_edges1/spatial_cut2_box.pdf")
fig0b.savefig(f"plots_{file}_edges1/spatial_cut2_box.png")
plt.close(fig0b)

# c1 = ROOT.TCanvas("c1", "d0_x and d0_y of channeled particles", 1000, 900)
# pad_center = ROOT.TPad("pad_center", "pad_center", 0, 0, 0.65, 0.65)
# pad_center.Draw()
# pad_top = ROOT.TPad("pad_top", "pad_top", 0.0, 0.60, 0.65, 1.0)
# pad_top.Draw()
# pad_right = ROOT.TPad("pad_right", "pad_right", 0.60, 0.0, 1.0, 0.65)
# pad_right.Draw()
# h_d0_Out_xy_ch.SetTitle(""); pad_center.cd(); h_d0_Out_xy_ch.Draw("COL"); box4.Draw("SAME")
# h_d0_Out_x_ch.SetTitle(""); h_d0_Out_x_ch.GetXaxis().SetTitle(""); pad_top.cd(); h_d0_Out_x_ch.SetFillColor(ROOT.kAzure-3); h_d0_Out_x_ch.Draw("BAR X+")
# l1 = ROOT.TLine(x_min+delta_x, 0, x_min+delta_x, h_d0_Out_x_ch.GetMaximum()); l1.SetLineColor(ROOT.kRed); l1.SetLineStyle(2); l1.SetLineWidth(2); l1.Draw("SAME")
# l2 = ROOT.TLine(x_max+delta_x, 0, x_max+delta_x, h_d0_Out_x_ch.GetMaximum()); l2.SetLineColor(ROOT.kRed); l2.SetLineStyle(2); l2.SetLineWidth(2); l2.Draw("SAME")
# h_d0_Out_y_ch.SetTitle(""); h_d0_Out_y_ch.GetXaxis().SetTitle(""); pad_right.cd(); h_d0_Out_y_ch.SetFillColor(ROOT.kAzure-3); h_d0_Out_y_ch.Draw("HBAR Y+")
# l3 = ROOT.TLine(0, y_min, h_d0_Out_y_ch.GetMaximum(), y_min); l3.SetLineColor(ROOT.kRed); l3.SetLineStyle(2); l3.SetLineWidth(2); l3.Draw("SAME")
# l4 = ROOT.TLine(0, y_max, h_d0_Out_y_ch.GetMaximum(), y_max); l4.SetLineColor(ROOT.kRed); l4.SetLineStyle(2); l4.SetLineWidth(2); l4.Draw("SAME")
# g_fit_out = ROOT.TGraph(); step = (y_max - y_min) / 1700
# for i in range(1700):
#     y_val = y_min + i * step
#     x_val = f_gaus_out.Eval(y_val) 
#     g_fit_out.SetPoint(i, x_val, y_val) 
# g_fit_out.SetLineColor(ROOT.kBlue+2); g_fit_out.SetLineWidth(2); g_fit_out.Draw("L SAME")
# c1.Update()

fig1 = plt.figure(figsize=(10, 10))
gs1 = gridspec.GridSpec(2, 2, width_ratios=(3, 1), height_ratios=(1, 3), wspace=0.1, hspace=0.1)

ax_main1 = fig1.add_subplot(gs1[1, 0])
ax_top1 = fig1.add_subplot(gs1[0, 0], sharex=ax_main1)
ax_right1 = fig1.add_subplot(gs1[1, 1], sharey=ax_main1)

# Mappa 2D principale
pu.plot_histo2d(plot_h_d0_Out_xy_ch, ax=ax_main1, xlabel="x [mm]", ylabel="y [mm]")
ax_main1.add_patch(Rectangle((x_min+delta_x, y_min), x_max-x_min, y_max-y_min, fill=False, edgecolor='red', lw=2))
# Proiezione X (Top)
pu.plot_histo1d(h_d0_Out_x_ch, ax=ax_top1, style="fill", color="darkblue", ylabel="Counts")
ax_top1.axvline(x_min+delta_x, color='red', ls='--', lw=2)
ax_top1.axvline(x_max+delta_x, color='red', ls='--', lw=2)
ax_top1.tick_params(labelbottom=False); ax_top1.set_xlabel("")
# Proiezione Y (Right) - Estratta manualmente per usare l'orientamento orizzontale
centers_y1, contents_y1, yerr_y1, edges_y1 = pu.th1_to_arrays(h_d0_Out_y_ch)
ax_right1.stairs(contents_y1, edges_y1, fill=True, color="darkblue", orientation='horizontal')
ax_right1.axhline(y_min, color='red', ls='--', lw=2)
ax_right1.axhline(y_max, color='red', ls='--', lw=2)
# Fit su asse Y destro invertendo X e Y per adattarsi all'orientamento orizzontale
fy_x, fy_y = pu.tf1_to_curve(f_gaus_out, y_min, y_max)
ax_right1.plot(fy_y, fy_x, color='darkblue', lw=2) 
ax_right1.tick_params(labelleft=False); ax_right1.set_xlabel("Counts")

fig1.savefig(f"plots_{file}_edges1/d0_Out_projxy_box.pdf")
fig1.savefig(f"plots_{file}_edges1/d0_Out_projxy_box.png")
plt.close(fig1)

# c2 = ROOT.TCanvas("c2", "d0_x and d0_y of channeled particles", 1000, 900)
# pad_center = ROOT.TPad("pad_center", "pad_center", 0, 0, 0.65, 0.65)
# pad_center.Draw()
# pad_top = ROOT.TPad("pad_top", "pad_top", 0.0, 0.60, 0.65, 1.0)
# pad_top.Draw()
# pad_right = ROOT.TPad("pad_right", "pad_right", 0.60, 0.0, 1.0, 0.65)
# pad_right.Draw()
# h_d0_xy_ch.SetTitle(""); pad_center.cd(); h_d0_xy_ch.Draw("COL"); box3.Draw("SAME")
# h_d0_x_ch.SetTitle(""); h_d0_x_ch.GetXaxis().SetTitle(""); pad_top.cd(); h_d0_x_ch.SetFillColor(ROOT.kAzure-3); h_d0_x_ch.Draw("BAR X+"); 
# l5 = ROOT.TLine(x_min, 0, x_min, h_d0_x_ch.GetMaximum()); l5.SetLineColor(ROOT.kRed); l5.SetLineStyle(2); l5.SetLineWidth(2); l5.Draw("SAME")
# l6 = ROOT.TLine(x_max, 0, x_max, h_d0_x_ch.GetMaximum()); l6.SetLineColor(ROOT.kRed); l6.SetLineStyle(2); l6.SetLineWidth(2); l6.Draw("SAME")
# h_d0_y_ch.SetTitle(""); h_d0_y_ch.GetXaxis().SetTitle(""); pad_right.cd(); h_d0_y_ch.SetFillColor(ROOT.kAzure-3); h_d0_y_ch.Draw("HBAR Y+");
# g_fit_in = ROOT.TGraph(); step = (y_max - y_min) / 1700
# for i in range(1700):
#     y_val = y_min + i * step
#     x_val = f_gaus_in.Eval(y_val) 
#     g_fit_in.SetPoint(i, x_val, y_val) 
# g_fit_in.SetLineColor(ROOT.kBlue+2); g_fit_in.SetLineWidth(2); g_fit_in.Draw("L SAME")
# l3.Draw("SAME"); l4.Draw("SAME")
# c2.Update()

fig2 = plt.figure(figsize=(10, 10))
gs2 = gridspec.GridSpec(2, 2, width_ratios=(3, 1), height_ratios=(1, 3), wspace=0.1, hspace=0.1)

ax_main2 = fig2.add_subplot(gs2[1, 0])
ax_top2 = fig2.add_subplot(gs2[0, 0], sharex=ax_main2)
ax_right2 = fig2.add_subplot(gs2[1, 1], sharey=ax_main2)

# Mappa 2D principale
pu.plot_histo2d(plot_h_d0_xy_ch, ax=ax_main2, xlabel="x [mm]", ylabel="y [mm]")
ax_main2.add_patch(Rectangle((x_min, y_min), x_max-x_min, y_max-y_min, fill=False, edgecolor='red', lw=2))
# Proiezione X (Top)
pu.plot_histo1d(h_d0_x_ch, ax=ax_top2, style="fill", color="darkblue", ylabel="Counts")
ax_top2.axvline(x_min, color='red', ls='--', lw=2)
ax_top2.axvline(x_max, color='red', ls='--', lw=2)
ax_top2.tick_params(labelbottom=False); ax_top2.set_xlabel("")
# Proiezione Y (Right) - Estratta manualmente per l'orientamento orizzontale
centers_y2, contents_y2, yerr_y2, edges_y2 = pu.th1_to_arrays(h_d0_y_ch)
ax_right2.stairs(contents_y2, edges_y2, fill=True, color="darkblue", orientation='horizontal')
ax_right2.axhline(y_min, color='red', ls='--', lw=2)
ax_right2.axhline(y_max, color='red', ls='--', lw=2)
# Fit su asse Y destro 
fy_x2, fy_y2 = pu.tf1_to_curve(f_gaus_in, y_min, y_max)
ax_right2.plot(fy_y2, fy_x2, color='darkblue', lw=2)
ax_right2.tick_params(labelleft=False); ax_right2.set_xlabel("Counts")

fig2.savefig(f"plots_{file}_edges1/d0_projxy_box.pdf")
fig2.savefig(f"plots_{file}_edges1/d0_projxy_box.png")
plt.close(fig2)

# now that we have the crystal area/position, we can apply the spatial cut to the entire dataframe
spatial_cut = f"Tracks.d0_x > {x_min} && Tracks.d0_x < {x_max} && Tracks.d0_y > {y_min} && Tracks.d0_y < {y_max}"
df_phys = df_phys.Filter(spatial_cut, "Spatial Cut (Crystal Area)")

print(f"All matplotlib plots saved under plots_{file}_edges1/")



