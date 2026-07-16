import ROOT
import math
import sys
import numpy as np
from array import array

ROOT.ROOT.EnableImplicitMT() 

# file = input("File number: ")
file = 8650
filename = "recoDataSimple_" + str(file) + "_xtalMerging.root"

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

if file in [8430, 8431, 8650]:
    theta_L = 12.9
    deflection_peak = 6015.0
    width = 8
    height = 2
    delta_x = 0.5*70*7*pow(10,-3) # 0.245 mm
    max_value = 8000
elif file in [8655, 8656]: 
    theta_L = 12.5
    deflection_peak = 13000
    width = 22.5
    height = 2
    delta_x = 0.5*70.5*13.3*pow(10,-3) # 0.468 mm
    max_value = 14000
else:
    print("Run not found...")
    sys.exit(1)

# ROOT DATA FRAME
df = ROOT.RDataFrame("simpleEvent", filename)
print("="*50)
print(f"Analyzing {filename} ...")

# FILTERING the data (single tracks and have entered the crystal) and conversion from rad to urad
# FILTER 1: SingleTrack == 1
df_phys = df.Filter("SingleTrack == 1")
P0 = df.Count().GetValue()
P1 = df_phys.Count().GetValue()
df_phys = df_phys.Define("thetaIn_x", "Tracks.thetaIn_x * 1e6").Define("Deltatheta_x", "(Tracks.thetaOut_x - Tracks.thetaIn_x) * 1e6")
print("="*50)
print(f"Filter 1 (SingleTrack == 1): {((P0 - P1)/P0*100):.2f}% events out of {P0} got discarded.")

# FILTER 2: good chi squared (chi2 < value)
good_chi2_value = 5
df_phys = df_phys.Filter(f"Tracks.chi2_x < {good_chi2_value} && Tracks.chi2_y < {good_chi2_value}")
P2 = df_phys.Count().GetValue()
print(f"Filter 2 (chi2 < {good_chi2_value}): {((P1 - P2)/P1*100):.2f}% events out of {P1} got discarded.")

# HISTOGRAMS PRE SPATIAL CUT
h_d0_xy = df_phys.Histo2D(("h_d0_xy", "Incoming d0_x vs d0_y of all particles; d0_x [mm]; d0_y [mm]", 200, -10, 20, 200, -20, 20), "Tracks.d0_x", "Tracks.d0_y")
h_d0_Out_xy = df_phys.Histo2D(("h_d0_Out_xy", "Outgoing d0_x vs d0_y of all particles; d0_x [mm]; d0_y [mm]", 200, -10, 20, 200, -20, 20), "Tracks.d0Out_x", "Tracks.d0Out_y")
c0 = ROOT.TCanvas("c0", "d0_x and d0_y of all particles", 1600, 900)
c0.Divide(2, 2)
c0.cd(1); h_d0_xy.Draw("COLZ")
x1, y1, x2, y2 = -2.0, -4.0, 2.0, 6.0
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
print(f"\tPreliminary cut on Deltatheta_x to select channeled particles: {preliminary_cut:.2f} urad (mean = {pre_mean:.2f}, sigma = {pre_sigma:.2f})")
df_cut = df_phys.Filter(f"Deltatheta_x > {preliminary_cut}", "Preliminary cut on Deltatheta_x to select channeled particles")

# i look at the beam profile of d0 of those who channeled, and find the region (width x length)
h_d0_xy_ch = df_cut.Histo2D(("h_d0_xy_ch", "Incoming beam - channeled particles; d0_x [mm]; d0_y [mm]", 400, -2, 2, 200, -4, 6), "Tracks.d0_x", "Tracks.d0_y")
h_d0_Out_xy_ch = df_cut.Histo2D(("h_d0_Out_xy_ch", "Outgoing beam - channeled particles; d0Out_x [mm]; d0Out_y [mm]", 400, -2, 2, 200, -4, 6), "Tracks.d0Out_x", "Tracks.d0Out_y")

h_d0_x_ch = df_cut.Histo1D(("h_d0_x_ch", "d0_x of channeled particles; d0_x [mm]; Counts", 400, -2, 2), "Tracks.d0_x")
h_d0_y_ch = df_cut.Histo1D(("h_d0_y_ch", "d0_y of channeled particles; d0_y [mm]; Counts", 200, -4, 6), "Tracks.d0_y")
h_d0_Out_x_ch = df_cut.Histo1D(("h_d0_Out_x_ch", "d0Out_x of channeled particles; d0Out_x [mm]; Counts", 400, -2, 2), "Tracks.d0Out_x")
h_d0_Out_y_ch = df_cut.Histo1D(("h_d0_Out_y_ch", "d0Out_y of channeled particles; d0Out_y [mm]; Counts", 200, -4, 6), "Tracks.d0Out_y")

y_c4 = h_d0_y_ch.GetMean(); y_c2 = h_d0_Out_y_ch.GetMean(); y_c = (y_c4 + y_c2) / 2.0 # centro del cristallo in y (mm)
y_min, y_max = y_c - width / 2.0, y_c + width / 2.0

# metodo sliding window che massimizza l'integrale per trovare margini esatti in x, piuttosto che la media essendo distribuzione asimmetrica
h_x = h_d0_Out_x_ch.GetValue()
max_particles = -1
best_x_min = 0
for i in range(1, h_x.GetNbinsX() + 1):

    current_x_start = h_x.GetBinLowEdge(i)
    current_x_end = current_x_start + height
    
    bin_start = i
    bin_end = h_x.FindBin(current_x_end)
    
    particles_in_window = h_x.Integral(bin_start, bin_end)
    
    if particles_in_window > max_particles:
        max_particles = particles_in_window
        best_x_min = current_x_start

best_x_max = best_x_min + height
print(f"\tTaglio ottimale in x trovato: [{best_x_min:.2f} mm, {best_x_max:.2f} mm]")

x_min, x_max = best_x_min - delta_x, best_x_max - delta_x
print(f"\tTaglio finale in x (dopo shift di {delta_x:.3f} mm): [{x_min:.2f} mm, {x_max:.2f} mm]")

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
h_d0_y_ch.SetTitle(""); h_d0_y_ch.GetXaxis().SetTitle(""); pad_right.cd(); h_d0_y_ch.SetFillColor(ROOT.kAzure-3); h_d0_y_ch.Draw("HBAR Y+"); l3.Draw("SAME"); l4.Draw("SAME")
c2.Update()

# now that we have the crystal area/position, we can apply the spatial cut to the entire dataframe
spatial_cut = f"Tracks.d0_x > {x_min} && Tracks.d0_x < {x_max} && Tracks.d0_y > {y_min} && Tracks.d0_y < {y_max}"
df_phys = df_phys.Filter(spatial_cut, "Spatial Cut (Crystal Area)")
P3 = df_phys.Count().GetValue()
print(f"Filter 3 (Spatial cut): {((P2 - P3)/P2*100):.2f}% events out of {P2} got discarded.")

# HISTOGRAMS PRE ANGLE CUT
h_theta_before = df_phys.Histo1D(("h_theta_before", "#theta_{in,x}; #theta_{x} [#murad]; Counts", 500, -150, 150), "thetaIn_x")
h_defl_before2 = df_phys.Histo1D(("h_defl_before2", "Angular Deflection; #Delta#theta_{x} [#murad]; Counts", 500, -2000, max_value), "Deltatheta_x") 
h_scan = df_phys.Histo2D(("h_scan", "Angular deflection of the particles as a function of the incident angle; Incident angle #theta_{In, x} [#murad]; Deflection #Delta#theta_{x} [#murad]", 500, -150, 150, 500, -2000, max_value), "thetaIn_x", "Deltatheta_x")

# FILTER 4: Lindhard cut (only keeping events that have entered the crystal with an angle within +/- 1/2 theta_L)
# FINDING THETA 0

h_theta_all = df_phys.Histo1D(("h_theta_all", "", 1000, -100, 100), "thetaIn_x")
h_theta_chan = df_phys.Filter(f"Deltatheta_x > {preliminary_cut}").Histo1D(("h_theta_chan", "", 1000, -100, 100), "thetaIn_x")
    
h_all = h_theta_all.GetValue()
h_chan = h_theta_chan.GetValue()

best_theta_0_raw = 0.0
max_eff_raw = -1.0
theta_vals = []
eff_vals = []

scan_min = -80
scan_max = 80
step = 0.2 # urad (resolution of the scan)

current_theta = scan_min
while current_theta <= scan_max:
    # cut_str = f"abs(thetaIn_x - ({current_theta})) <= {theta_L / 2.0}"
    # df_test = df_phys.Filter(cut_str)
    # n_tot_test = df_test.Count().GetValue()

    bin_min = h_all.FindBin(current_theta - theta_L / 2.0)
    bin_max = h_all.FindBin(current_theta + theta_L / 2.0)
    n_tot_test = h_all.Integral(bin_min, bin_max)
    n_ch_test  = h_chan.Integral(bin_min, bin_max)
    
    eff_test = 0.0
    if n_tot_test > 0:
        # i just count how many are in the peak, then i will do the complete computation with the final theta_0
        # n_ch_test = df_test.Filter(f"Deltatheta_x > {preliminary_cut}").Count().GetValue()
        eff_test = n_ch_test / n_tot_test * 100.0
        
        if eff_test > max_eff_raw:
            max_eff_raw = eff_test
            best_theta_0_raw = current_theta

    theta_vals.append(current_theta)
    eff_vals.append(eff_test)
            
    current_theta += step

# create graph of efficiency vs theta_in
n_points = len(theta_vals)
arr_theta = array('d', theta_vals)
arr_eff = array('d', eff_vals)

gr_eff = ROOT.TGraph(n_points, arr_theta, arr_eff)
gr_eff.SetTitle("Angular Acceptance Curve;Incoming Angle #theta_{in,x} [#murad];Channeling Efficiency [%]")
gr_eff.SetMarkerStyle(20); gr_eff.SetMarkerSize(0.6); gr_eff.SetMarkerColor(ROOT.kBlue+2)

fit_min = best_theta_0_raw - 5.0
fit_max = best_theta_0_raw + 20.0
gaus_eff = ROOT.TF1("gaus_eff", "gaus", fit_min, fit_max)

gaus_eff.SetParameters(max_eff_raw, best_theta_0_raw, 5.0)
gaus_eff.SetLineColor(ROOT.kRed)
gaus_eff.SetLineWidth(2)
gr_eff.Fit(gaus_eff, "RQ")

best_theta_0 = gaus_eff.GetParameter(1)
best_max_eff = gaus_eff.GetParameter(0)

print(f"\tMiglior theta_0 trovato: {best_theta_0:.2f} urad")
print(f"\tEfficienza Globale Massima: {(best_max_eff):.2f}%")

c3 = ROOT.TCanvas("c3", "Angular Acceptance Curve", 1000, 700)
gr_eff.Draw("AP")
gaus_eff.Draw("SAME")
leg = ROOT.TPaveText(0.55, 0.75, 0.90, 0.90, "NDC")
leg.SetBorderSize(1)
leg.SetFillColor(ROOT.kWhite)
leg.AddText(f"Max Efficiency #epsilon_{{ch}} = {best_max_eff:.1f} %"); leg.Draw("SAME")
c3.Update()

cut = f"abs(thetaIn_x - ({best_theta_0})) <= {theta_L / 2.0}"
df_phys = df_phys.Filter(cut)
P4 = df_phys.Count().GetValue()
print(f"Filter 4 (Lindhard cut): {((P3 - P4)/P3*100):.2f}% events out of {P3} got discarded.")
print("="*50)



# N TOT
N_tot = df_phys.Count().GetValue()
h_theta_cut = df_phys.Histo1D(("h_theta_cut", "#theta_{in,x} after filtering ; #theta_{x} [#murad]; Counts", 500, -150, 150), "thetaIn_x")
# bin width circa 20 urad
h_defl_cut = df_phys.Histo1D(("h_defl_cut", "Angular Deflection cut at #pm #theta_{L}/2; #Delta#theta_{x} [#murad]; No. particles", 5000, -2000, max_value), "Deltatheta_x")
h_cut_value = h_defl_cut.GetValue()

c4 = ROOT.TCanvas("c4", "Channeling Efficiency Analysis", 1400, 900)
pad_center = ROOT.TPad("pad_center", "pad_center", 0, 0, 0.65, 1.0)
pad_center.Draw()
pad_right = ROOT.TPad("pad_right", "pad_right", 0.60, 0.0, 1.0, 1.0)
pad_right.Draw()
pad_center.cd(); h_scan.Draw("COL")
line1 = ROOT.TLine(best_theta_0 - theta_L/2.0, -2000, best_theta_0 - theta_L/2.0, max_value)
line2 = ROOT.TLine(best_theta_0 + theta_L/2.0, -2000, best_theta_0 + theta_L/2.0, max_value)
line1.SetLineColor(ROOT.kOrange+10); line1.SetLineStyle(2); line1.SetLineWidth(2); line1.Draw("SAME")
line2.SetLineColor(ROOT.kOrange+10); line2.SetLineStyle(2); line2.SetLineWidth(2); line2.Draw("SAME")
pad_right.cd(); h_cut_value.SetFillColor(ROOT.kOrange+10); h_cut_value.GetXaxis().SetTitle(""); h_cut_value.Draw("HBAR Y+")
c4.Update()


# GAUSSIAN FIT
eff_ch = 0.0
eff_err = 0.0
N_ch = 0.0
fit_mean = 0.0
fit_sigma = 0.0
    
if N_tot > 0:
    # fitting only right side of the peak (cleanest one)
    fit_min = deflection_peak - 15
    fit_max = deflection_peak + 500
    # fit_max = max_value
    gaus_fit = ROOT.TF1("gaus_fit", "gaus", fit_min, fit_max)
    gaus_fit.SetLineColor(ROOT.kRed)
    
    # R = usa il range specificato, Q = modalità silenziosa, 0 = non disegnare subito
    h_cut_value.Fit(gaus_fit, "RQ0")

    fit_mean = gaus_fit.GetParameter(1)
    fit_sigma = gaus_fit.GetParameter(2)
    
    # COMPUTE N CH FROM GAUSSIAN INTEGRAL
    bin_width = h_cut_value.GetBinWidth(1)
    # method 1: integral of the gaussian fit
    # N_ch = gaus_fit.Integral(5950, max_value) / bin_width
    # method 2: integral of the histogram
    bin_min = h_cut_value.FindBin(fit_mean - 3.0 * fit_sigma)
    # bin_max = h_cut_value.FindBin(fit_mean + 4.0 * fit_sigma)
    bin_max = h_cut_value.FindBin(max_value)
    N_ch = h_cut_value.Integral(bin_min, bin_max)
    # print(f"Number of channeled particles (N_ch) = {N_ch:.0f} (from bin {bin_min} to {bin_max})")
    
    # EFFICIENCY WITH ITS ERROR (binomial distribution)
    eff_ch = (N_ch / N_tot) * 100.0
    eff_err = math.sqrt(eff_ch/100.0 * (1.0 - eff_ch/100.0) / N_tot) * 100.0


c5 = ROOT.TCanvas("c5", "Channeling Efficiency Fit", 1400, 900)
h_cut_value.SetFillColorAlpha(ROOT.kOrange, 0.6)
h_cut_value.SetLineColor(ROOT.kOrange)
h_cut_value.Draw("HIST")
c5.SetLogy()
c5.Update()

c6 = ROOT.TCanvas("c6", "Channeling Efficiency Fit", 1400, 900)
h_cut_value.SetFillColorAlpha(ROOT.kOrange, 0.6)
h_cut_value.SetLineColor(ROOT.kOrange)
h_cut_value.GetXaxis().SetRangeUser(5900, 6150)
h_cut_value.Draw("HIST")
if N_tot > 0:
    gaus_fit.Draw("SAME")

# line3 = ROOT.TLine(fit_mean - 3*fit_sigma, 0, fit_mean - 3*fit_sigma, h_cut_value.GetMaximum())
# line3.SetLineColor(ROOT.kRed); line3.SetLineStyle(2); line3.SetLineWidth(2); line3.Draw()
# line4 = ROOT.TLine(fit_mean + 4*fit_sigma, 0, fit_mean + 4*fit_sigma, h_cut_value.GetMaximum())
# line4.SetLineColor(ROOT.kRed); line4.SetLineStyle(2); line4.SetLineWidth(2); line4.Draw()

legend = ROOT.TPaveText(0.75, 0.75, 0.90, 0.90, "NDC")
legend.SetBorderSize(1)
legend.SetFillColor(ROOT.kWhite)
legend.SetTextAlign(12)
legend.AddText(f"Cut at #pm #theta_{{L}}/2 (#theta_{{0}} = {best_theta_0:.2f} #murad)")
legend.AddText(f"#epsilon_{{ch}} = {eff_ch:.1f} #pm {eff_err:.1f} %")
legend.AddText(f"Fit mean: {fit_mean:.1f} #murad")
legend.Draw()
c6.Update()



# c0.SaveAs(f"plots_{file}/spatial_cut_box.pdf")
# c0.SaveAs(f"plots_{file}/spatial_cut_box.png")

# c0b.SaveAs(f"plots_{file}/spatial_cut2_box.pdf")
# c0b.SaveAs(f"plots_{file}/spatial_cut2_box.png")

# c1.SaveAs(f"plots_{file}/d0_Out_projxy_box.pdf")
# c1.SaveAs(f"plots_{file}/d0_Out_projxy_box.png")

# c2.SaveAs(f"plots_{file}/d0_projxy_box.pdf")
# c2.SaveAs(f"plots_{file}/d0_projxy_box.png")

# c3.SaveAs(f"plots_{file}/angular_acceptance_curve.pdf")
# c3.SaveAs(f"plots_{file}/angular_acceptance_curve.png")

# c4.SaveAs(f"plots_{file}/channeling_efficiency_scan.pdf")
# c4.SaveAs(f"plots_{file}/channeling_efficiency_scan.png")

# c5.SaveAs(f"plots_{file}/channeling_efficiency_fit.pdf")
# c5.SaveAs(f"plots_{file}/channeling_efficiency_fit.png")

# c6.SaveAs(f"plots_{file}/channeling_efficiency_fit_zoom.pdf")
# c6.SaveAs(f"plots_{file}/channeling_efficiency_fit_zoom.png")




