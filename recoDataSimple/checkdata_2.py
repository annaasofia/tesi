import ROOT
import math
import sys

ROOT.ROOT.EnableImplicitMT() 

# file = input("File number: ")
file = 8430
filename = "recoDataSimple_" + str(file) + "_xtalMerging.root"

ROOT.gStyle.SetOptStat(1111)
ROOT.gStyle.SetPalette(ROOT.kBird)

# VARIABLES
theta_L = 0 # Lindhard angle for TCCP and TCCPA in urad
width = 0 # width of the crystal in mm (for spatial cut)
height = 0 # height of the crystal in mm (for spatial cut)
delta_x = 0 # shift in x (for spatial cut)
theta_0 = 0 # centro della tua isola di channeling - looking at 2D "h_scan" find on which x (thetaIn) the channeling isle is centered
max_value = 0 # histograms range
deflection_peak = 0 # urad

if file in [8430, 8431, 8650]:
    theta_L = 12.9
    deflection_peak = 6010.0
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
h_d0_xy_ch = df_cut.Histo2D(("h_d0_xy_ch", "d0_x vs d0_y of channeled particles; d0_x [mm]; d0_y [mm]", 300, -2, 2, 160, -4, 6), "Tracks.d0_x", "Tracks.d0_y")
h_d0_Out_xy_ch = df_cut.Histo2D(("h_d0_Out_xy_ch", "d0Out_x vs d0Out_y of channeled particles; d0Out_x [mm]; d0Out_y [mm]", 300, -2, 2, 160, -4, 6), "Tracks.d0Out_x", "Tracks.d0Out_y")

h_d0_x_ch = df_cut.Histo1D(("h_d0_x_ch", "d0_x of channeled particles; d0_x [mm]; Counts", 300, -1.5, 1.5), "Tracks.d0_x")
h_d0_y_ch = df_cut.Histo1D(("h_d0_y_ch", "d0_y of channeled particles; d0_y [mm]; Counts", 100, -4, 6), "Tracks.d0_y")
h_d0_Out_x_ch = df_cut.Histo1D(("h_d0_Out_x_ch", "d0Out_x of channeled particles; d0Out_x [mm]; Counts", 300, -1.5, 1.5), "Tracks.d0Out_x")
h_d0_Out_y_ch = df_cut.Histo1D(("h_d0_Out_y_ch", "d0Out_y of channeled particles; d0Out_y [mm]; Counts", 100, -4, 6), "Tracks.d0Out_y")

y_c1 = h_d0_y_ch.GetMean(); y_c2 = h_d0_Out_y_ch.GetMean(); y_c = (y_c1 + y_c2) / 2.0 # centro del cristallo in y (mm)
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

# now that we have the crystal area/position, we can apply the spatial cut to the entire dataframe
spatial_cut = f"Tracks.d0_x > {x_min} && Tracks.d0_x < {x_max} && Tracks.d0_y > {y_min} && Tracks.d0_y < {y_max}"
df_phys = df_phys.Filter(spatial_cut, "Spatial Cut (Crystal Area)")
P3 = df_phys.Count().GetValue()
print(f"Filter 3 (Spatial cut): {((P2 - P3)/P2*100):.2f}% events out of {P2} got discarded.")

# HISTOGRAMS PRE ANGLE CUT
h_theta_before = df_phys.Histo1D(("h_theta_before", "#theta_{in,x}; #theta_{x} [#murad]; Counts", 500, -150, 150), "thetaIn_x")
h_defl_before = df_phys.Histo1D(("h_defl_before", "Angular Deflection; #Delta#theta_{x} [#murad]; Counts", 500, -2000, max_value), "Deltatheta_x") 
h_scan = df_phys.Histo2D(("h_scan", "Angular deflection of the particles as a function of the incident angle; Incident angle #theta_{In, x} [#murad]; Deflection #Delta#theta_{x} [#murad]", 500, -150, 150, 500, -2000, max_value), "thetaIn_x", "Deltatheta_x")

# FILTER 4: Lindhard cut (only keeping events that have entered the crystal with an angle within +/- 1/2 theta_L)
# FINDING THETA 0
h_theta0 = df_cut.Histo1D(("h_theta0", "#theta_{in,x}; #theta_{x} [#murad]; Counts", 500, -150, 150), "thetaIn_x")
theta_0 = h_theta0.GetMean()

cut = f"abs(thetaIn_x - ({theta_0})) <= {theta_L / 2.0}"
df_phys = df_phys.Filter(cut)
P4 = df_phys.Count().GetValue()
print(f"Filter 4 (Lindhard cut): {((P3 - P4)/P3*100):.2f}% events out of {P3} got discarded.")
    
# N TOT
N_tot = df_phys.Count().GetValue() # so actually is P4 
print("="*50)

h_theta_cut = df_phys.Histo1D(("h_theta_cut", "#theta_{in,x} after filtering ; #theta_{x} [#murad]; Counts", 500, -150, 150), "thetaIn_x")
# bin width circa 20 urad
h_defl_cut = df_phys.Histo1D(("h_defl_cut", "Angular Deflection cut at #pm #theta_{L}/2; #Delta#theta_{x} [#murad]; No. particles", 5000, -2000, max_value), "Deltatheta_x")
h_cut_value = h_defl_cut.GetValue()

# GAUSSIAN FIT
eff_ch = 0.0
eff_err = 0.0
N_ch = 0.0
fit_mean = 0.0
fit_sigma = 0.0
    
if N_tot > 0:
    # fitting only right side of the peak (cleanest one)
    fit_min = deflection_peak - 25
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
    print(f"Number of channeled particles (N_ch) = {N_ch:.0f} (from bin {bin_min} to {bin_max})")
    
    # EFFICIENCY WITH ITS ERROR (binomial distribution)
    eff_ch = (N_ch / N_tot) * 100.0
    eff_err = math.sqrt(eff_ch/100.0 * (1.0 - eff_ch/100.0) / N_tot) * 100.0

# DRAWINGS

c1 = ROOT.TCanvas("c1", "Channeling Efficiency Analysis", 1800, 700)
c1.Divide(3, 1) 

# Pad Angular deflection (Log scale)
c1.cd(1)
ROOT.gPad.SetLogy()
# h_defl_before.SetLineWidth(2)
h_defl_before.Draw("HIST")
# h_defl_cut.SetLineWidth(2)
# h_defl_cut.SetLineColor(ROOT.kRed)
# h_defl_cut.Draw("HIST SAME")

# Pad Scan Angolare 2D
c1.cd(2)
h_scan.Draw("COLZ TICKXY")
    
# Disegniamo due linee verticali per mostrare dove stiamo tagliando
line1 = ROOT.TLine(theta_0 - theta_L/2.0, -2000, theta_0 - theta_L/2.0, max_value)
line2 = ROOT.TLine(theta_0 + theta_L/2.0, -2000, theta_0 + theta_L/2.0, max_value)
line1.SetLineColor(ROOT.kRed); line1.SetLineStyle(2); line1.SetLineWidth(2); line1.Draw()
line2.SetLineColor(ROOT.kRed); line2.SetLineStyle(2); line2.SetLineWidth(2); line2.Draw()

# Pad Gaussian Fit after cut
c1.cd(3)
h_cut_value.SetFillColorAlpha(ROOT.kGreen+1, 0.6)
h_cut_value.SetLineColor(ROOT.kGreen+2)
h_cut_value.Draw("HIST")
if N_tot > 0:
    gaus_fit.Draw("SAME")
line3 = ROOT.TLine(fit_mean - 3*fit_sigma, 0, fit_mean - 3*fit_sigma, h_cut_value.GetMaximum())
line3.SetLineColor(ROOT.kRed); line3.SetLineStyle(2); line3.SetLineWidth(2); line3.Draw()
line4 = ROOT.TLine(fit_mean + 4*fit_sigma, 0, fit_mean + 4*fit_sigma, h_cut_value.GetMaximum())
line4.SetLineColor(ROOT.kRed); line4.SetLineStyle(2); line4.SetLineWidth(2); line4.Draw()

legend = ROOT.TPaveText(0.45, 0.75, 0.90, 0.90, "NDC")
legend.SetBorderSize(1)
legend.SetFillColor(ROOT.kWhite)
legend.SetTextAlign(12)
legend.AddText(f"Cut at #pm #theta_{{L}}/2 (#theta_{{0}} = {theta_0:.2f} #murad)")
legend.AddText(f"#epsilon_{{ch}} = {eff_ch:.1f} #pm {eff_err:.1f} %")
legend.AddText(f"Fit mean: {fit_mean:.1f} #murad")
legend.Draw()

print(f"#sigma = {fit_sigma:.2f}, so bin width should be = {3.49*fit_sigma/math.pow(N_tot,1/3)}")

c1.Update()