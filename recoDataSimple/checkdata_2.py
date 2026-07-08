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
theta_0 = 0 # centro della tua isola di channeling - looking at 2D "h_scan" find on which x (thetaIn) the channeling isle is centered
max_value = 0 # histograms range
deflection_peak = 6010.0 # urad

if file in [8430, 8431, 8650]:
    theta_L = 12.9
    width = 8
    height = 2
    max_value = 8000
elif file in [8655, 8656]: 
    theta_L = 12.5
    width = 22.5
    height = 2
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

# FILTER 3: spatial cut (d0_x and d0_y within the crystal area)
# i look at the beam profile of d0 of those who channeled, and find the region width x length
# ---------------------------------------------------------------------------------------------
# d0_x = df_phys.Histo1D(("d0_x", "d0_x; d0_x [mm]; Counts", 100, -20, 30), "Tracks.d0_x")
# d0_y = df_phys.Histo1D(("d0_y", "d0_y; d0_y [mm]; Counts", 100, -20, 30), "Tracks.d0_y")
# x_c = d0_x.GetMean() # centro del cristallo in x (mm)
# y_c = d0_y.GetMean() # centro del cristallo in y (mm)
# x_min, x_max = x_c - width / 2.0, x_c + width / 2.0
# y_min, y_max = y_c - height / 2.0, y_c + height / 2.0
# spatial_cut = f"Tracks.d0_x > {x_min} && Tracks.d0_x < {x_max} && Tracks.d0_y > {y_min} && Tracks.d0_y < {y_max}"
# df_phys = df_phys.Filter(spatial_cut, "Spatial Cut (Crystal Area)")
# print(f"{((P0 - P2)/P0*100):.2f}% events out of {P0} got discarded (SingleTrack = 0, Spatial Cut).")

# HISTOGRAMS PRE CUT
h_theta_before = df_phys.Histo1D(("h_theta_before", "#theta_{in,x}; #theta_{x} [#murad]; Counts", 500, -150, 150), "thetaIn_x")
h_defl_before = df_phys.Histo1D(("h_defl_before", "Angular Deflection; #Delta#theta_{x} [#murad]; Counts", 500, -2000, max_value), "Deltatheta_x") 
h_scan = df_phys.Histo2D(("h_scan", "Angular deflection of the particles as a function of the incident angle; Incident angle #theta_{In, x} [#murad]; Deflection #Delta#theta_{x} [#murad]", 500, -150, 150, 500, -2000, max_value), "thetaIn_x", "Deltatheta_x")

# FINDING THETA 0
df_theta0a = df_phys.Filter("Deltatheta_x > 5500")
h_theta0a = df_theta0a.Histo1D(("h_theta0a", "#theta_{in,x}; #theta_{x} [#murad]; Counts", 500, -150, 150), "thetaIn_x")
df_theta0b = df_phys.Filter("Deltatheta_x < 5000")
h_theta0b = df_theta0b.Histo1D(("h_theta0b", "#theta_{in,x}; #theta_{x} [#murad]; Counts", 500, -150, 150), "thetaIn_x")

# SPATIAL CUT
h_d0_x_ch = df_theta0a.Histo1D(("h_d0_x_ch", "d0_x of channeled particles; d0_x [mm]; Counts", 100, -20, 20), "Tracks.d0_x")
h_d0_y_ch = df_theta0a.Histo1D(("h_d0_y_ch", "d0_y of channeled particles; d0_y [mm]; Counts", 100, -20, 20), "Tracks.d0_y")
h_d0_xy_ch = df_theta0a.Histo2D(("h_d0_xy_ch", "d0_x vs d0_y of channeled particles; d0_x [mm]; d0_y [mm]", 200, -5, 5, 200, -10, 10), "Tracks.d0_x", "Tracks.d0_y")
h_d0_xy = df_phys.Histo2D(("h_d0_xy", "d0_x vs d0_y of all particles; d0_x [mm]; d0_y [mm]", 200, -5, 5, 200, -10, 10), "Tracks.d0_x", "Tracks.d0_y")

# LINDHARD CUT
# only keeping events whose theta_in is within +/- 0.5 * theta_L
theta_0 = h_theta0a.GetMean()
cut = f"abs(thetaIn_x - ({theta_0})) <= {theta_L / 2.0}"
df_cut = df_phys.Filter(cut)
    
# how many particles survive the cut (N_tot)
N_tot = df_cut.Count().GetValue()
print(f"{(N_tot/df_phys.Count().GetValue()*100):.2f}% events selected after the +/- 1/2 theta_L cut ({(N_tot/df.Count().GetValue()*100):.2f}% of the total).")
print("="*50)

h_theta_cut = df_cut.Histo1D(("h_theta_cut", "#theta_{in,x} after filtering ; #theta_{x} [#murad]; Counts", 500, -150, 150), "thetaIn_x")
# bin width circa 20 urad
h_defl_cut = df_cut.Histo1D(("h_defl_cut", "Angular Deflection cut at #pm #theta_{L}/2; #Delta#theta_{x} [#murad]; No. particles", 5000, -2000, max_value), "Deltatheta_x")
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
c0 = ROOT.TCanvas("c0", "theta_0", 1200, 600)
c0.Divide(1,1)
# c0.cd(1)
# h_theta0b.Draw()
# h_theta0a.Draw("SAME")
# h_theta0a.SetLineColor(ROOT.kRed)

c0.cd(1)
h_theta_before.Draw()
h_theta_cut.SetLineColor(ROOT.kBlack)
h_theta0b.Draw("SAME"); h_theta0b.SetLineColor(ROOT.kGreen)
h_theta0a.Draw("SAME"); h_theta0a.SetLineColor(ROOT.kBlue)
h_theta_cut.Draw("SAME"); h_theta_cut.SetLineColor(ROOT.kRed)

c0.Update()

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

c2 = ROOT.TCanvas("c2", "d0_x and d0_y of channeled particles", 1800, 900)
c2.Divide(2, 2)
c2.cd(1)
h_d0_x_ch.Draw()
c2.cd(2)
h_d0_y_ch.Draw()    
c2.cd(3)
h_d0_xy_ch.Draw("COLZ")
c2.cd(4)   
h_d0_xy.Draw("COLZ")