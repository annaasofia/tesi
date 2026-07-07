import ROOT
import math
import sys

ROOT.ROOT.EnableImplicitMT() 

# file = input("File number: ")
file = 8430
filename = "recoDataSimple_" + str(file) + "_xtalMerging.root"

ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetPalette(ROOT.kBird)

# VARIABLES
theta_L = 0 # Lindhard angle for TCCP and TCCPA in urad
theta_0 = 0 # centro della tua isola di channeling - looking at 2D "h_scan" find on which x (thetaIn) the channeling isle is centered
max_value = 0 # histograms range
deflection_peak = 6010.0 # urad

if file in [8430, 8431, 8650]:
    theta_L = 12.9
    max_value = 8000
elif file in [8655, 8656]: 
    theta_L = 12.5
    max_value = 14000
else:
    print("Run not found...")
    sys.exit(1)

# ROOT DATA FRAME
df = ROOT.RDataFrame("simpleEvent", filename)
print("="*50)
print(f"Analyzing {filename} ...")

# FILTERING the data and conversion from rad to urad
df_clean = df.Filter("SingleTrack == 1")
print("="*50)
print(f"{((df.Count().GetValue() - df_clean.Count().GetValue())/df.Count().GetValue()*100):.2f}% events out of {df.Count().GetValue()} got discarded (SingleTrack = 0).")
df_phys = df_clean.Define("thetaIn_x", "Tracks.thetaIn_x * 1e6").Define("Deltatheta_x", "(Tracks.thetaOut_x - Tracks.thetaIn_x) * 1e6")

# HISTOGRAMS PRE CUT
h_theta_before = df_phys.Histo1D(("h_theta_before", "#theta_{in,x}; #theta_{x} [#murad]; Counts", 200, -150, 150), "thetaIn_x")
h_defl_before = df_phys.Histo1D(("h_defl_before", "Angular Deflection; #Delta#theta_{x} [#murad]; Counts", 200, -2000, max_value), "Deltatheta_x") 
h_scan = df_phys.Histo2D(("h_scan", "Angular deflection of the particles as a function of the incident angle; Incident angle #theta_{In, x} [#murad]; Deflection #Delta#theta_{x} [#murad]", 500, -150, 150, 500, -2000, max_value), "thetaIn_x", "Deltatheta_x")

# FINDING THETA 0
df_theta0a = df_phys.Filter("Deltatheta_x > 5000")
h_theta0a = df_theta0a.Histo1D(("h_theta0a", "#theta_{in,x}; #theta_{x} [#murad]; Counts", 200, -150, 150), "thetaIn_x")
df_theta0b = df_phys.Filter("Deltatheta_x < 5000")
h_theta0b = df_theta0b.Histo1D(("h_theta0b", "#theta_{in,x}; #theta_{x} [#murad]; Counts", 200, -150, 150), "thetaIn_x")

# LINDHARD CUT
# only keeping events whose theta_in is within +/- 0.5 * theta_L
theta_0 = h_theta0a.GetMean()
cut = f"abs(thetaIn_x - ({theta_0})) <= {theta_L / 2.0}"
df_cut = df_phys.Filter(cut)
    
# how many particles survive the cut (N_tot)
N_tot = df_cut.Count().GetValue()
print(f"{(N_tot/df_clean.Count().GetValue()*100):.2f}% events selected after the +/- 1/2 theta_L cut ({(N_tot/df.Count().GetValue()*100):.2f}% of the total).")
print("="*50)

# bin width circa 50 urad
h_theta_cut = df_cut.Histo1D(("h_theta_cut", "#theta_{in,x} after filtering ; #theta_{x} [#murad]; Counts", 200, -150, 150), "thetaIn_x")
h_defl_cut = df_cut.Histo1D(("h_defl_cut", "Angular Deflection cut at #pm #theta_{L}/2; #Delta#theta_{x} [#murad]; No. particles", 500, -2000, max_value), "Deltatheta_x")
h_cut_value = h_defl_cut.GetValue()

# GAUSSIAN FIT
eff_ch = 0.0
eff_err = 0.0
N_ch = 0.0
fit_mean = 0.0
fit_sigma = 0.0
    
if N_tot > 0:
    # fitting only right side of the peak (cleanest one)
    fit_min = deflection_peak - 500
    # fit_max = deflection_peak + 500
    fit_max = max_value
    gaus_fit = ROOT.TF1("gaus_fit", "gaus", fit_min, fit_max)
    gaus_fit.SetLineColor(ROOT.kRed)
    
    # R = usa il range specificato, Q = modalità silenziosa, 0 = non disegnare subito
    h_cut_value.Fit(gaus_fit, "RQ0")
    
    # COMPUTE N CH FROM GAUSSIAN INTEGRAL
    bin_width = h_cut_value.GetBinWidth(1)
    N_ch = gaus_fit.Integral(fit_min, fit_max) / bin_width
    
    fit_mean = gaus_fit.GetParameter(1)
    fit_sigma = gaus_fit.GetParameter(2)
    
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

legend = ROOT.TPaveText(0.45, 0.75, 0.90, 0.90, "NDC")
legend.SetBorderSize(1)
legend.SetFillColor(ROOT.kWhite)
legend.SetTextAlign(12)
legend.AddText(f"Cut at #pm #theta_{{L}}/2 (#theta_{{0}} = {theta_0:.2f} #murad)")
legend.AddText(f"#epsilon_{{ch}} = {eff_ch:.1f} #pm {eff_err:.1f} %")
legend.AddText(f"Fit mean: {fit_mean:.1f} #murad")
legend.Draw()

c1.Update()