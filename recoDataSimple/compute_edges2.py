import ROOT
import math
import sys
import numpy as np
from array import array



ROOT.ROOT.EnableImplicitMT() 
ROOT.gStyle.SetOptStat(0)

file = 8430
filename = "recoDataSimple_" + str(file) + "_xtalMerging.root"
files = ["recoDataSimple_8430_xtalMerging.root", "recoDataSimple_8431_xtalMerging.root"]

ROOT.gStyle.SetPalette(ROOT.kBird)

# VARIABLES
width = 12.8 # width of the crystal in mm (for spatial cut)
height = 2 # height of the crystal in mm (for spatial cut)
max_value = 8000 # histograms range
deflection_peak = 6010 # urad
fit_range = 150 # urad
minimum_entries = 30

# ROOT DATA FRAME
df = ROOT.RDataFrame("simpleEvent", filename)
# df = ROOT.RDataFrame("simpleEvent", files)
print("="*50)
print(f"Analyzing {filename} ...")

# FILTERING the data: single tracks and conversion from rad to urad
df_phys = df.Filter("SingleTrack == 1")
df_phys = df_phys.Define("thetaIn_x", "Tracks.thetaIn_x * 1e6")\
    .Define("Deltatheta_x", "(Tracks.thetaOut_x - Tracks.thetaIn_x) * 1e6")\
    .Define("Deltatheta_y", "(Tracks.thetaOut_y - Tracks.thetaIn_y) * 1e6")\
    .Define("DeltathetaErr_x", "sqrt(Tracks.thetaInErr_x * Tracks.thetaInErr_x + Tracks.thetaOutErr_x * Tracks.thetaOutErr_x) * 1e6")\
    .Define("DeltathetaErr_y", "sqrt(Tracks.thetaInErr_y * Tracks.thetaInErr_y + Tracks.thetaOutErr_y * Tracks.thetaOutErr_y) * 1e6")
print("="*50)

def book_scan_histogram(df_in, scan_var, delta_var, scan_min, scan_max, n_bins, slice_var=None, slice_min=None, slice_max=None, n_dtheta_bins=300):

    df_slice = df_in
    if slice_var is not None:
        df_slice = df_slice.Filter(f"{slice_var} > {slice_min} && {slice_var} < {slice_max}")
 
    h2 = df_slice.Histo2D((f"h2_scat_{scan_var}", "", n_bins, scan_min, scan_max, n_dtheta_bins, -fit_range, fit_range), scan_var, delta_var)
    return h2

def extract_widths_from_h2(h2_lazy, range_fit_gaus=fit_range, min_entries=minimum_entries):
    h2v = h2_lazy.GetValue()
    centers, sigmas, sigma_errs = [], [], []
 
    n_bins = h2v.GetNbinsX()
    for ix in range(1, n_bins + 1):
        h1 = h2v.ProjectionY(f"py_{h2v.GetName()}_{ix}", ix, ix)
 
        if h1.GetEntries() < min_entries:
            continue
 
        f = ROOT.TF1(f"f_{h2v.GetName()}_{ix}", "gaus", -range_fit_gaus, range_fit_gaus)
        f.SetParameters(h1.GetMaximum(), h1.GetMean(), max(h1.GetRMS(), 5.0))
        h1.Fit(f, "RQ0")
 
        centers.append(h2v.GetXaxis().GetBinCenter(ix))
        sigmas.append(f.GetParameter(2))
        sigma_errs.append(f.GetParError(2))
 
    return centers, sigmas, sigma_errs
 
def fit_step_edges(centers, sigmas, sigma_errs, edge_lo_guess, edge_hi_guess, baseline_guess=None, amplitude_guess=None, transition_guess=0.05):

    n = len(centers)
    if n == 0:
        print("ERROR: Zero data points! The previous cut filtered out all events.")
        return 0, 0, 0, 0, None, None

    gr = ROOT.TGraphErrors(n, array('d', centers), array('d', sigmas), array('d', [0.0] * n), array('d', sigma_errs))
 
    # 1. Rilevamento Intelligente: Guardiamo i dati per capire la forma
    mid_val = sigmas[n // 2] # Valore al centro
    edge_val = (sigmas[0] + sigmas[-1]) / 2.0 # Valore medio ai bordi esterni
 
    if baseline_guess is None:
        baseline_guess = edge_val
    if amplitude_guess is None:
        # Se è una buca (mid_val < edge_val), l'ampiezza diventa automaticamente negativa!
        amplitude_guess = mid_val - edge_val 
 
    f = ROOT.TF1("f_step_edges", "[0] + [1]*0.5*(TMath::Erf((x-[2])/[4]) - TMath::Erf((x-[3])/[4]))", min(centers), max(centers))
    f.SetParameters(baseline_guess, amplitude_guess, edge_lo_guess, edge_hi_guess, transition_guess)
    f.SetParLimits(4, 0.001, (max(centers) - min(centers)) / 4.0) 
 
    gr.Fit(f, "RQ0")
 
    edge_lo = f.GetParameter(2)
    edge_hi = f.GetParameter(3)
    edge_lo_err = f.GetParError(2)
    edge_hi_err = f.GetParError(3)
 
    # Se il fit si è "girato", rimettiamo i bordi nell'ordine giusto (minore, maggiore)
    if edge_lo > edge_hi:
        edge_lo, edge_hi = edge_hi, edge_lo
        edge_lo_err, edge_hi_err = edge_hi_err, edge_lo_err
 
    return edge_lo, edge_hi, edge_lo_err, edge_hi_err, f, gr

 
# y edges (scan y, using the full x range -- or restrict to rough x window from the old method / previous run)
h2_y = book_scan_histogram(df_phys, "Tracks.d0Out_y", "Deltatheta_y", scan_min=-5, scan_max=8, n_bins=200, slice_var="Tracks.d0_x", slice_min=-1.5, slice_max=1.5) # 0.01 mm precision
y_centers, y_sigmas, y_sigma_errs = extract_widths_from_h2(h2_y)
 
y_lo, y_hi, y_lo_err, y_hi_err, f_y, gr_y = fit_step_edges(y_centers, y_sigmas, y_sigma_errs, edge_lo_guess=-1.0, edge_hi_guess=11.0)
 
# x edges (scan x, restricted to the just-found y window)
h2_x = book_scan_histogram(df_phys, "Tracks.d0Out_x", "Deltatheta_x", scan_min=-3, scan_max=4, n_bins=140, slice_var="Tracks.d0_y", slice_min=y_lo, slice_max=y_hi)
x_centers, x_sigmas, x_sigma_errs = extract_widths_from_h2(h2_x)
 
x_lo, x_hi, x_lo_err, x_hi_err, f_x, gr_x = fit_step_edges(x_centers, x_sigmas, x_sigma_errs, edge_lo_guess=-1.0, edge_hi_guess=1.0, transition_guess=0.02)

h2_y = book_scan_histogram(df_phys, "Tracks.d0Out_y", "Deltatheta_y", scan_min=-5, scan_max=8, n_bins=200, slice_var="Tracks.d0_x", slice_min=x_lo, slice_max=x_hi)
y_centers, y_sigmas, y_sigma_errs = extract_widths_from_h2(h2_y)
y_lo, y_hi, y_lo_err, y_hi_err, f_y, gr_y = fit_step_edges(y_centers, y_sigmas, y_sigma_errs, edge_lo_guess=y_lo, edge_hi_guess=y_hi)

print(f"X edges (scattering method): [{x_lo:.4f} ± {x_lo_err:.4f}, {x_hi:.4f} ± {x_hi_err:.4f}] mm")
print(f"Y edges (scattering method): [{y_lo:.4f} ± {y_lo_err:.4f}, {y_hi:.4f} ± {y_hi_err:.4f}] mm")


mean_res_x = df_phys.Mean("DeltathetaErr_x").GetValue()
mean_res_y = df_phys.Mean("DeltathetaErr_y").GetValue()

sigma_bsl_x, sigma_jump_x = f_x.GetParameter(0), f_x.GetParameter(1)
sigma_bsl_y, sigma_jump_y = f_y.GetParameter(0), f_y.GetParameter(1)

sigma_mcs_out_x = math.sqrt(max(sigma_bsl_x**2 - mean_res_x**2, 0))
sigma_mcs_in_x  = math.sqrt(max((sigma_bsl_x + sigma_jump_x)**2 - mean_res_x**2, 0))
sigma_mcs_out_y = math.sqrt(max(sigma_bsl_y**2 - mean_res_y**2, 0))
sigma_mcs_in_y  = math.sqrt(max((sigma_bsl_y + sigma_jump_y)**2 - mean_res_y**2, 0))

print(f"theta_mcs (x): baseline {sigma_mcs_out_x:.2f} -> in-crystal {sigma_mcs_in_x:.2f} urad")
print(f"theta_mcs (y): baseline {sigma_mcs_out_y:.2f} -> in-crystal {sigma_mcs_in_y:.2f} urad")

mean_d0err_x = df_phys.Mean("Tracks.d0Err_x").GetValue()
mean_d0err_y = df_phys.Mean("Tracks.d0Err_y").GetValue()
print(f"fitted transition x = {f_x.GetParameter(4):.4f} mm vs mean d0Err_x = {mean_d0err_x:.4f} mm")
print(f"fitted transition y = {f_y.GetParameter(4):.4f} mm vs mean d0Err_y = {mean_d0err_y:.4f} mm")

c_y = ROOT.TCanvas("c_y", "Scattering width vs y", 900, 600)
gr_y.SetTitle("Local scattering width vs y; d0_y [mm]; #sigma(#Delta#theta_{y}) [#murad]")
gr_y.SetMarkerStyle(20); gr_y.SetMarkerSize(0.6)
gr_y.Draw("AP")
f_y.SetLineColor(ROOT.kRed)
f_y.Draw("SAME")
c_y.Update()

c_x = ROOT.TCanvas("c_x", "Scattering width vs x", 900, 600)
gr_x.SetTitle("Local scattering width vs x; d0_x [mm]; #sigma(#Delta#theta_{x}) [#murad]")
gr_x.SetMarkerStyle(20); gr_x.SetMarkerSize(0.6)
gr_x.Draw("AP")
f_x.SetLineColor(ROOT.kRed)
f_x.Draw("SAME")
c_x.Update()

 
print("=" * 50)
print(f"Crystal footprint from scattering method:")
print(f"\tx = [{x_lo:.4f}({x_lo_err*1e4:.0f}), {x_hi:.4f}({x_hi_err*1e4:.0f})] mm  (width = {x_hi - x_lo:.4f} ± {pow(x_lo_err*x_lo_err + x_hi_err*x_hi_err,0.5):.4f} mm)")
print(f"\ty = [{y_lo:.4f}({y_lo_err*1e4:.0f}), {y_hi:.4f}({y_hi_err*1e4:.0f})] mm  (width = {y_hi - y_lo:.4f} ± {pow(y_lo_err*y_lo_err + y_hi_err*y_hi_err,0.5):.4f} mm)")
print(f"\tbaseline sigma x  = {f_x.GetParameter(0):.2f} ± {f_x.GetParError(0):.2f} urad")
print(f"\tin-crystal jump x = {f_x.GetParameter(1):.2f} ± {f_x.GetParError(1):.2f} urad")
print(f"\tbaseline sigma y  = {f_y.GetParameter(0):.2f} ± {f_y.GetParError(0):.2f} urad")
print(f"\tin-crystal jump y = {f_y.GetParameter(1):.2f} ± {f_y.GetParError(1):.2f} urad")

# Plot della mappa 2D Posizione Y vs Deflessione (h2_scat_y)
c_2d_y = ROOT.TCanvas("c_2d_y", "Position Y vs Deflection", 900, 600)
h2_y_val = h2_y.GetValue()
h2_y_val.SetTitle("Posizione Y vs Deflessione #Delta#theta_{y}; d0_y [mm]; #Delta#theta_{x} [#murad]")
h2_y_val.Draw("colz")
c_2d_y.Update()

# Plot della mappa 2D Posizione X vs Deflessione (h2_scat_x)
c_2d_x = ROOT.TCanvas("c_2d_x", "Position X vs Deflection", 900, 600)
h2_x_val = h2_x.GetValue()
h2_x_val.SetTitle("Posizione X vs Deflessione #Delta#theta_{x}; d0_x [mm]; #Delta#theta_{x} [#murad]")
h2_x_val.Draw("colz")
c_2d_x.Update()

# Estrazione e plot di DUE singole fettine (Una DENTRO e una FUORI dal cristallo)
# Scegliamo due bin a caso basandoci sui centri Y trovati
bin_dentro = h2_y_val.GetXaxis().FindBin((y_lo + y_hi) / 2.0) # Esattamente a metà cristallo
bin_fuori = h2_y_val.GetXaxis().FindBin(y_lo - 2.0)           # 2 mm fuori dal bordo inferiore

# Estraiamo gli istogrammi 1D
h1_dentro = h2_y_val.ProjectionY("h1_dentro", bin_dentro, bin_dentro)
h1_fuori = h2_y_val.ProjectionY("h1_fuori", bin_fuori, bin_fuori)

# Fit e Plot per la zona DENTRO il cristallo
c_slice_in = ROOT.TCanvas("c_slice_in", "Fit Slice (INSIDE Crystal)", 800, 600)
h1_dentro.SetTitle(f"Fettina DENTRO il cristallo (Y = {h2_y_val.GetXaxis().GetBinCenter(bin_dentro):.2f} mm); #Delta#theta_{{x}} [#murad]; Counts")
h1_dentro.SetLineColor(ROOT.kBlue); h1_dentro.SetLineWidth(2)
f_in = ROOT.TF1("f_in", "gaus", -fit_range, fit_range)
f_in.SetLineColor(ROOT.kRed)
h1_dentro.Fit(f_in, "RQ")
h1_dentro.Draw("HIST")
f_in.Draw("SAME")

leg_in = ROOT.TLegend(0.65, 0.75, 0.88, 0.88); leg_in.SetBorderSize(0)
leg_in.AddEntry(f_in, f"#sigma = {f_in.GetParameter(2):.1f} #murad", "l")
leg_in.Draw("SAME")
c_slice_in.Update()

# Fit e Plot per la zona FUORI dal cristallo
c_slice_out = ROOT.TCanvas("c_slice_out", "Fit Slice (OUTSIDE Crystal)", 800, 600)
h1_fuori.SetTitle(f"Fettina FUORI dal cristallo (Y = {h2_y_val.GetXaxis().GetBinCenter(bin_fuori):.2f} mm); #Delta#theta_{{x}} [#murad]; Counts")
h1_fuori.SetLineColor(ROOT.kMagenta); h1_fuori.SetLineWidth(2)
f_out = ROOT.TF1("f_out", "gaus", -fit_range, fit_range)
f_out.SetLineColor(ROOT.kRed)
h1_fuori.Fit(f_out, "RQ")
h1_fuori.Draw("HIST")
f_out.Draw("SAME")

leg_out = ROOT.TLegend(0.65, 0.75, 0.88, 0.88); leg_out.SetBorderSize(0)
leg_out.AddEntry(f_out, f"#sigma = {f_out.GetParameter(2):.1f} #murad", "l")
leg_out.Draw("SAME")
c_slice_out.Update()

# Scegliamo un bin al centro del cristallo e uno 1.5 mm fuori dal bordo destro
bin_dentro_x = h2_x_val.GetXaxis().FindBin((x_lo + x_hi) / 2.0) 
bin_fuori_x  = h2_x_val.GetXaxis().FindBin(x_hi + 1.5)           

# Estraiamo gli istogrammi 1D
h1_dentro_x = h2_x_val.ProjectionY("h1_dentro_x", bin_dentro_x, bin_dentro_x)
h1_fuori_x  = h2_x_val.ProjectionY("h1_fuori_x", bin_fuori_x, bin_fuori_x)

# Fit e Plot per la zona DENTRO il cristallo (Asse X)
c_slice_in_x = ROOT.TCanvas("c_slice_in_x", "Fit Slice X (INSIDE Crystal)", 800, 600)
h1_dentro_x.SetTitle(f"Fettina X DENTRO il cristallo (X = {h2_x_val.GetXaxis().GetBinCenter(bin_dentro_x):.2f} mm); #Delta#theta_{{x}} [#murad]; Counts")
h1_dentro_x.SetLineColor(ROOT.kBlue+2); h1_dentro_x.SetLineWidth(2)

f_in_x = ROOT.TF1("f_in_x", "gaus", -fit_range, fit_range)
f_in_x.SetLineColor(ROOT.kRed)
h1_dentro_x.Fit(f_in_x, "RQ")
h1_dentro_x.Draw("HIST")
f_in_x.Draw("SAME")

leg_in_x = ROOT.TLegend(0.65, 0.75, 0.88, 0.88); leg_in_x.SetBorderSize(0)
leg_in_x.AddEntry(f_in_x, f"#sigma = {f_in_x.GetParameter(2):.1f} #murad", "l")
leg_in_x.Draw("SAME")
c_slice_in_x.Update()

# Fit e Plot per la zona FUORI dal cristallo (Asse X)
c_slice_out_x = ROOT.TCanvas("c_slice_out_x", "Fit Slice X (OUTSIDE Crystal)", 800, 600)
h1_fuori_x.SetTitle(f"Fettina X FUORI dal cristallo (X = {h2_x_val.GetXaxis().GetBinCenter(bin_fuori_x):.2f} mm); #Delta#theta_{{x}} [#murad]; Counts")
h1_fuori_x.SetLineColor(ROOT.kMagenta+2); h1_fuori_x.SetLineWidth(2)

f_out_x = ROOT.TF1("f_out_x", "gaus", -fit_range, fit_range)
f_out_x.SetLineColor(ROOT.kRed)
h1_fuori_x.Fit(f_out_x, "RQ")
h1_fuori_x.Draw("HIST")
f_out_x.Draw("SAME")

leg_out_x = ROOT.TLegend(0.65, 0.75, 0.88, 0.88); leg_out_x.SetBorderSize(0)
leg_out_x.AddEntry(f_out_x, f"#sigma = {f_out_x.GetParameter(2):.1f} #murad", "l")
leg_out_x.Draw("SAME")
c_slice_out_x.Update()

for obj in [c_2d_y, c_2d_x, c_slice_in, c_slice_out, h1_dentro, h1_fuori, f_in, f_out, leg_in, leg_out,
            c_slice_in_x, c_slice_out_x, h1_dentro_x, h1_fuori_x, f_in_x, f_out_x, leg_in_x, leg_out_x]:
    ROOT.SetOwnership(obj, False)


 
