import ROOT
import math
import sys
import os
import numpy as np
from array import array
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator

def root_graph_to_arrays(gr):
    n = gr.GetN()
    x = np.array([gr.GetPointX(i) for i in range(n)])
    y = np.array([gr.GetPointY(i) for i in range(n)])
    ex = np.array([gr.GetErrorX(i) for i in range(n)])
    ey = np.array([gr.GetErrorY(i) for i in range(n)])
    return x, y, ex, ey

def root_tf1_to_curve(f, xmin, xmax, n=800):
    xs = np.linspace(xmin, xmax, n)
    ys = np.array([f.Eval(x) for x in xs])
    return xs, ys

def root_hist_to_arrays(h):
    n = h.GetNbinsX()
    centers = np.array([h.GetBinCenter(i) for i in range(1, n + 1)])
    contents = np.array([h.GetBinContent(i) for i in range(1, n + 1)])
    errors = np.array([h.GetBinError(i) for i in range(1, n + 1)])
    edges = np.array([h.GetBinLowEdge(i) for i in range(1, n + 2)])
    return centers, contents, errors, edges

def plot_scan_with_fit(gr, f, xlabel, ylabel, title, outpath):
    x, y, ex, ey = root_graph_to_arrays(gr)
    xmin, xmax = f.GetXmin(), f.GetXmax()
    xf, yf = root_tf1_to_curve(f, xmin, xmax)

    fig, ax = plt.subplots(figsize=(9, 6))
    
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    # Opzionale: migliora lo stile (es. tacche rivolte verso l'interno, visibili su tutti i lati)
    ax.tick_params(which='both', direction='in', top=True, right=True)
    ax.tick_params(which='minor', length=3)
    ax.tick_params(which='major', length=5)
    ax.errorbar(x, y, yerr=ey, xerr=ex if ex.any() else None,
                fmt='o', ms=3, color='black', ecolor='gray',
                elinewidth=0.8, capsize=0, label='data')
    ax.plot(xf, yf, color='red', lw=1.8, label='fit')
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(outpath, dpi=150)
    plt.close(fig)
    return fig

def plot_slice_with_gaus_fit(h, f, xlabel, title, outpath, color='tab:blue'):
    centers, contents, errors, edges = root_hist_to_arrays(h)
    xmin, xmax = f.GetXmin(), f.GetXmax()
    xf, yf = root_tf1_to_curve(f, xmin, xmax)

    fig, ax = plt.subplots(figsize=(8, 6))

    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    # Opzionale: migliora lo stile (es. tacche rivolte verso l'interno, visibili su tutti i lati)
    ax.tick_params(which='both', direction='in', top=True, right=True)
    ax.tick_params(which='minor', length=3)
    ax.tick_params(which='major', length=5)
    ax.stairs(contents, edges, fill=True, color=color, alpha=0.55, label='data')
    ax.plot(xf, yf, color='red', lw=1.8,
            label=f"gaus fit ($\\sigma$ = {f.GetParameter(2):.1f} $\\mu$rad)")
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Counts")
    ax.set_title(title)
    ax.legend(frameon=False, loc='upper right')
    fig.tight_layout()
    fig.savefig(outpath, dpi=150)
    plt.close(fig)
    return fig



ROOT.ROOT.EnableImplicitMT() 
ROOT.gStyle.SetOptStat(0)

file = 8656
filename = "recoDataSimple_" + str(file) + "_xtalMerging.root"
files = ["recoDataSimple_8430_xtalMerging.root", "recoDataSimple_8431_xtalMerging.root"]

ROOT.gStyle.SetPalette(ROOT.kBird)

# VARIABLES
width = 12.8 # width of the crystal in mm (for spatial cut)
height = 2 # height of the crystal in mm (for spatial cut)
max_value = 8000 # histograms range
deflection_peak = 6010 # urad
fit_range = 150 # urad
minimum_entries = 500

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
    centers, sigmas, sigma_errs, entries = [], [], [], []
 
    n_bins = h2v.GetNbinsX()
    for ix in range(1, n_bins + 1):
        h1 = h2v.ProjectionY(f"py_{h2v.GetName()}_{ix}", ix, ix)

        e = h1.GetEntries()
        if e < min_entries:
            continue
        f = ROOT.TF1(f"f_{h2v.GetName()}_{ix}", "gaus", -range_fit_gaus, range_fit_gaus)
        f.SetParameters(h1.GetMaximum(), h1.GetMean(), max(h1.GetRMS(), 5.0))
        h1.Fit(f, "RQ0")
 
        centers.append(h2v.GetXaxis().GetBinCenter(ix))
        sigmas.append(f.GetParameter(2))
        sigma_errs.append(f.GetParError(2))
        entries.append(e)

    return centers, sigmas, sigma_errs, entries
 
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

def fit_four_step_edges(centers, sigmas, sigma_errs,
                         e1_guess, e2_guess, e3_guess, e4_guess,
                         baseline_guess=None, clamp_level_guess=None, crystal_level_guess=None,
                         transition_guess=0.1):
    n = len(centers)
    if n == 0:
        print("ERROR: Zero data points! The previous cut filtered out all events.")
        return None, None, None, None
    gr = ROOT.TGraphErrors(n, array('d', centers), array('d', sigmas),
                            array('d', [0.0] * n), array('d', sigma_errs))

    if baseline_guess is None:
        baseline_guess = min(sigmas)
    if clamp_level_guess is None:
        clamp_level_guess = max(sigmas)
    if crystal_level_guess is None:
        crystal_level_guess = sorted(sigmas)[len(sigmas) // 2]

    formula = ("[0]"
               " + ([1]-[0])*0.5*(1+TMath::Erf((x-[3])/[7]))"   # e1: vuoto->clamp
               " + ([2]-[1])*0.5*(1+TMath::Erf((x-[4])/[7]))"   # e2: clamp->crystal
               " + ([1]-[2])*0.5*(1+TMath::Erf((x-[5])/[7]))"   # e3: crystal->clamp
               " + ([0]-[1])*0.5*(1+TMath::Erf((x-[6])/[7]))")  # e4: clamp->vuoto

    f = ROOT.TF1("f_four_step", formula, min(centers), max(centers))
    f.SetParameters(baseline_guess, clamp_level_guess, crystal_level_guess,
                     e1_guess, e2_guess, e3_guess, e4_guess, transition_guess)
    f.SetParLimits(7, 0.001, (max(centers) - min(centers)) / 4.0)
    for i in (3, 4, 5, 6):
        f.SetParLimits(i, min(centers), max(centers))

    gr.Fit(f, "RQ")

    edges = {"e1": (f.GetParameter(3), f.GetParError(3)),
             "e2": (f.GetParameter(4), f.GetParError(4)),
             "e3": (f.GetParameter(5), f.GetParError(5)),
             "e4": (f.GetParameter(6), f.GetParError(6))}
    levels = {"baseline": (f.GetParameter(0), f.GetParError(0)),
              "clamp": (f.GetParameter(1), f.GetParError(1)),
              "crystal": (f.GetParameter(2), f.GetParError(2))}
    return edges, levels, f, gr

 
# # y edges (scan y, using the full x range -- or restrict to rough x window from the old method / previous run)
# h2_y = book_scan_histogram(df_phys, "Tracks.d0Out_y", "Deltatheta_y", scan_min=-15, scan_max=9, n_bins=150, slice_var="Tracks.d0_x", slice_min=-1.5, slice_max=1.5) # 0.01 mm precision
# y_centers, y_sigmas, y_sigma_errs, y_entries = extract_widths_from_h2(h2_y)
# y_edges, y_levels, f_y, gr_y = fit_four_step_edges(y_centers, y_sigmas, y_sigma_errs, e1_guess=-11.5, e2_guess=-6.0, e3_guess=7.0, e4_guess=8.5)
 
# # x edges (scan x, restricted to the just-found y window)
# h2_x = book_scan_histogram(df_phys, "Tracks.d0Out_x", "Deltatheta_x", scan_min=-3, scan_max=4, n_bins=140, slice_var="Tracks.d0_y", slice_min=y_edges["e2"][0], slice_max=y_edges["e3"][0])
# x_centers, x_sigmas, x_sigma_errs, _ = extract_widths_from_h2(h2_x)
# x_lo, x_hi, x_lo_err, x_hi_err, f_x, gr_x = fit_step_edges(x_centers, x_sigmas, x_sigma_errs, edge_lo_guess=-1.0, edge_hi_guess=1.0, transition_guess=0.02)

# h2_y = book_scan_histogram(df_phys, "Tracks.d0Out_y", "Deltatheta_y", scan_min=-15, scan_max=9, n_bins=150, slice_var="Tracks.d0_x", slice_min=x_lo, slice_max=x_hi)
# y_centers, y_sigmas, y_sigma_errs, _ = extract_widths_from_h2(h2_y)
# y_edges, y_levels, f_y, gr_y = fit_four_step_edges(y_centers, y_sigmas, y_sigma_errs, e1_guess=y_edges["e1"][0], e2_guess=y_edges["e2"][0], e3_guess=y_edges["e3"][0], e4_guess=y_edges["e4"][0])

# h2_x = book_scan_histogram(df_phys, "Tracks.d0Out_x", "Deltatheta_x", scan_min=-3, scan_max=4, n_bins=140, slice_var="Tracks.d0_y", slice_min=y_edges["e2"][0], slice_max=y_edges["e3"][0])
# x_centers, x_sigmas, x_sigma_errs, _ = extract_widths_from_h2(h2_x)
# x_lo, x_hi, x_lo_err, x_hi_err, f_x, gr_x = fit_step_edges(x_centers, x_sigmas, x_sigma_errs, edge_lo_guess=-1.0, edge_hi_guess=1.0, transition_guess=0.02)

# y_lo, y_lo_err = y_edges["e2"]
# y_hi, y_hi_err = y_edges["e3"]

x_guess = (-1.0, 1.0)
y_guess = (-11.5, -6.0, 7.0, 8.5)

def iterate_until_converged(df_phys, x_guess, y_guess, tol=1e-3, max_iter=8):
    x_lo, x_hi = x_guess
    y_e1, y_e2, y_e3, y_e4 = y_guess
    h2_x, h2_y = None, None

    best_shift = float("inf")
    best_state = None
    prev_shift = None
    diverging_count = 0

    for it in range(max_iter):
        h2_x = book_scan_histogram(df_phys, "Tracks.d0Out_x", "Deltatheta_x", scan_min=-3, scan_max=4, n_bins=140,
                                    slice_var="Tracks.d0_y", slice_min=y_e2, slice_max=y_e3)
        xc, xs, xe, _ = extract_widths_from_h2(h2_x)
        if len(xc) == 0:
            print(f"iter {it}: x-scan produced zero bins (slice y=[{y_e2:.3f},{y_e3:.3f}]) -- stopping, keeping last good state.")
            break

        x_lo_new, x_hi_new, x_lo_err, x_hi_err, f_x, gr_x = fit_step_edges(xc, xs, xe, edge_lo_guess=x_lo, edge_hi_guess=x_hi, transition_guess=0.02)
        if f_x is None:
            print(f"iter {it}: x fit failed -- stopping, keeping last good state.")
            break
        if file in [8430, 8431, 8650]:
            scan_min_y, scan_max_y = -15, 9
        elif file in [8655, 8656]:
            scan_min_y, scan_max_y = -15, 15
        h2_y = book_scan_histogram(df_phys, "Tracks.d0Out_y", "Deltatheta_y", scan_min=scan_min_y, scan_max=scan_max_y, n_bins=150,
                                    slice_var="Tracks.d0_x", slice_min=x_lo_new, slice_max=x_hi_new)
        yc, ys, yerr, _ = extract_widths_from_h2(h2_y)
        if len(yc) == 0:
            print(f"iter {it}: y-scan produced zero bins (slice x=[{x_lo_new:.3f},{x_hi_new:.3f}]) -- stopping, keeping last good state.")
            break

        y_edges, y_levels, f_y, gr_y = fit_four_step_edges(yc, ys, yerr, e1_guess=y_e1, e2_guess=y_e2, e3_guess=y_e3, e4_guess=y_e4)
        if y_edges is None:
            print(f"iter {it}: y fit failed -- stopping, keeping last good state.")
            break

        y_e2_new, y_e3_new = y_edges["e2"][0], y_edges["e3"][0]

        # sanity: ordering must make physical sense
        if not (y_edges["e1"][0] < y_e2_new < y_e3_new < y_edges["e4"][0]):
            print(f"iter {it}: WARNING -- edges out of order (e1={y_edges['e1'][0]:.3f}, "
                  f"e2={y_e2_new:.3f}, e3={y_e3_new:.3f}, e4={y_edges['e4'][0]:.3f}) -- stopping.")
            break
        if not (x_lo_new < x_hi_new):
            print(f"iter {it}: WARNING -- x edges out of order (x_lo={x_lo_new:.3f}, x_hi={x_hi_new:.3f}) -- stopping.")
            break

        shift = max(abs(x_lo_new - x_lo), abs(x_hi_new - x_hi),
                    abs(y_e2_new - y_e2), abs(y_e3_new - y_e3))
        print(f"iter {it}: max shift = {shift:.5f} mm")

        # keep the best (smallest-shift) state seen so far as a fallback
        if shift < best_shift:
            best_shift = shift
            best_state = (x_lo_new, x_hi_new, x_lo_err, x_hi_err, f_x, gr_x, h2_x,
                          y_edges, y_levels, f_y, gr_y, h2_y)

        # divergence detection: shift growing for 2 iterations in a row
        if prev_shift is not None and shift > prev_shift:
            diverging_count += 1
            if diverging_count >= 2:
                print(f"iter {it}: shift increasing for {diverging_count} iterations in a row -- "
                      f"diverging. Falling back to best state (shift={best_shift:.5f} at that point).")
                break
        else:
            diverging_count = 0
        prev_shift = shift

        x_lo, x_hi = x_lo_new, x_hi_new
        y_e1, y_e2, y_e3, y_e4 = y_edges["e1"][0], y_e2_new, y_e3_new, y_edges["e4"][0]

        if shift < tol:
            print(f"Converged after {it+1} iterations.")
            best_state = (x_lo, x_hi, x_lo_err, x_hi_err, f_x, gr_x, h2_x, y_edges, y_levels, f_y, gr_y, h2_y)
            break

    if best_state is None:
        raise RuntimeError("iterate_until_converged: never obtained a valid fit state.")

    return best_state

x_lo, x_hi, x_lo_err, x_hi_err, f_x, gr_x, h2_x, y_edges, y_levels, f_y, gr_y, h2_y = iterate_until_converged(df_phys, x_guess, y_guess)

y_lo, y_lo_err = y_edges["e2"]
y_hi, y_hi_err = y_edges["e3"]

mean_res_x = df_phys.Mean("DeltathetaErr_x").GetValue() #they are all the same but safer to take the mean of the distribution
mean_res_y = df_phys.Mean("DeltathetaErr_y").GetValue()
print("="*50)
print(f"Mean resolution (x): {mean_res_x:.2f} urad")
print(f"Mean resolution (y): {mean_res_y:.2f} urad")

sigma_bsl_x, sigma_jump_x = f_x.GetParameter(0), f_x.GetParameter(1)
sigma_bsl_y, sigma_clamp_y, sigma_crystal_y = f_y.GetParameter(0), f_y.GetParameter(1), f_y.GetParameter(2)

theta_crystal_x = math.sqrt(max((sigma_bsl_x + sigma_jump_x)**2 - sigma_bsl_x**2, 0))
theta_crystal_y = math.sqrt(max(sigma_crystal_y**2 - sigma_bsl_y**2, 0))
theta_clamp_y   = math.sqrt(max(sigma_clamp_y**2   - sigma_bsl_y**2, 0))

err_bsl_x, err_jump_x = f_x.GetParError(0), f_x.GetParError(1)
err_bsl_y, err_clamp_y, err_crystal_y = f_y.GetParError(0), f_y.GetParError(1), f_y.GetParError(2)

err_theta_crystal_x = 0.0
if theta_crystal_x > 0:
    deriv_bsl_x  = sigma_jump_x / theta_crystal_x
    deriv_jump_x = (sigma_bsl_x + sigma_jump_x) / theta_crystal_x
    err_theta_crystal_x = math.sqrt((deriv_bsl_x * err_bsl_x)**2 + (deriv_jump_x * err_jump_x)**2)
err_theta_crystal_y = 0.0
if theta_crystal_y > 0:
    deriv_crystal_y = sigma_crystal_y / theta_crystal_y
    deriv_bsl_y     = sigma_bsl_y / theta_crystal_y
    err_theta_crystal_y = math.sqrt((deriv_crystal_y * err_crystal_y)**2 + (deriv_bsl_y * err_bsl_y)**2)
err_theta_clamp_y = 0.0
if theta_clamp_y > 0:
    deriv_clamp_y     = sigma_clamp_y / theta_clamp_y
    deriv_bsl_y_clamp = sigma_bsl_y / theta_clamp_y
    err_theta_clamp_y = math.sqrt((deriv_clamp_y * err_clamp_y)**2 + (deriv_bsl_y_clamp * err_bsl_y)**2)

# we are comapring the fitted erf transition point with the mean of the d0Err_x/y distributions, which is a good sanity check
mean_d0err_x = df_phys.Mean("Tracks.d0Err_x").GetValue()
mean_d0err_y = df_phys.Mean("Tracks.d0Err_y").GetValue()

print(f"fitted transition x = {f_x.GetParameter(4):.4f} ± {f_x.GetParError(4):.4f} mm vs mean d0Err_x = {mean_d0err_x:.4f} mm")
print(f"fitted transition y = {f_y.GetParameter(7):.4f} ± {f_y.GetParError(7):.4f} mm vs mean d0Err_y = {mean_d0err_y:.4f} mm")
print("="*50)

print("x:")
print(f"\tbaseline sigma x  = {sigma_bsl_x:.2f} ± {f_x.GetParError(0):.2f} urad")
print(f"\tjump sigma x      = {sigma_jump_x:.2f} ± {f_x.GetParError(1):.2f} urad")
print(f"\tsigma mcs x = {theta_crystal_x:.2f} ± {err_theta_crystal_x:.2f} urad")
print('='*50)
print("y:")
print(f"\tbaseline sigma y  = {sigma_bsl_y:.2f} ± {f_y.GetParError(0):.2f} urad")
print(f"\tsigma mcs y = {theta_crystal_y:.2f} ± {err_theta_crystal_y:.2f} urad")
print(f"\tsigma mcs clamp y = {theta_clamp_y:.2f} ± {err_theta_clamp_y:.2f} urad")

os.makedirs(f"plots_{file}_edges", exist_ok=True)

fig_y = plot_scan_with_fit(
    gr_y, f_y,
    xlabel="impact position y [mm]", ylabel=r"$\sigma(\Delta\theta_y)$ [$\mu$rad]",
    title="Local scattering width vs y",
    outpath=f"plots_{file}_edges/scattering_width_vs_y.png"
)

fig_x = plot_scan_with_fit(
    gr_x, f_x,
    xlabel="impact position x [mm]", ylabel=r"$\sigma(\Delta\theta_x)$ [$\mu$rad]",
    title="Local scattering width vs x",
    outpath=f"plots_{file}_edges/scattering_width_vs_x.png"
)

 
print("=" * 50)
print(f"Crystal footprint from scattering method:")
print(f"\tx = [{x_lo:.4f} ± {x_lo_err:.4f}, {x_hi:.4f} ± {x_hi_err:.4f}] mm  (width = {x_hi - x_lo:.4f} ± {pow(x_lo_err*x_lo_err + x_hi_err*x_hi_err,0.5):.4f} mm)")
print(f"\ty = [{y_lo:.4f} ± {y_lo_err:.4f}, {y_hi:.4f} ± {y_hi_err:.4f}] mm  (width = {y_hi - y_lo:.4f} ± {pow(y_lo_err*y_lo_err + y_hi_err*y_hi_err,0.5):.4f} mm)")
print(f"\ty clamp position = {y_edges['e1'][0]:.4f} ± {y_edges['e1'][1]:.4f} mm")

# Plot della mappa 2D Posizione x/y vs deflessione 
# c_2d_y = ROOT.TCanvas("c_2d_y", "Position Y vs Deflection", 900, 600)
h2_y_val = h2_y.GetValue()
h2_y_val.SetTitle("Posizione Y vs Deflessione #Delta#theta_{y}; impact position y [mm]; #Delta#theta_{x} [#murad]")
# h2_y_val.Draw("colz")
# c_2d_y.Update()

# c_2d_x = ROOT.TCanvas("c_2d_x", "Position X vs Deflection", 900, 600)
h2_x_val = h2_x.GetValue()
h2_x_val.SetTitle("Posizione X vs Deflessione #Delta#theta_{x}; impact position x [mm]; #Delta#theta_{x} [#murad]")
# h2_x_val.Draw("colz")
# c_2d_x.Update()

# Estrazione e plot di DUE singole fettine (Una DENTRO e una FUORI dal cristallo)
# Scegliamo due bin a caso basandoci sui centri Y trovati
# bin_dentro = h2_y_val.GetXaxis().FindBin((y_lo + y_hi) / 2.0) # Esattamente a metà cristallo
bin_dentro = h2_y_val.GetXaxis().FindBin(0)
bin_clamp = h2_y_val.GetXaxis().FindBin(-10)
bin_fuori = h2_y_val.GetXaxis().FindBin(y_edges["e1"][0] - 1.0)
# bin_fuori = h2_y_val.GetXaxis().FindBin(-12.0)

# Estraiamo gli istogrammi 1D
h1_dentro = h2_y_val.ProjectionY("h1_dentro", bin_dentro, bin_dentro)
h1_clamp = h2_y_val.ProjectionY("h1_clamp", bin_clamp, bin_clamp)
h1_clamp.Rebin(5) # number of bin gets divided by this factor
h1_fuori = h2_y_val.ProjectionY("h1_fuori", bin_fuori, bin_fuori)

# Fit gaussiano sulle 3 fette (dentro/clamp/fuori) -- fit resta in ROOT come sempre
f_in = ROOT.TF1("f_in", "gaus", -fit_range, fit_range)
h1_dentro.Fit(f_in, "RQ0")

f_clamp = ROOT.TF1("f_clamp", "gaus", -fit_range, fit_range)
h1_clamp.Fit(f_clamp, "RQ0")

f_out = ROOT.TF1("f_out", "gaus", -fit_range, fit_range)
h1_fuori.Fit(f_out, "RQ0")

fig_in = plot_slice_with_gaus_fit(
    h1_dentro, f_in, xlabel=r"$\Delta\theta_x$ [$\mu$rad]",
    title=f"Crystal slice (Y = {h2_y_val.GetXaxis().GetBinCenter(bin_dentro):.2f} mm)",
    outpath=f"plots_{file}_edges/slice_inside_y.png", color='tab:blue'
)

fig_clamp = plot_slice_with_gaus_fit(
    h1_clamp, f_clamp, xlabel=r"$\Delta\theta_x$ [$\mu$rad]",
    title=f"Crystal slice - clamp (Y = {h2_y_val.GetXaxis().GetBinCenter(bin_clamp):.2f} mm)",
    outpath=f"plots_{file}_edges/slice_clamp_y.png", color='tab:green'
)

fig_out = plot_slice_with_gaus_fit(
    h1_fuori, f_out, xlabel=r"$\Delta\theta_x$ [$\mu$rad]",
    title=f"Outside slice (Y = {h2_y_val.GetXaxis().GetBinCenter(bin_fuori):.2f} mm)",
    outpath=f"plots_{file}_edges/slice_outside_y.png", color='tab:pink'
)

# Scegliamo un bin al centro del cristallo e uno 1.5 mm fuori dal bordo destro
bin_dentro_x = h2_x_val.GetXaxis().FindBin((x_lo + x_hi) / 2.0) 
bin_fuori_x  = h2_x_val.GetXaxis().FindBin(x_hi + 1.5)           

# Estraiamo gli istogrammi 1D
h1_dentro_x = h2_x_val.ProjectionY("h1_dentro_x", bin_dentro_x, bin_dentro_x)
h1_fuori_x  = h2_x_val.ProjectionY("h1_fuori_x", bin_fuori_x, bin_fuori_x)

# Fit gaussiano sulle 2 fette lungo x (dentro/fuori)
f_in_x = ROOT.TF1("f_in_x", "gaus", -fit_range, fit_range)
h1_dentro_x.Fit(f_in_x, "RQ0")

f_out_x = ROOT.TF1("f_out_x", "gaus", -fit_range, fit_range)
h1_fuori_x.Fit(f_out_x, "RQ0")

fig_in_x = plot_slice_with_gaus_fit(
    h1_dentro_x, f_in_x, xlabel=r"$\Delta\theta_x$ [$\mu$rad]",
    title=f"Crystal slice (X = {h2_x_val.GetXaxis().GetBinCenter(bin_dentro_x):.2f} mm)",
    outpath=f"plots_{file}_edges/slice_inside_x.png", color='tab:blue'
)

fig_out_x = plot_slice_with_gaus_fit(
    h1_fuori_x, f_out_x, xlabel=r"$\Delta\theta_x$ [$\mu$rad]",
    title=f"Outside slice (X = {h2_x_val.GetXaxis().GetBinCenter(bin_fuori_x):.2f} mm)",
    outpath=f"plots_{file}_edges/slice_outside_x.png", color='tab:pink'
)

print(f"\nAll matplotlib plots saved under plots_{file}_edges/")