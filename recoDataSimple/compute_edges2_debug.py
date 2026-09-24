import ROOT
import math
import os
import numpy as np
import matplotlib.pyplot as plt
from array import array
from plotting_utils import plot_histo1d, tf1_to_curve

ROOT.ROOT.EnableImplicitMT()
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetPalette(ROOT.kBird)

file = 8430
filename = "data/recoDataSimple_" + str(file) + "_xtalMerging.root"

# VARIABLES
fit_range = 150       # urad, range of the gaussian fit on each slice
minimum_entries = 500 # min entries in a slice to attempt a fit

os.makedirs(f"plots_{file}_edges2_debug", exist_ok=True)

# ROOT DATA FRAME
df = ROOT.RDataFrame("simpleEvent", filename)
print("=" * 50)
print(f"Analyzing {filename} ...")

# FILTERING the data: single tracks and conversion from rad to urad
df_phys = df.Filter("SingleTrack == 1")
df_phys = df_phys.Define("thetaIn_x", "Tracks.thetaIn_x * 1e6") \
    .Define("Deltatheta_x", "(Tracks.thetaOut_x - Tracks.thetaIn_x) * 1e6") \
    .Define("Deltatheta_y", "(Tracks.thetaOut_y - Tracks.thetaIn_y) * 1e6") \
    .Define("DeltathetaErr_x", "sqrt(Tracks.thetaInErr_x * Tracks.thetaInErr_x + Tracks.thetaOutErr_x * Tracks.thetaOutErr_x) * 1e6") \
    .Define("DeltathetaErr_y", "sqrt(Tracks.thetaInErr_y * Tracks.thetaInErr_y + Tracks.thetaOutErr_y * Tracks.thetaOutErr_y) * 1e6")
print("=" * 50)


def book_scan_histogram(df_in, scan_var, delta_var, scan_min, scan_max, n_bins,
                         slice_var=None, slice_min=None, slice_max=None, n_dtheta_bins=300):
    df_slice = df_in
    if slice_var is not None:
        df_slice = df_slice.Filter(f"{slice_var} > {slice_min} && {slice_var} < {slice_max}")

    h2 = df_slice.Histo2D((f"h2_scat_{scan_var}", "", n_bins, scan_min, scan_max,
                           n_dtheta_bins, -fit_range, fit_range), scan_var, delta_var)
    return h2


def extract_and_plot_slices(h2_lazy, tag, axis_label, rebin=1,
                             range_fit_gaus=fit_range, min_entries=minimum_entries, save=None,
                             double_peak=None):
    """
    For every bin along the scan axis of a booked TH2 (h2_lazy), project onto the
    delta-theta axis and fit it, then draw every slice's histogram + fit(s) into one
    grid figure (as many subplots as slices with enough statistics).

    By default a single gaussian is fitted over the full [-range_fit_gaus, range_fit_gaus]
    range. Pass `double_peak` to instead fit TWO independent gaussians, each restricted
    to its own side of a separating cut, for slices whose center falls inside a given
    scan-position window (used where MCS and volume-reflection peaks overlap):

        double_peak = {
            "range":      (pos_min, pos_max),  # scan-position window [mm] to apply this to
            "cut":        urad_value,          # boundary between the two peaks
            "guess_low":  mean guess for the fit on [-range_fit_gaus, cut]
            "guess_high": mean guess for the fit on [cut, +range_fit_gaus]
            "primary":    "low" or "high"      # which peak's sigma feeds the returned sigmas[]
        }

    Returns: centers, sigmas, sigma_errs, entries
    (sigmas/sigma_errs come from the single fit, or from the "primary" peak when double_peak applies)
    """
    h2v = h2_lazy.GetValue() if hasattr(h2_lazy, "GetValue") else h2_lazy
    n_bins = h2v.GetNbinsX()

    centers, sigmas, sigma_errs, entries = [], [], [], []
    slice_histos = []  # (center, h1, [(TF1, (xmin, xmax)), ...]) kept only for fitted slices

    for ix in range(1, n_bins + 1):
        h1 = h2v.ProjectionY(f"py_{h2v.GetName()}_{ix}", ix, ix)
        if rebin and rebin > 1:
            h1.Rebin(rebin)

        e = h1.GetEntries()
        if e < min_entries:
            continue

        center = h2v.GetXaxis().GetBinCenter(ix)

        use_double = (double_peak is not None and
                      double_peak["range"][0] <= center <= double_peak["range"][1])

        if use_double:
            cut = double_peak["cut"]

            f_lo = ROOT.TF1(f"f_{h2v.GetName()}_{ix}_lo", "gaus", -range_fit_gaus, cut)
            f_lo.SetParameters(h1.GetMaximum(), double_peak["guess_low"], 15.0)
            h1.Fit(f_lo, "RQ0")

            f_hi = ROOT.TF1(f"f_{h2v.GetName()}_{ix}_hi", "gaus", cut, range_fit_gaus)
            f_hi.SetParameters(h1.GetMaximum(), double_peak["guess_high"], 15.0)
            h1.Fit(f_hi, "RQ0")

            fits = [(f_lo, (-range_fit_gaus, cut)), (f_hi, (cut, range_fit_gaus))]
            primary = f_lo if double_peak["primary"] == "low" else f_hi
            sigma, sigma_err = primary.GetParameter(2), primary.GetParError(2)
        else:
            f = ROOT.TF1(f"f_{h2v.GetName()}_{ix}", "gaus", -range_fit_gaus, range_fit_gaus)
            f.SetParameters(h1.GetMaximum(), h1.GetMean(), max(h1.GetRMS(), 5.0))
            h1.Fit(f, "RQ0")

            fits = [(f, (-range_fit_gaus, range_fit_gaus))]
            sigma, sigma_err = f.GetParameter(2), f.GetParError(2)

        centers.append(center)
        sigmas.append(sigma)
        sigma_errs.append(sigma_err)
        entries.append(e)
        slice_histos.append((center, h1, fits))

    if save is not None:
        _plot_slice_grid(slice_histos, axis_label, save)

    return centers, sigmas, sigma_errs, entries


def _plot_slice_grid(slice_histos, axis_label, save):
    """Draw one subplot per slice (histogram + gaussian fit(s)), arranged in a
    roughly square grid, reusing plot_histo1d for the histogram and overlaying
    each fit's curve on its own restricted range."""
    n = len(slice_histos)
    if n == 0:
        print(f"WARNING: no slices to plot for '{save}' (all below min_entries).")
        return

    ncols = math.ceil(math.sqrt(n))
    nrows = math.ceil(n / ncols)

    fig, axs = plt.subplots(nrows, ncols, figsize=(3.2 * ncols, 2.6 * nrows), squeeze=False)
    axs_flat = axs.flatten()

    fit_colors = ["red", "darkorange"]

    for ax, (center, h1, fits) in zip(axs_flat, slice_histos):
        plot_histo1d(h1, fit_func=None, ax=ax, style="fill", color="tab:blue",
                     xlabel="", ylabel="")

        for (f, (xmin, xmax)), color in zip(fits, fit_colors):
            fx, fy = tf1_to_curve(f, xmin, xmax)
            ax.plot(fx, fy, color=color, lw=1.5)

        if len(fits) == 2:
            s_lo = fits[0][0].GetParameter(2)
            s_hi = fits[1][0].GetParameter(2)
            title = f"{axis_label} = {center:.2f}   $\\sigma_{{lo}}$={s_lo:.1f}  $\\sigma_{{hi}}$={s_hi:.1f}"
        else:
            title = f"{axis_label} = {center:.2f}   $\\sigma$ = {fits[0][0].GetParameter(2):.1f} $\\mu$rad"

        ax.set_title(title, fontsize=9)

    for ax in axs_flat[n:]:
        ax.axis("off")

    fig.tight_layout()
    fig.savefig(save, dpi=150)
    # fig.savefig(save.replace(".pdf", ".png"), dpi=150)
    plt.close(fig)
    print(f"Saved slice grid ({n} panels, {nrows}x{ncols}) -> {save}")


# ===================== SCANS =====================
# x-scan: Deltatheta_x vs Tracks.d0Out_x, restricted to a rough y window around the crystal
h2_x = book_scan_histogram(df_phys, "Tracks.d0Out_x", "Deltatheta_x", scan_min=-3, scan_max=4, n_bins=140,
                            slice_var="Tracks.d0_y", slice_min=-6.0, slice_max=7.0)
# in [-1.1, 1] mm (inside the crystal) Deltatheta_x shows two peaks: MCS around zero
# and volume reflection around -40/-50 urad. Separate them at -25 urad: the low-side
# fit (far-left tail, up to -25) catches the VR peak, the high-side fit (from -25 up)
# catches the zero-centered MCS peak we actually care about.
double_peak_x = {
    "range": (-1.1, 1.0),
    "cut": -25.0,
    "guess_low": -45.0,
    "guess_high": 0.0,
    "primary": "high",  # the zero-centered (MCS) peak
}

x_centers, x_sigmas, x_sigma_errs, _ = extract_and_plot_slices(
    h2_x, tag="x", axis_label="x", rebin=3,
    save=f"plots_{file}_edges2_debug/slices_grid_x.pdf",
    double_peak=double_peak_x)

# y-scan: Deltatheta_y vs Tracks.d0Out_y, over the full illuminated range
if file in [8430, 8431, 8650]:
    scan_min_y, scan_max_y = -15, 9
elif file in [8655, 8656]:
    scan_min_y, scan_max_y = -15, 15
else:
    scan_min_y, scan_max_y = -15, 15

h2_y = book_scan_histogram(df_phys, "Tracks.d0Out_y", "Deltatheta_y", scan_min=scan_min_y, scan_max=scan_max_y,
                            n_bins=100, slice_var="Tracks.d0_x", slice_min=-1.0, slice_max=1.0)
# in [-1, 4] mm (inside the crystal) Deltatheta_y shows two peaks: MCS around zero and
# volume reflection around +50 urad. Separate them at +25 urad: the low-side fit (up to
# +25) catches the zero-centered MCS peak we care about, the high-side fit (from +25 up)
# catches the VR peak.
double_peak_y = {
    "range": (-1.0, 4.0),
    "cut": 25.0,
    "guess_low": 0.0,
    "guess_high": 50.0,
    "primary": "low",  # the zero-centered (MCS) peak
}

y_centers, y_sigmas, y_sigma_errs, _ = extract_and_plot_slices(
    h2_y, tag="y", axis_label="y", rebin=3,
    save=f"plots_{file}_edges2_debug/slices_grid_y.pdf",
    double_peak=double_peak_y)

# ===================== RESOLUTIONS (kept, unrelated to edge fitting) =====================
mean_res_x = df_phys.Mean("DeltathetaErr_x").GetValue()
mean_res_y = df_phys.Mean("DeltathetaErr_y").GetValue()
print("=" * 50)
print(f"Mean resolution (x): {mean_res_x:.2f} urad - (y): {mean_res_y:.2f} urad")
print(f"Mean resolution combined: {np.sqrt(mean_res_y**2 + mean_res_x**2):.2f} urad")

# ===================== sigma(x) and sigma(y) profile plots (data points only, no step fit) =====================
fig_x, ax_x = plt.subplots(figsize=(8, 6))
ax_x.errorbar(x_centers, x_sigmas, yerr=x_sigma_errs, fmt='o', color='black',
              ecolor='gray', markersize=3, elinewidth=0.8)
ax_x.set_xlabel("x [mm]")
ax_x.set_ylabel(r"$\sigma(\Delta\theta_x)$ [$\mu$rad]")
ax_x.set_title("Local scattering width vs x")
fig_x.tight_layout()
# fig_x.savefig(f"plots_{file}_edges2_debug/scattering_width_vs_x.pdf")
fig_x.savefig(f"plots_{file}_edges2_debug/scattering_width_vs_x.png")
plt.close(fig_x)

fig_y, ax_y = plt.subplots(figsize=(8, 6))
ax_y.errorbar(y_centers, y_sigmas, yerr=y_sigma_errs, fmt='o', color='black',
              ecolor='gray', markersize=3, elinewidth=0.8)
ax_y.set_xlabel("y [mm]")
ax_y.set_ylabel(r"$\sigma(\Delta\theta_y)$ [$\mu$rad]")
ax_y.set_title("Local scattering width vs y")
fig_y.tight_layout()
# fig_y.savefig(f"plots_{file}_edges2_debug/scattering_width_vs_y.pdf")
fig_y.savefig(f"plots_{file}_edges2_debug/scattering_width_vs_y.png")
plt.close(fig_y)

print(f"\nAll matplotlib plots saved under plots_{file}_edges2_debug/")