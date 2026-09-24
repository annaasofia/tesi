import ROOT
import math
import os
import numpy as np
import matplotlib.pyplot as plt
from plotting_utils import plot_histo1d

ROOT.ROOT.EnableImplicitMT()
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetPalette(ROOT.kBird)

file = 8430
filename = "data/recoDataSimple_" + str(file) + "_xtalMerging.root"

# VARIABLES
fit_range = 150       # urad, range of the gaussian fit on each slice (Deltatheta_x/y)
theta_L = 13/2 # urad
minimum_entries = 500 # min entries in a position slice to attempt a fit

PLOT_DIR = f"plots_{file}_edges2_debug"
os.makedirs(PLOT_DIR, exist_ok=True)

# ROOT DATA FRAME
df = ROOT.RDataFrame("simpleEvent", filename)
print("=" * 50)
print(f"Analyzing {filename} ...")

# FILTERING the data: single tracks and conversion from rad to urad
df_phys = df.Filter("SingleTrack == 1")
df_phys = df_phys.Define("thetaIn_x", "Tracks.thetaIn_x * 1e6")
df_phys = df_phys.Define("Deltatheta_x", "(Tracks.thetaOut_x - Tracks.thetaIn_x) * 1e6") \
    .Define("Deltatheta_y", "(Tracks.thetaOut_y - Tracks.thetaIn_y) * 1e6") \
    .Define("DeltathetaErr_x", "sqrt(Tracks.thetaInErr_x * Tracks.thetaInErr_x + Tracks.thetaOutErr_x * Tracks.thetaOutErr_x) * 1e6") \
    .Define("DeltathetaErr_y", "sqrt(Tracks.thetaInErr_y * Tracks.thetaInErr_y + Tracks.thetaOutErr_y * Tracks.thetaOutErr_y) * 1e6")
df_filtered = df_phys.Filter(f"abs(thetaIn_x) < {theta_L}")

print("=" * 50)


# ===================== GENERIC HELPERS =====================

def book_scan_histogram(df_in, scan_var, delta_var, scan_min, scan_max, n_bins,
                         slice_var=None, slice_min=None, slice_max=None,
                         delta_min=-fit_range, delta_max=fit_range, n_delta_bins=300):
    """Book a lazy TH2(scan_var, delta_var), optionally restricted to a slice_var window.
    Generic: works both for a position scan (scan_var=position, delta_var=Deltatheta)
    and for a Deltatheta scan (scan_var=Deltatheta_x, delta_var=thetaIn_x, see idea 2)."""
    df_slice = df_in
    if slice_var is not None:
        df_slice = df_slice.Filter(f"{slice_var} > {slice_min} && {slice_var} < {slice_max}")

    h2 = df_slice.Histo2D((f"h2_scat_{scan_var}", "", n_bins, scan_min, scan_max,
                           n_delta_bins, delta_min, delta_max), scan_var, delta_var)
    return h2


def extract_and_plot_slices(h2_lazy, axis_label, save, rebin=1,
                             fit_range_local=fit_range, min_entries=minimum_entries,
                             guess_mean=None, guess_sigma=15.0, fit_gaus=True,
                             xlabel_hist=""):
    """
    For every bin along the scan axis (X) of a booked TH2, project onto Y, optionally
    fit a single gaussian over the FULL range, and plot every slice's histogram (+fit)
    into one grid figure.

    Returns: centers, means, sigmas, mean_errs, sigma_errs, entries, slice_data
    where slice_data = [(center, h1, fit_or_None), ...] for the slices actually plotted.
    """
    h2v = h2_lazy.GetValue() if hasattr(h2_lazy, "GetValue") else h2_lazy
    n_bins = h2v.GetNbinsX()

    centers, means, sigmas, mean_errs, sigma_errs, entries = [], [], [], [], [], []
    slice_data = []

    for ix in range(1, n_bins + 1):
        h1 = h2v.ProjectionY(f"py_{h2v.GetName()}_{ix}", ix, ix)
        if rebin and rebin > 1:
            h1.Rebin(rebin)

        e = h1.GetEntries()
        if e < min_entries:
            continue

        center = h2v.GetXaxis().GetBinCenter(ix)
        f = None
        mean, mean_err, sigma, sigma_err = h1.GetMean(), 0.0, h1.GetRMS(), 0.0

        if fit_gaus:
            m0 = guess_mean if guess_mean is not None else h1.GetMean()
            f = ROOT.TF1(f"f_{h2v.GetName()}_{ix}", "gaus", -fit_range_local, fit_range_local)
            f.SetParameters(h1.GetMaximum(), m0, max(h1.GetRMS(), guess_sigma))
            h1.Fit(f, "RQ0")
            mean, mean_err = f.GetParameter(1), f.GetParError(1)
            sigma, sigma_err = f.GetParameter(2), f.GetParError(2)

        centers.append(center)
        means.append(mean); mean_errs.append(mean_err)
        sigmas.append(sigma); sigma_errs.append(sigma_err)
        entries.append(e)
        slice_data.append((center, h1, f))

    if save is not None:
        _plot_slice_grid(slice_data, axis_label, save, xlabel_hist=xlabel_hist)

    return centers, means, sigmas, mean_errs, sigma_errs, entries, slice_data


def _plot_slice_grid(slice_data, axis_label, save, xlabel_hist=""):
    """Draw one subplot per slice (histogram + gaussian fit, if any), arranged in a
    roughly square grid, reusing plot_histo1d for each panel."""
    n = len(slice_data)
    if n == 0:
        print(f"WARNING: no slices to plot for '{save}' (all below min_entries).")
        return

    ncols = math.ceil(math.sqrt(n))
    nrows = math.ceil(n / ncols)

    fig, axs = plt.subplots(nrows, ncols, figsize=(3.2 * ncols, 2.6 * nrows), squeeze=False)
    axs_flat = axs.flatten()

    for ax, (center, h1, f) in zip(axs_flat, slice_data):
        plot_histo1d(h1, fit_func=f, ax=ax, style="fill", color="tab:blue",
                     xlabel=xlabel_hist, ylabel="")
        if f is not None:
            title = f"{axis_label} = {center:.2f}   $\\sigma$ = {f.GetParameter(2):.1f}"
        else:
            title = f"{axis_label} = {center:.2f}   N = {h1.GetEntries():.0f}"
        ax.set_title(title, fontsize=9)

    for ax in axs_flat[n:]:
        ax.axis("off")

    fig.tight_layout()
    fig.savefig(save, dpi=150)
    fig.savefig(save.replace(".pdf", ".png"), dpi=150)
    plt.close(fig)
    print(f"Saved slice grid ({n} panels, {nrows}x{ncols}) -> {save}")


# # x-scan: Deltatheta_x vs Tracks.d0Out_x, restricted to a rough y window around the crystal
# h2_x = book_scan_histogram(df_phys, "Tracks.d0Out_x", "Deltatheta_x", scan_min=-3, scan_max=4, n_bins=140,
#                             slice_var="Tracks.d0_y", slice_min=-6.0, slice_max=7.0)
# x_centers, x_means, x_sigmas, _, x_sigma_errs, _, x_slices = extract_and_plot_slices(
#     h2_x, axis_label="x", save=f"{PLOT_DIR}/slices_grid_x.pdf", rebin=3,
#     xlabel_hist=r"$\Delta\theta_x$ [$\mu$rad]")

# # y-scan: Deltatheta_y vs Tracks.d0Out_y, over the full illuminated range
# if file in [8430, 8431, 8650]:
#     scan_min_y, scan_max_y = -15, 9
# elif file in [8655, 8656]:
#     scan_min_y, scan_max_y = -15, 15
# else:
#     scan_min_y, scan_max_y = -15, 15

# h2_y = book_scan_histogram(df_phys, "Tracks.d0Out_y", "Deltatheta_y", scan_min=scan_min_y, scan_max=scan_max_y,
#                             n_bins=100, slice_var="Tracks.d0_x", slice_min=-1.0, slice_max=1.0)
# y_centers, y_means, y_sigmas, _, y_sigma_errs, _, y_slices = extract_and_plot_slices(
#     h2_y, axis_label="y", save=f"{PLOT_DIR}/slices_grid_y.pdf", rebin=3,
#     xlabel_hist=r"$\Delta\theta_y$ [$\mu$rad]")

# # sigma(x) / sigma(y) profile points, no edge fit overlay
# fig_x, ax_x = plt.subplots(figsize=(8, 6))
# ax_x.errorbar(x_centers, x_sigmas, yerr=x_sigma_errs, fmt='o', color='black', ecolor='gray', markersize=3)
# ax_x.set_xlabel("x [mm]"); ax_x.set_ylabel(r"$\sigma(\Delta\theta_x)$ [$\mu$rad]")
# ax_x.set_title("Local scattering width vs x")
# fig_x.tight_layout(); fig_x.savefig(f"{PLOT_DIR}/scattering_width_vs_x.png"); plt.close(fig_x)

# fig_y, ax_y = plt.subplots(figsize=(8, 6))
# ax_y.errorbar(y_centers, y_sigmas, yerr=y_sigma_errs, fmt='o', color='black', ecolor='gray', markersize=3)
# ax_y.set_xlabel("y [mm]"); ax_y.set_ylabel(r"$\sigma(\Delta\theta_y)$ [$\mu$rad]")
# ax_y.set_title("Local scattering width vs y")
# fig_y.tight_layout(); fig_y.savefig(f"{PLOT_DIR}/scattering_width_vs_y.png"); plt.close(fig_y)

# Definisci i due dataset da processare
datasets = {
    "filtered": df_filtered,
    "unfiltered": df_phys
}

for label, current_df in datasets.items():
    print(f"\nProcessing dataset: {label.upper()}")
    
    # x-scan: Deltatheta_x vs Tracks.d0Out_x, restricted to a rough y window around the crystal
    h2_x = book_scan_histogram(current_df, "Tracks.d0Out_x", "Deltatheta_x", scan_min=-3, scan_max=4, n_bins=140,
                                slice_var="Tracks.d0_y", slice_min=-6.0, slice_max=7.0)
    x_centers, x_means, x_sigmas, _, x_sigma_errs, _, x_slices = extract_and_plot_slices(
        h2_x, axis_label="x", save=f"{PLOT_DIR}/slices_grid_x_{label}.pdf", rebin=3,
        xlabel_hist=r"$\Delta\theta_x$ [$\mu$rad]")

    # y-scan: Deltatheta_y vs Tracks.d0Out_y, over the full illuminated range
    if file in [8430, 8431, 8650]:
        scan_min_y, scan_max_y = -15, 9
    elif file in [8655, 8656]:
        scan_min_y, scan_max_y = -15, 15
    else:
        scan_min_y, scan_max_y = -15, 15

    h2_y = book_scan_histogram(current_df, "Tracks.d0Out_y", "Deltatheta_y", scan_min=scan_min_y, scan_max=scan_max_y,
                                n_bins=100, slice_var="Tracks.d0_x", slice_min=-1.0, slice_max=1.0)
    y_centers, y_means, y_sigmas, _, y_sigma_errs, _, y_slices = extract_and_plot_slices(
        h2_y, axis_label="y", save=f"{PLOT_DIR}/slices_grid_y_{label}.pdf", rebin=3,
        xlabel_hist=r"$\Delta\theta_y$ [$\mu$rad]")

    # sigma(x) / sigma(y) profile points, no edge fit overlay
    fig_x, ax_x = plt.subplots(figsize=(8, 6))
    ax_x.errorbar(x_centers, x_sigmas, yerr=x_sigma_errs, fmt='o', color='black', ecolor='gray', markersize=3)
    ax_x.set_xlabel("x [mm]"); ax_x.set_ylabel(r"$\sigma(\Delta\theta_x)$ [$\mu$rad]")
    ax_x.set_title(f"Local scattering width vs x ({label})")
    fig_x.tight_layout(); fig_x.savefig(f"{PLOT_DIR}/scattering_width_vs_x_{label}.png"); plt.close(fig_x)

    fig_y, ax_y = plt.subplots(figsize=(8, 6))
    ax_y.errorbar(y_centers, y_sigmas, yerr=y_sigma_errs, fmt='o', color='black', ecolor='gray', markersize=3)
    ax_y.set_xlabel("y [mm]"); ax_y.set_ylabel(r"$\sigma(\Delta\theta_y)$ [$\mu$rad]")
    ax_y.set_title(f"Local scattering width vs y ({label})")
    fig_y.tight_layout(); fig_y.savefig(f"{PLOT_DIR}/scattering_width_vs_y_{label}.png"); plt.close(fig_y)

