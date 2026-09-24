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
df_phys = df_phys.Filter(f"abs(thetaIn_x) < {theta_L}")
df_phys = df_phys.Define("Deltatheta_x", "(Tracks.thetaOut_x - Tracks.thetaIn_x) * 1e6") \
    .Define("Deltatheta_y", "(Tracks.thetaOut_y - Tracks.thetaIn_y) * 1e6") \
    .Define("DeltathetaErr_x", "sqrt(Tracks.thetaInErr_x * Tracks.thetaInErr_x + Tracks.thetaOutErr_x * Tracks.thetaOutErr_x) * 1e6") \
    .Define("DeltathetaErr_y", "sqrt(Tracks.thetaInErr_y * Tracks.thetaInErr_y + Tracks.thetaOutErr_y * Tracks.thetaOutErr_y) * 1e6")
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
    plt.close(fig)
    print(f"Saved slice grid ({n} panels, {nrows}x{ncols}) -> {save}")


# ===================== BASELINE: single gaussian per position slice =====================
# x-scan: Deltatheta_x vs Tracks.d0Out_x, restricted to a rough y window around the crystal
h2_x = book_scan_histogram(df_phys, "Tracks.d0Out_x", "Deltatheta_x", scan_min=-3, scan_max=4, n_bins=140,
                            slice_var="Tracks.d0_y", slice_min=-6.0, slice_max=7.0)
x_centers, x_means, x_sigmas, _, x_sigma_errs, _, x_slices = extract_and_plot_slices(
    h2_x, axis_label="x", save=f"{PLOT_DIR}/slices_grid_x.pdf", rebin=3,
    xlabel_hist=r"$\Delta\theta_x$ [$\mu$rad]")

# y-scan: Deltatheta_y vs Tracks.d0Out_y, over the full illuminated range
if file in [8430, 8431, 8650]:
    scan_min_y, scan_max_y = -15, 9
elif file in [8655, 8656]:
    scan_min_y, scan_max_y = -15, 15
else:
    scan_min_y, scan_max_y = -15, 15

h2_y = book_scan_histogram(df_phys, "Tracks.d0Out_y", "Deltatheta_y", scan_min=scan_min_y, scan_max=scan_max_y,
                            n_bins=100, slice_var="Tracks.d0_x", slice_min=-1.0, slice_max=1.0)
y_centers, y_means, y_sigmas, _, y_sigma_errs, _, y_slices = extract_and_plot_slices(
    h2_y, axis_label="y", save=f"{PLOT_DIR}/slices_grid_y.pdf", rebin=3,
    xlabel_hist=r"$\Delta\theta_y$ [$\mu$rad]")

# sigma(x) / sigma(y) profile points, no edge fit overlay
fig_x, ax_x = plt.subplots(figsize=(8, 6))
ax_x.errorbar(x_centers, x_sigmas, yerr=x_sigma_errs, fmt='o', color='black', ecolor='gray', markersize=3)
ax_x.set_xlabel("x [mm]"); ax_x.set_ylabel(r"$\sigma(\Delta\theta_x)$ [$\mu$rad]")
ax_x.set_title("Local scattering width vs x")
fig_x.tight_layout(); fig_x.savefig(f"{PLOT_DIR}/scattering_width_vs_x.png"); plt.close(fig_x)

fig_y, ax_y = plt.subplots(figsize=(8, 6))
ax_y.errorbar(y_centers, y_sigmas, yerr=y_sigma_errs, fmt='o', color='black', ecolor='gray', markersize=3)
ax_y.set_xlabel("y [mm]"); ax_y.set_ylabel(r"$\sigma(\Delta\theta_y)$ [$\mu$rad]")
ax_y.set_title("Local scattering width vs y")
fig_y.tight_layout(); fig_y.savefig(f"{PLOT_DIR}/scattering_width_vs_y.png"); plt.close(fig_y)


# ===================== IDEA 1: subtract the well-fitted (VR) peak, refit the residual =====================
# Only for x, only inside the crystal: the low-side gaussian (the "red" one, guess mean
# ~-45 urad) always fits cleanly. Fit ONLY that one, subtract it bin-by-bin from the
# full histogram, and refit a fresh gaussian on what's left -- hopefully a clean
# zero-centered peak. y is left untouched for now.

def subtract_secondary_peak_and_refit(h2_lazy, x_range, cut,
                                       guess_secondary_mean, guess_secondary_sigma=15.0,
                                       guess_primary_mean=0.0, guess_primary_sigma=15.0,
                                       fit_range_local=fit_range, min_entries=minimum_entries,
                                       rebin=1, save_prefix=None, clip_negative=True,
                                       xlabel_hist=""):
    """
    For slices of h2_lazy (X axis) whose center falls within x_range:
      1) fit a single gaussian ONLY on the well-behaved secondary peak, restricted to
         its own side of `cut` (the side where guess_secondary_mean sits)
      2) subtract that gaussian (evaluated bin-by-bin) from the full histogram
      3) fit a fresh, independent gaussian on the residual

    Produces two grid figures (if save_prefix is given): the original histograms with
    the secondary-peak fit overlaid, and the residual histograms with the refitted
    primary peak overlaid.

    Returns: centers, primary_means, primary_mean_errs, primary_sigmas, primary_sigma_errs
    """
    h2v = h2_lazy.GetValue() if hasattr(h2_lazy, "GetValue") else h2_lazy
    n_bins = h2v.GetNbinsX()
    secondary_is_low = guess_secondary_mean < cut

    centers, prim_means, prim_mean_errs, prim_sigmas, prim_sigma_errs = [], [], [], [], []
    orig_slices, resid_slices = [], []

    for ix in range(1, n_bins + 1):
        center = h2v.GetXaxis().GetBinCenter(ix)
        if not (x_range[0] <= center <= x_range[1]):
            continue

        h1 = h2v.ProjectionY(f"py_sub_{h2v.GetName()}_{ix}", ix, ix)
        if rebin and rebin > 1:
            h1.Rebin(rebin)
        if h1.GetEntries() < min_entries:
            continue

        # 1) fit the well-behaved secondary peak on its own side of the cut only
        lo, hi = (-fit_range_local, cut) if secondary_is_low else (cut, fit_range_local)
        f_sec = ROOT.TF1(f"f_sec_{ix}", "gaus", lo, hi)
        f_sec.SetParameters(h1.GetMaximum(), guess_secondary_mean, guess_secondary_sigma)
        h1.Fit(f_sec, "RQ0")

        # 2) subtract it bin-by-bin over the FULL range (its tail can leak across `cut`)
        h_sub = h1.Clone(f"{h1.GetName()}_sub")
        for i in range(1, h_sub.GetNbinsX() + 1):
            x_bin = h_sub.GetBinCenter(i)
            new_val = h1.GetBinContent(i) - f_sec.Eval(x_bin)
            if clip_negative and new_val < 0:
                new_val = 0.0
            h_sub.SetBinContent(i, new_val)

        # 3) refit the residual, independently, over the full range
        f_prim = ROOT.TF1(f"f_prim_{ix}", "gaus", -fit_range_local, fit_range_local)
        f_prim.SetParameters(h_sub.GetMaximum(), guess_primary_mean, guess_primary_sigma)
        h_sub.Fit(f_prim, "RQ0")

        centers.append(center)
        prim_means.append(f_prim.GetParameter(1)); prim_mean_errs.append(f_prim.GetParError(1))
        prim_sigmas.append(f_prim.GetParameter(2)); prim_sigma_errs.append(f_prim.GetParError(2))
        orig_slices.append((center, h1, f_sec))
        resid_slices.append((center, h_sub, f_prim))

    if save_prefix is not None:
        _plot_slice_grid(orig_slices, "x", f"{save_prefix}_original_and_secondary_fit.pdf",
                          xlabel_hist=xlabel_hist)
        _plot_slice_grid(resid_slices, "x", f"{save_prefix}_residual_and_primary_fit.pdf",
                          xlabel_hist=xlabel_hist + " (secondary peak subtracted)")

    return centers, prim_means, prim_mean_errs, prim_sigmas, prim_sigma_errs


sub_x_centers, sub_x_means, sub_x_mean_errs, sub_x_sigmas, sub_x_sigma_errs = subtract_secondary_peak_and_refit(
    h2_x, x_range=(-1.1, 1.0), cut=-25.0,
    guess_secondary_mean=-45.0, guess_secondary_sigma=15.0,
    guess_primary_mean=0.0, guess_primary_sigma=15.0,
    rebin=3, save_prefix=f"{PLOT_DIR}/x_vr_subtraction",
    xlabel_hist=r"$\Delta\theta_x$ [$\mu$rad]")

# for name, sigma, sigma_err in zip(sub_x_centers, sub_x_sigmas, sub_x_sigma_errs):
#     print(f"x = {name:.2f} mm  ->  sigma (post-subtraction) = {sigma:.1f} +/- {sigma_err:.1f} urad")


# ===================== IDEA 2: thetaIn_x distribution per Deltatheta_x bin, one fixed (x,y) slice =====================
# Pick one narrow spatial bin (x, y) where the two peaks overlap, and for each bin of
# Deltatheta_x look at the distribution of the *incoming* angle thetaIn_x of the
# particles that ended up there. This tells us whether the secondary peak comes from a
# distinct population in incoming angle.

x_slice_center, x_slice_halfwidth = -0.52, 0.03   # mm -- matches the ~0.05 mm x-scan bin width
y_slice_center, y_slice_halfwidth = 1.20, 0.12    # mm -- matches the ~0.24 mm y-scan bin width

df_pos_slice = df_phys.Filter(
    f"Tracks.d0Out_x > {x_slice_center - x_slice_halfwidth} && "
    f"Tracks.d0Out_x < {x_slice_center + x_slice_halfwidth} && "
    f"Tracks.d0Out_y > {y_slice_center - y_slice_halfwidth} && "
    f"Tracks.d0Out_y < {y_slice_center + y_slice_halfwidth}")

n_events_slice = df_pos_slice.Count().GetValue()
print("=" * 50)
print(f"Events in the (x={x_slice_center}, y={y_slice_center}) mm slice: {n_events_slice}")

thetain_min, thetain_max = -60, 60   # urad -- adjust to where thetaIn_x actually lives
n_dtheta_bins_2d = 20                # coarser binning along Deltatheta_x, for a readable grid

h2_dtheta_thetain = book_scan_histogram(
    df_pos_slice, "Deltatheta_x", "thetaIn_x",
    scan_min=-fit_range, scan_max=fit_range, n_bins=n_dtheta_bins_2d,
    delta_min=thetain_min, delta_max=thetain_max, n_delta_bins=60)


def _safe_tag(v):
    s = f"{v:+.2f}"
    return s.replace("+", "p").replace("-", "m").replace(".", "")


tag_slice = f"x{_safe_tag(x_slice_center)}_y{_safe_tag(y_slice_center)}"

# min_entries here must be much lower: we've already cut down to one narrow spatial bin,
# then split further by Deltatheta_x bin, so each panel has far fewer events.
extract_and_plot_slices(
    h2_dtheta_thetain, axis_label=r"$\Delta\theta_x$", rebin=1, min_entries=20,
    fit_gaus=True, guess_mean=0.0, guess_sigma=10.0,
    xlabel_hist=r"$\theta_{in,x}$ [$\mu$rad]",
    save=f"{PLOT_DIR}/thetaIn_per_deltatheta_bin_{tag_slice}.pdf")

print(f"\nAll matplotlib plots saved under {PLOT_DIR}/")