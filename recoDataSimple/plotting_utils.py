"""
plotting_utils.py

Convert ROOT objects (TH1, TH2, TH3, TGraph/TGraphErrors, TF1) into numpy
arrays and matplotlib plots, without leaving your existing PyROOT/RDataFrame
analysis pipeline.

Typical usage inside one of your analysis scripts (e.g. compute_channeling.py):

    import plotting_utils as pu

    # you already have a lazy RDataFrame result, e.g.:
    h_defl = df.Histo1D((...), "Deltatheta_x")

    # 1) plot the raw histogram
    pu.plot_histo1d(h_defl, xlabel=r"$\\Delta\\theta_x$ [$\\mu$rad]",
                     ylabel="Events", save="deflection_x.pdf")

    # 2) plot histogram + an existing ROOT TF1 fit (e.g. gaus_fit from
    #    channeling_efficiency()) overlaid
    pu.plot_histo1d(h_defl, fit_func=gaus_fit,
                     xlabel=r"$\\Delta\\theta_x$ [$\\mu$rad]", ylabel="Events",
                     fit_range=(fit_min, fit_max), save="deflection_fit.pdf")

    # 3) torsion map (TH2D)
    pu.plot_histo2d(h2_torsion_map, xlabel="x [mm]", ylabel="y [mm]",
                     zlabel=r"$\\theta_0$ [$\\mu$rad]", save="torsion_map.pdf")

    # 4) scattering-width edge scan: TGraphErrors + TF1 step function
    pu.plot_graph_with_fit(gr_y, f_y, xlabel="d0_y [mm]",
                            ylabel=r"$\\sigma(\\Delta\\theta_y)$ [$\\mu$rad]",
                            save="edge_scan_y.pdf")

All functions accept either the "lazy" RResultPtr you get straight out of
RDataFrame (Histo1D/Histo2D/Histo3D return these) or an already-materialized
ROOT object (TH1, TH2, TGraph, ...) — GetValue() is called automatically
when needed.
"""

import numpy as np
import matplotlib.pyplot as plt

# ------------------------------------------------------------------
# low-level: ROOT object -> numpy arrays
# ------------------------------------------------------------------

def _resolve(obj):
    """Materialize an RDataFrame lazy result (RResultPtr) if needed."""
    if hasattr(obj, "GetValue"):
        return obj.GetValue()
    return obj


def th1_to_arrays(h):
    """TH1 -> (bin_centers, contents, errors, bin_half_width)."""
    h = _resolve(h)
    n = h.GetNbinsX()
    x = np.array([h.GetBinCenter(i) for i in range(1, n + 1)])
    y = np.array([h.GetBinContent(i) for i in range(1, n + 1)])
    yerr = np.array([h.GetBinError(i) for i in range(1, n + 1)])
    xerr = h.GetBinWidth(1) / 2.0
    return x, y, yerr, xerr


def th2_to_arrays(h2):
    """TH2 -> (xedges, yedges, z) ready for pcolormesh (z has shape (ny, nx))."""
    h2 = _resolve(h2)
    nx, ny = h2.GetNbinsX(), h2.GetNbinsY()
    z = np.array([[h2.GetBinContent(ix, iy) for iy in range(1, ny + 1)]
                  for ix in range(1, nx + 1)])
    xedges = np.array([h2.GetXaxis().GetBinLowEdge(i) for i in range(1, nx + 2)])
    yedges = np.array([h2.GetYaxis().GetBinLowEdge(i) for i in range(1, ny + 2)])
    return xedges, yedges, z.T  # transpose: pcolormesh wants (ny, nx)


def th2_errors_to_array(h2):
    """TH2 -> 2D array of bin errors, same orientation as th2_to_arrays' z."""
    h2 = _resolve(h2)
    nx, ny = h2.GetNbinsX(), h2.GetNbinsY()
    zerr = np.array([[h2.GetBinError(ix, iy) for iy in range(1, ny + 1)]
                     for ix in range(1, nx + 1)])
    return zerr.T


def th3_to_arrays(h3):
    """TH3 -> (xedges, yedges, zedges, values) with values.shape = (nx, ny, nz).

    Useful mainly to project/slice afterwards in numpy (e.g. sum over one
    axis) since matplotlib has no native 3D-histogram volume plot worth
    using for physics data. For a specific 2D slice use
    plot_histo3d_slice() below instead of pulling the whole cube.
    """
    h3 = _resolve(h3)
    nx, ny, nz = h3.GetNbinsX(), h3.GetNbinsY(), h3.GetNbinsZ()
    xedges = np.array([h3.GetXaxis().GetBinLowEdge(i) for i in range(1, nx + 2)])
    yedges = np.array([h3.GetYaxis().GetBinLowEdge(i) for i in range(1, ny + 2)])
    zedges = np.array([h3.GetZaxis().GetBinLowEdge(i) for i in range(1, nz + 2)])
    values = np.empty((nx, ny, nz))
    for ix in range(1, nx + 1):
        for iy in range(1, ny + 1):
            for iz in range(1, nz + 1):
                values[ix - 1, iy - 1, iz - 1] = h3.GetBinContent(ix, iy, iz)
    return xedges, yedges, zedges, values


def project_th3_to_1d(h3, axis, bin_range_other1=None, bin_range_other2=None):
    """Project a TH3 onto one axis (ROOT's own ProjectionX/Y/Z), then convert
    to arrays. axis: 'x', 'y', or 'z'. bin_range_other1/2 restrict the other
    two axes to (bin_lo, bin_hi) if given (bin numbers, 1-indexed).
    Mirrors the pattern already used in torsion_map() for h3_all/h3_chan.
    """
    h3 = _resolve(h3)
    proj_map = {"x": h3.ProjectionX, "y": h3.ProjectionY, "z": h3.ProjectionZ}
    args = []
    if bin_range_other1 is not None:
        args += list(bin_range_other1)
    if bin_range_other2 is not None:
        args += list(bin_range_other2)
    h1 = proj_map[axis](f"{h3.GetName()}_proj_{axis}", *args)
    return th1_to_arrays(h1)


def tgraph_to_arrays(g):
    """TGraph or TGraphErrors -> (x, y, xerr, yerr). xerr/yerr are zero
    arrays if the graph has no errors (plain TGraph)."""
    g = _resolve(g)
    n = g.GetN()
    x = np.array([g.GetPointX(i) for i in range(n)])
    y = np.array([g.GetPointY(i) for i in range(n)])
    try:
        xerr = np.array([g.GetErrorX(i) for i in range(n)])
        yerr = np.array([g.GetErrorY(i) for i in range(n)])
    except AttributeError:
        xerr = np.zeros(n)
        yerr = np.zeros(n)
    return x, y, xerr, yerr


def tf2_to_grid(f, x_range, y_range, nx=120, ny=120):
    """Evaluate a ROOT TF2 (e.g. the fitted torsion surface torsion_fit_2d)
    over a regular grid -> (xedges, yedges, z) ready for pcolormesh, the
    matplotlib equivalent of f.Draw('surf2')."""
    x_min, x_max = x_range
    y_min, y_max = y_range
    xcenters = np.linspace(x_min, x_max, nx)
    ycenters = np.linspace(y_min, y_max, ny)
    z = np.array([[f.Eval(xc, yc) for xc in xcenters] for yc in ycenters])

    dx = (x_max - x_min) / (nx - 1) if nx > 1 else 0.0
    dy = (y_max - y_min) / (ny - 1) if ny > 1 else 0.0
    xedges = np.linspace(x_min - dx / 2.0, x_max + dx / 2.0, nx + 1)
    yedges = np.linspace(y_min - dy / 2.0, y_max + dy / 2.0, ny + 1)
    return xedges, yedges, z


def plot_tf2(f, x_range, y_range, nx=120, ny=120, ax=None,
             xlabel="", ylabel="", zlabel="", title="",
             cmap='viridis', save=None, colorbar=True):
    """Plot a ROOT TF2 as a smooth pcolormesh over the given range (replaces
    f.Draw('surf2') / f.Draw('colz')), e.g. the continuous fitted torsion
    surface torsion_fit_2d."""
    xedges, yedges, z = tf2_to_grid(f, x_range, y_range, nx=nx, ny=ny)

    if ax is None:
        fig, ax = plt.subplots(figsize=(7, 6))

    mesh = ax.pcolormesh(xedges, yedges, z, cmap=cmap, shading='auto')
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    if title:
        ax.set_title(title)
    if colorbar:
        cb = plt.colorbar(mesh, ax=ax)
        cb.set_label(zlabel)

    if save:
        plt.savefig(save, bbox_inches='tight')
    return ax


def tf1_to_curve(f, x_min=None, x_max=None, n_points=500):
    """Evaluate a ROOT TF1 over a range -> (x, y) smooth curve for overlay
    on a matplotlib plot (replaces f.Draw('SAME'))."""
    if x_min is None:
        x_min = f.GetXmin()
    if x_max is None:
        x_max = f.GetXmax()
    x = np.linspace(x_min, x_max, n_points)
    y = np.array([f.Eval(xi) for xi in x])
    return x, y


# ------------------------------------------------------------------
# high-level plotting wrappers
# ------------------------------------------------------------------

_LOC_ANCHORS = {
    # (x, y) in axes fraction, plus matching ha/va — mimics ROOT's
    # TPaveText(x1, y1, x2, y2, "NDC") corners, but with the text
    # actually centered inside its own box (ha='center'), not left-aligned.
    "upper right": (0.97, 0.94, "right", "top"),
    "upper left":  (0.03, 0.94, "left",  "top"),
    "upper center": (0.50, 0.94, "center", "top"),
    "lower right": (0.97, 0.06, "right", "bottom"),
    "lower left":  (0.03, 0.06, "left",  "bottom"),
    "lower center": (0.50, 0.06, "center", "bottom"),
}


def add_info_box(ax, text, loc="upper right", fontsize=10):
    """Draw a text box with properly CENTERED text (multi-line included),
    the matplotlib equivalent of a ROOT TPaveText with SetTextAlign(22).

    text: a single string, or a list of strings (one per line, e.g.
    [f"eff = {eff:.1f} +/- {err:.1f} %", f"peak = {mean:.1f} urad"]).
    loc: one of 'upper right', 'upper left', 'upper center',
         'lower right', 'lower left', 'lower center'.

    Unlike passing fit_label into ax.legend() (which left-aligns text next
    to the line swatch), this centers every line of text horizontally
    within its own bounding box, matching how you used TPaveText/AddText.
    """
    if isinstance(text, (list, tuple)):
        text = "\n".join(text)

    x, y, ha, va = _LOC_ANCHORS[loc]
    ax.text(x, y, text, transform=ax.transAxes,
             ha="center", va=va,
             multialignment="center", fontsize=fontsize,
             bbox=dict(boxstyle="round,pad=0.4", facecolor="white",
                       edgecolor="black", linewidth=0.8))
    return ax


def plot_histo1d(h, fit_func=None, fit_range=None, ax=None,
                  xlabel="", ylabel="Events", title="",
                  data_label="data", fit_label="fit",
                  info_text=None, info_loc="upper right",
                  logy=False, save=None, show_errors=True):
    """Plot a TH1 as points+errorbars, optionally overlay a TF1 fit curve.

    fit_range: (x_min, x_max) to restrict where the fit curve is drawn;
    defaults to the TF1's own range if omitted.
    """
    x, y, yerr, xerr = th1_to_arrays(h)

    if ax is None:
        fig, ax = plt.subplots(figsize=(7, 5))

    if show_errors:
        ax.errorbar(x, y, yerr=yerr, xerr=xerr, fmt='.', color='darkorange',
                     ecolor='darkorange', markersize=4, label=data_label)
    else:
        ax.step(x, y, where='mid', color='darkorange', label=data_label)

    if fit_func is not None:
        fr = fit_range if fit_range is not None else (None, None)
        fx, fy = tf1_to_curve(fit_func, *fr)
        # only put the fit in the legend if there's no dedicated centered
        # info box (avoids the same info appearing twice)
        ax.plot(fx, fy, color='crimson', lw=1.8,
                 label=None if info_text else fit_label)

    if logy:
        ax.set_yscale('log')

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    if title:
        ax.set_title(title)
    if fit_func is not None and not info_text:
        ax.legend(frameon=False)
    if info_text is not None:
        add_info_box(ax, info_text, loc=info_loc)

    if save:
        plt.savefig(save, bbox_inches='tight')
    return ax


def plot_histo2d(h2, ax=None, xlabel="", ylabel="", zlabel="",
                  title="", cmap='viridis', save=None, colorbar=True):
    """Plot a TH2 (e.g. torsion map, efficiency map) as a pcolormesh."""
    xedges, yedges, z = th2_to_arrays(h2)

    if ax is None:
        fig, ax = plt.subplots(figsize=(7, 6))

    mesh = ax.pcolormesh(xedges, yedges, z, cmap=cmap, shading='flat')
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    if title:
        ax.set_title(title)
    if colorbar:
        cb = plt.colorbar(mesh, ax=ax)
        cb.set_label(zlabel)

    if save:
        plt.savefig(save, bbox_inches='tight')
    return ax


def plot_graph_with_fit(g, fit_func=None, fit_range=None, ax=None,
                         xlabel="", ylabel="", title="",
                         data_label="data", fit_label="fit",
                         info_text=None, info_loc="upper right",
                         connect_points=False, save=None):
    """Plot a TGraphErrors as points+errorbars, optionally overlay a TF1
    (e.g. the erf-based step function from fit_step_edges() in
    compute_edges2.py, or a gaussian efficiency-scan fit).

    connect_points: draw a thin line through the data points too, matching
    ROOT's gr.Draw("APL") style used e.g. for sliding-window scans (no fit
    involved, just points connected in order).
    """
    x, y, xerr, yerr = tgraph_to_arrays(g)

    if ax is None:
        fig, ax = plt.subplots(figsize=(7, 5))

    fmt = 'o-' if connect_points else 'o'
    ax.errorbar(x, y, yerr=yerr, xerr=xerr, fmt=fmt, color='steelblue',
                 ecolor='steelblue', markersize=4, lw=1.2, label=data_label)

    if fit_func is not None:
        fr = fit_range if fit_range is not None else (None, None)
        fx, fy = tf1_to_curve(fit_func, *fr)
        ax.plot(fx, fy, color='crimson', lw=1.8,
                 label=None if info_text else fit_label)

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    if title:
        ax.set_title(title)
    if fit_func is not None and not info_text:
        ax.legend(frameon=False)
    if info_text is not None:
        add_info_box(ax, info_text, loc=info_loc)

    if save:
        plt.savefig(save, bbox_inches='tight')
    return ax


def plot_histo3d_slice(h3, axis_fixed, bin_fixed, ax=None,
                        xlabel="", ylabel="", zlabel="",
                        title="", cmap='viridis', save=None):
    """Plot a 2D slice of a TH3 at a fixed bin along one axis
    (e.g. a fixed d0_y bin of the thetaIn_x/d0_x/d0_y cube used in
    torsion_map()). axis_fixed: 'x', 'y', or 'z'."""
    h3 = _resolve(h3)

    if axis_fixed == "x":
        h3.GetXaxis().SetRange(bin_fixed, bin_fixed)
        h2 = h3.Project3D("zy")
    elif axis_fixed == "y":
        h3.GetYaxis().SetRange(bin_fixed, bin_fixed)
        h2 = h3.Project3D("zx")
    elif axis_fixed == "z":
        h3.GetZaxis().SetRange(bin_fixed, bin_fixed)
        h2 = h3.Project3D("yx")
    else:
        raise ValueError("axis_fixed must be 'x', 'y', or 'z'")

    return plot_histo2d(h2, ax=ax, xlabel=xlabel, ylabel=ylabel,
                          zlabel=zlabel, title=title, cmap=cmap, save=save)