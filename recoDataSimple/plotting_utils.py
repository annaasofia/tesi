"""
Typical usage inside one of your analysis scripts:

    import plotting_utils as pu

    # you already have a lazy RDataFrame result, e.g.:
    h_defl = df.Histo1D((...), "Deltatheta_x")

    # 1) Plot raw histogram as points with errors
    pu.plot_histo1d(h_defl, xlabel=r"$\\Delta\\theta_x$ [$\\mu$rad]", 
                    style="errorbar", save="deflection_x.pdf")

    # 2) Plot histogram as filled steps + existing TF1 fit overlaid
    pu.plot_histo1d(h_defl, fit_func=gaus_fit, style="fill",
                    xlabel=r"$\\Delta\\theta_x$ [$\\mu$rad]",
                    fit_label=f"gaus ($\\sigma$={gaus_fit.GetParameter(2):.1f})",
                    save="deflection_fit.pdf")

    # 3) Torsion map (TH2D)
    pu.plot_histo2d(h2_torsion_map, xlabel="x [mm]", ylabel="y [mm]",
                    zlabel=r"$\\theta_0$ [$\\mu$rad]", save="torsion_map.pdf")

    # 4) Scattering-width scan: TGraphErrors + TF1
    pu.plot_graph_with_fit(gr_y, f_y, xlabel="d0_y [mm]",
                           ylabel=r"$\\sigma(\\Delta\\theta_y)$ [$\\mu$rad]",
                           save="edge_scan_y.pdf")
"""


import matplotlib
matplotlib.use("Agg")
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator, LogLocator


# ------------------------------------------------------------------
# low-level: ROOT object -> numpy arrays
# ------------------------------------------------------------------

def _resolve(obj):
    """Materialize an RDataFrame lazy result (RResultPtr) if needed."""
    if hasattr(obj, "GetValue"):
        return obj.GetValue()
    return obj
 
 
def _style_axes(ax):
    """Shared cosmetic touch-up: minor ticks between the major ones,
    all ticks pointing inward, drawn on all four sides of the plot."""
    
    # Asse X: usa LogLocator se la scala è log, altrimenti AutoMinorLocator
    if ax.get_xscale() == 'log':
        ax.xaxis.set_minor_locator(LogLocator(base=10.0, subs=np.arange(2, 10)*0.1, numticks=100))
    else:
        ax.xaxis.set_minor_locator(AutoMinorLocator())
        
    # Asse Y: usa LogLocator se la scala è log, altrimenti AutoMinorLocator
    if ax.get_yscale() == 'log':
        ax.yaxis.set_minor_locator(LogLocator(base=10.0, subs=np.arange(2, 10)*0.1, numticks=100))
    else:
        ax.yaxis.set_minor_locator(AutoMinorLocator())

    ax.tick_params(which='both', direction='in', top=True, right=True)
    ax.tick_params(which='minor', length=3)
    ax.tick_params(which='major', length=5)


def th1_to_arrays(h):
    """TH1 -> (centers, contents, errors, edges)."""
    h = _resolve(h)
    n = h.GetNbinsX()
    centers = np.array([h.GetBinCenter(i) for i in range(1, n + 1)])
    contents = np.array([h.GetBinContent(i) for i in range(1, n + 1)])
    errors = np.array([h.GetBinError(i) for i in range(1, n + 1)])
    edges = np.array([h.GetBinLowEdge(i) for i in range(1, n + 2)])
    return centers, contents, errors, edges


def th2_to_arrays(h2):
    """TH2 -> (xedges, yedges, z) ready for pcolormesh (z has shape (ny, nx))."""
    h2 = _resolve(h2)
    nx, ny = h2.GetNbinsX(), h2.GetNbinsY()
    z = np.array([[h2.GetBinContent(ix, iy) for iy in range(1, ny + 1)]
                  for ix in range(1, nx + 1)])
    xedges = np.array([h2.GetXaxis().GetBinLowEdge(i) for i in range(1, nx + 2)])
    yedges = np.array([h2.GetYaxis().GetBinLowEdge(i) for i in range(1, ny + 2)])
    return xedges, yedges, z.T  # Transpose: pcolormesh wants (ny, nx)


def th2_errors_to_array(h2):
    """TH2 -> 2D array of bin errors, same orientation as th2_to_arrays' z."""
    h2 = _resolve(h2)
    nx, ny = h2.GetNbinsX(), h2.GetNbinsY()
    zerr = np.array([[h2.GetBinError(ix, iy) for iy in range(1, ny + 1)]
                     for ix in range(1, nx + 1)])
    return zerr.T


def th3_to_arrays(h3):
    """TH3 -> (xedges, yedges, zedges, values) with values.shape = (nx, ny, nz)."""
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
    """Project a TH3 onto one axis, then convert to arrays."""
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
    """TGraph or TGraphErrors -> (x, y, xerr, yerr). xerr/yerr are zero arrays for plain TGraph."""
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


def tf1_to_curve(f, x_min=None, x_max=None, n_points=800):
    """Evaluate a ROOT TF1 over a range -> (x, y) smooth curve for overlay."""
    if x_min is None: x_min = f.GetXmin()
    if x_max is None: x_max = f.GetXmax()
    x = np.linspace(x_min, x_max, n_points)
    y = np.array([f.Eval(xi) for xi in x])
    return x, y

def tf2_to_grid(f, x_range, y_range, nx=120, ny=120):
    """Evaluate a ROOT TF2 over a regular grid -> (xedges, yedges, z) ready for pcolormesh."""
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

# ------------------------------------------------------------------
# high-level plotting wrappers
# ------------------------------------------------------------------

_LOC_ANCHORS = {
    "upper right":  (0.85, 0.94, "right", "top"),
    "upper left":   (0.03, 0.94, "left",  "top"),
    "upper center": (0.50, 0.94, "center", "top"),
    "lower right":  (0.97, 0.06, "right", "bottom"),
    "lower left":   (0.03, 0.06, "left",  "bottom"),
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
            ha="center", va=va, multialignment="center", fontsize=fontsize,
            bbox=dict(boxstyle="round,pad=0.4", facecolor="white",
                    edgecolor="black", linewidth=0.8, alpha=0.9))
    return ax


def plot_histo1d(h, fit_func=None, fit_range=None, ax=None,
                 xlabel="", ylabel="Events", title="",
                 data_label="data", fit_label="fit",
                 info_text=None, info_loc="upper right",
                 style="errorbar", color="tab:blue", logy=False, save=None):
    """
    Plot a TH1 with an optional TF1 fit.
    style: "errorbar" (points+errors), "step" (outline), or "fill" (filled steps).
    """
    centers, contents, yerr, edges = th1_to_arrays(h)
    
    owns_axes = False
    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 6))
        owns_axes = True
    else:
        fig = ax.figure

    if style == "errorbar":
        xerr = (edges[1:] - edges[:-1]) / 2.0
        ax.errorbar(centers, contents, yerr=yerr, xerr=xerr, fmt='.', color=color,
                    ecolor=color, markersize=4, label=data_label)
    elif style == "step":
        ax.stairs(contents, edges, fill=False, color=color, lw=1.5, label=data_label)
    elif style == "fill":
        ax.stairs(contents, edges, fill=True, color=color, alpha=1.0, label=data_label)
    else:
        raise ValueError("style must be 'errorbar', 'step', or 'fill'")

    if fit_func is not None:
        fr = fit_range if fit_range is not None else (None, None)
        fx, fy = tf1_to_curve(fit_func, *fr)
        ax.plot(fx, fy, color='red', lw=1.8, label=None if info_text else fit_label)

    if logy:
        ax.set_yscale('log')

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    if title:
        ax.set_title(title)
    if fit_func is not None and not info_text:
        ax.legend(frameon=False, loc='upper right')
    if info_text is not None:
        add_info_box(ax, info_text, loc=info_loc)

    _style_axes(ax)

    if save:
        fig.tight_layout()
        fig.savefig(save, dpi=150)
        
        if save.endswith(".pdf"):
            fig.savefig(save.replace(".pdf", ".png"), dpi=150)
            
        if owns_axes:
            plt.close(fig)
            
    return ax


def plot_histo2d(h2, ax=None, xlabel="", ylabel="", zlabel="",
                 title="", cmap='viridis', save=None, colorbar=True):
    """Plot a TH2 (e.g. torsion map, efficiency map) as a pcolormesh."""
    xedges, yedges, z = th2_to_arrays(h2)

    owns_axes = False
    if ax is None:
        fig, ax = plt.subplots(figsize=(7, 6))
        owns_axes = True
    else:
        fig = ax.figure

    mesh = ax.pcolormesh(xedges, yedges, z, cmap=cmap, shading='flat')
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    if title:
        ax.set_title(title)
    if colorbar:
        cb = plt.colorbar(mesh, ax=ax)
        cb.set_label(zlabel)

    _style_axes(ax)

    if save:
        fig.tight_layout()
        fig.savefig(save, dpi=150)
        
        if save.endswith(".pdf"):
            fig.savefig(save.replace(".pdf", ".png"), dpi=150)
            
        if owns_axes:
            plt.close(fig)
            
    return ax


def plot_graph_with_fit(g, fit_func=None, fit_range=None, ax=None,
                        xlabel="", ylabel="", title="",
                        data_label="data", fit_label="fit",
                        info_text=None, info_loc="upper right",
                        connect_points=False, save=None):
    """Plot a TGraphErrors as points+errorbars, optionally overlay a TF1."""
    x, y, xerr, yerr = tgraph_to_arrays(g)

    owns_axes = False
    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 6))
        owns_axes = True
    else:
        fig = ax.figure

    fmt = 'o-' if connect_points else 'o'
    ax.errorbar(x, y, yerr=yerr, xerr=xerr if xerr.any() else None, 
                fmt=fmt, color='black', ecolor='gray', markersize=3, 
                elinewidth=0.8, capsize=0, label=data_label)

    if fit_func is not None:
        fr = fit_range if fit_range is not None else (None, None)
        fx, fy = tf1_to_curve(fit_func, *fr)
        ax.plot(fx, fy, color='red', lw=1.8, label=None if info_text else fit_label)

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    if title:
        ax.set_title(title)
    if fit_func is not None and not info_text:
        ax.legend(frameon=False)
    if info_text is not None:
        add_info_box(ax, info_text, loc=info_loc)

    _style_axes(ax)

    if save:
        fig.tight_layout()
        fig.savefig(save, dpi=150)
        
        if save.endswith(".pdf"):
            fig.savefig(save.replace(".pdf", ".png"), dpi=150)
            
        if owns_axes:
            plt.close(fig)
            
    return ax
 

def plot_histo3d_slice(h3, axis_fixed, bin_fixed, ax=None,
                       xlabel="", ylabel="", zlabel="",
                       title="", cmap='viridis', save=None):
    """Plot a 2D slice of a TH3 at a fixed bin along one axis ('x', 'y', or 'z')."""
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


def plot_tf2(f, x_range, y_range, nx=120, ny=120, ax=None,
             xlabel="", ylabel="", zlabel="", title="",
             cmap='viridis', save=None, colorbar=True):
    """Plot a ROOT TF2 as a smooth pcolormesh over the given range."""
    xedges, yedges, z = tf2_to_grid(f, x_range, y_range, nx=nx, ny=ny)

    owns_axes = False
    if ax is None:
        fig, ax = plt.subplots(figsize=(7, 6))
        owns_axes = True
    else:
        fig = ax.figure

    mesh = ax.pcolormesh(xedges, yedges, z, cmap=cmap, shading='auto')
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    if title:
        ax.set_title(title)
    if colorbar:
        cb = plt.colorbar(mesh, ax=ax)
        cb.set_label(zlabel)

    _style_axes(ax)

    if save:
        fig.tight_layout()
        fig.savefig(save, dpi=150)

        if save.endswith(".pdf"):
            fig.savefig(save.replace(".pdf", ".png"), dpi=150)
            
        if owns_axes:
            plt.close(fig)
            
    return ax