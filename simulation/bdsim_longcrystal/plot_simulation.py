import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from scipy.optimize import curve_fit
import os

plt.rcParams.update({
    'font.size': 16,          # Dimensione generale del font
    'axes.titlesize': 18,     # Dimensione dei titoli dei singoli assi/subplot
    'axes.labelsize': 16,     # Dimensione delle etichette (asse x, y)
    'xtick.labelsize': 14,    # Dimensione dei numeri sull'asse x
    'ytick.labelsize': 14,    # Dimensione dei numeri sull'asse y
    'legend.fontsize': 14,    # Dimensione del testo nella legenda
    'figure.titlesize': 20    # Dimensione del suptitle (titolo globale della figura)
})

viridis_yellow = "#fde725"
viridis_green = "#5ec962"
viridis_teal = "#21918c"
viridis_blue = "#3b528b"
viridis_purple = "#440154"

ROUTINE = sys.argv[1] if len(sys.argv) > 1 else 'bdsim'

if ROUTINE == 'bdsim':
    filename = './runs/cry2_merged6_slim.npz'
    output_folder = 'plots6'
    dark_color = "#982916"
    light_color = "#DA7762"
elif ROUTINE == 'xcoll':
    # filename = '../xcoll_longcrystal/runs/cry2_20261021_16.npz'
    filename = '../xcoll_longcrystal/runs/cry2_20261112_10.npz'
    output_folder = '../xcoll_longcrystal/plots'
    dark_color = "#68238B"
    light_color = "#A76ED1"
else:
    raise ValueError(f"Unknown routine: {ROUTINE}")

data = np.load(filename)
os.makedirs(output_folder, exist_ok=True)

# reconstruct the arrays from the saved data
x_in = data['x_in']; theta_x_in = data['px_in']
y_in = data['y_in']; theta_y_in = data['py_in']
# zeta_in = data['zeta_in']; delta_in = data['delta_in']; ptau_in = data['ptau_in']

x_out=data['x_out']; theta_x_out=data['px_out']
y_out=data['y_out']; theta_y_out=data['py_out']
# zeta_out=data['zeta_out']; delta_out=data['delta_out']; ptau_out=data['ptau_out']

particle_id_in = data['particle_id_in']
# state_in = data['state_in']#; pdg_id_in = data['pdg_id_in']
particle_id_out = data['particle_id_out']; state_out = data['state_out']
#  pdg_id_out = data['pdg_id_out']; at_element_out = data['at_element_out']

# single values
npart = int(data['npart'])
crystal_length = float(data['crystal_length']); print(f"crystal_length = {crystal_length*1e3:.2f} mm")
crystal_x = float(data['crystal_x']); print(f"crystal_x = {crystal_x*1e3:.2f} mm")
crystal_y = float(data['crystal_y']); print(f"crystal_y = {crystal_y*1e3:.2f} mm")
bending_angle = float(data['bending_angle']); print(f"bending_angle = {bending_angle*1e3:.2f} mrad")
crystal_material = str(data['crystal_material'])
p0c = float(data['p0c'])
mass0 = float(data['mass0'])
label = str(data['label'])

theta_L = (1-180/5.7*bending_angle/(crystal_length*100))*np.sqrt(32/(180*10**9))
print(f"theta_L = {theta_L*1e6:.2f} urad")  # rad Lindhard critical angle for Si(110) at 180 GeV

if ROUTINE == 'bdsim':
    if 'state_in' in data.files:
        valid_in = data['state_in'] == 1
        x_in, y_in = x_in[valid_in], y_in[valid_in]
        theta_x_in, theta_y_in = theta_x_in[valid_in], theta_y_in[valid_in]
        particle_id_in = particle_id_in[valid_in]
    # altrimenti: file slim, già filtrato a monte
    # print(f"valid in: {particle_id_in.size} / npart = {npart}")

order_in = np.argsort(particle_id_in)
sorted_ids = particle_id_in[order_in]
assert np.unique(particle_id_in).size == particle_id_in.size

# perse
lost_mask = state_out == -333
idx_lost = particle_id_out[lost_mask]
src_lost = order_in[np.searchsorted(sorted_ids, idx_lost)]
x_in_lost, y_in_lost = x_in[src_lost], y_in[src_lost]
assert np.all(particle_id_in[src_lost] == idx_lost)

# sopravvissute
if ROUTINE == 'bdsim':
    survived = state_out == 1
if ROUTINE == 'xcoll':
    finite_out = (np.isfinite(x_out) & np.isfinite(y_out) &
              np.isfinite(theta_x_out) & np.isfinite(theta_y_out))
    # print(f"Non-finite outputs among state==1: {((state_out == 1) & ~finite_out).sum()}")
    survived = (state_out == 1) & finite_out
idx_out = particle_id_out[survived]
idx = np.searchsorted(sorted_ids, idx_out)
src = order_in[idx]                       # indice diretto, niente copie
assert np.all(particle_id_in[src] == idx_out)

x_in_survived, y_in_survived = x_in[src], y_in[src]
theta_x_in_survived, theta_y_in_survived = theta_x_in[src], theta_y_in[src]
x_out_survived, y_out_survived = x_out[survived], y_out[survived]
theta_x_out_survived = theta_x_out[survived]
theta_y_out_survived = theta_y_out[survived]
dtheta_x = theta_x_out_survived - theta_x_in_survived
dtheta_y = theta_y_out_survived - theta_y_in_survived

survived_data = {
    'x_in': x_in_survived,
    'y_in': y_in_survived,
    'x_out': x_out_survived,
    'y_out': y_out_survived,
    'theta_x_in': theta_x_in_survived,
    'theta_y_in': theta_y_in_survived,
    'theta_x_out': theta_x_out_survived,
    'theta_y_out': theta_y_out_survived,
    'dtheta_x': dtheta_x,
    'dtheta_y': dtheta_y,
    'particle_id': idx_out,
}

print(f"Total particles: {len(x_in)}")

surv_in = np.zeros(particle_id_in.size, dtype=bool)
surv_in[src] = True
lost_in = (~surv_in)
# controllo: le non sopravvissute sono tutte state == -333?
print("Not survived:", lost_in.sum(), "| state -333:", lost_mask.sum())
# print(f"Survived particles: {np.sum(survived)} = {np.sum(survived)/len(x_in)*100:.2f}%")

# ---------- DETECTOR SMEAR (dataset parallelo, gli originali restano intatti) ----------
rng_smear = np.random.default_rng(seed=83)

def highland(p_GeV, x_over_X0, beta=1.0, z=1):
    if x_over_X0 <= 0:
        return 0.0
    return (13.6e-3 / (beta * p_GeV)) * z * np.sqrt(x_over_X0) \
           * (1 + 0.038 * np.log(x_over_X0 * z**2 / beta**2))

X0_air      = 303.9      # m
L_up        = 0.0        # m aria tra ultimo piano a monte e cristallo  
L_down      = 0.0        # m aria tra cristallo e primo piano a valle   
sigma_track = 9e-6       # rad                                          
p_GeV       = p0c / 1e9
sigma_position = 1e-4  # m, smearing posizione (x,y) dei piani di tracking

s_in  = np.hypot(sigma_track, highland(p_GeV, L_up   / X0_air))
s_out = np.hypot(sigma_track, highland(p_GeV, L_down / X0_air))

# 1) θ_in misurato per TUTTE le particelle valide (un'estrazione sola)
theta_x_in_meas_all = theta_x_in + rng_smear.standard_normal(theta_x_in.size, dtype=np.float32) * s_in
theta_y_in_meas_all = theta_y_in + rng_smear.standard_normal(theta_y_in.size, dtype=np.float32) * s_in

x_in_meas_all = x_in + rng_smear.standard_normal(x_in.size, dtype=np.float32) * sigma_position
y_in_meas_all = y_in + rng_smear.standard_normal(y_in.size, dtype=np.float32) * sigma_position

# 2) stesso valore propagato alle sopravvissute (stessa mappatura già usata sopra)
theta_x_in_meas_surv = theta_x_in_meas_all[src]
theta_y_in_meas_surv = theta_y_in_meas_all[src]
x_in_meas_surv = x_in_meas_all[src]
y_in_meas_surv = y_in_meas_all[src]

# 3) θ_out misurato (estrazioni indipendenti, solo sulle sopravvissute)
theta_x_out_meas_surv = theta_x_out_survived + rng_smear.standard_normal(theta_x_out_survived.size, dtype=np.float32) * s_out
theta_y_out_meas_surv = theta_y_out_survived + rng_smear.standard_normal(theta_y_out_survived.size, dtype=np.float32) * s_out

survived_data_smeared = dict(survived_data)
survived_data_smeared.update({
    'x_in':        x_in_meas_surv,
    'y_in':        y_in_meas_surv,
    'theta_x_in':  theta_x_in_meas_surv,
    'theta_y_in':  theta_y_in_meas_surv,
    'theta_x_out': theta_x_out_meas_surv,
    'theta_y_out': theta_y_out_meas_surv,
    'dtheta_x':    theta_x_out_meas_surv - theta_x_in_meas_surv,
    'dtheta_y':    theta_y_out_meas_surv - theta_y_in_meas_surv,
})
# print(f"Smear: s_in = {s_in*1e6:.2f} urad, s_out = {s_out*1e6:.2f} urad")

# CHANNELED
d  = (crystal_x/2 - survived_data['x_in']) * 1e3          # mm, distanza dal bordo destro
dx = (survived_data['x_out'] - survived_data['x_in']) * 1e3
ok = (np.abs(survived_data['y_in']) < crystal_y/2) & (d > 0) & (d < 0.6)

fig, ax = plt.subplots(figsize=(6, 5))
ax.hist2d(d[ok], dx[ok], bins=[120, 150], cmap='viridis', norm=LogNorm())
ax.axhline(crystal_length*bending_angle/2*1e3, color='red', ls='--', label='Δx = θ_b L/2')
ax.plot([0, 0.6], [0, 0.6], 'w:', label='x_out = +edge')
ax.set_xlabel('d dal bordo destro [mm]'); ax.set_ylabel('x_out − x_in [mm]')
ax.legend(); fig.savefig(f'{output_folder}/edge_dx_test.png', dpi=200)

def _grid_edges(crystal_x, crystal_y, nx_bins, ny_bins, distance_unit='mm'):
    ds = 1e3 if distance_unit == 'mm' else 1.0
    xe = np.linspace(-crystal_x / 2, crystal_x / 2 - EDGE_CUT_RIGHT, nx_bins + 1) * ds
    ye = np.linspace(-crystal_y / 2, crystal_y / 2, ny_bins + 1) * ds
    return xe, ye, ds

def efficiency_maps(x_den, y_den, x_num, y_num, dtheta_num, mu, sigma, crystal_x, crystal_y,
                    nx_bins=6, ny_bins=7, n_sigma=3, min_den=50, min_num=10, tag='smeared',
                    distance_unit='mm', angle_unit='urad'):
    """
    x_den, y_den : posizioni di impatto di TUTTE le particelle 'attempted'
                   (geometric hit + Lindhard), incluse le perse -> denominatore
    x_num, y_num, dtheta_num : particelle sopravvissute dopo gli stessi tagli (final_data*)
    mu, sigma    : picco di canalizzazione dal fit globale (in rad)
    Canalizzata  = |dtheta - mu| < n_sigma * sigma  (conteggio, niente fit per bin)
    Griglia: nx_bins lungo x, ny_bins lungo y.
    """
    xe, ye, ds = _grid_edges(crystal_x, crystal_y, nx_bins, ny_bins, distance_unit)
    ascale = 1e6 if angle_unit == 'urad' else 1.0
    unit_str = r'\mu rad' if angle_unit == 'urad' else 'rad'
 
    in_win = np.abs(dtheta_num - mu) < n_sigma * sigma
    xc, yc, dc = x_num[in_win] * ds, y_num[in_win] * ds, dtheta_num[in_win]
 
    den, _, _ = np.histogram2d(x_den * ds, y_den * ds, bins=[xe, ye])
    num, _, _ = np.histogram2d(xc, yc, bins=[xe, ye])
 
    with np.errstate(divide='ignore', invalid='ignore'):
        eff = np.where(den >= min_den, num / den, np.nan)
 
    # peak per bin (mediana dei canalizzati)
    ix = np.digitize(xc, xe) - 1
    iy = np.digitize(yc, ye) - 1
    ok = (ix >= 0) & (ix < nx_bins) & (iy >= 0) & (iy < ny_bins)
    peak = np.full((nx_bins, ny_bins), np.nan)
    for i in range(nx_bins):
        for j in range(ny_bins):
            sel = ok & (ix == i) & (iy == j)
            if sel.sum() >= min_num:
                peak[i, j] = np.median(dc[sel]) * ascale
 
    # --- figura 2x1: efficiency | deflection peak ---
    fig, axs = plt.subplots(1, 2, figsize=(13, 6))
    panels = [
        (eff * 100, f'Channelling efficiency [%] - {ROUTINE.upper()}', axs[0]),
        (peak, rf'Deflection peak [${unit_str}$]', axs[1]),
    ]
    for M, title, ax in panels:
        pc = ax.pcolormesh(xe, ye, M.T, cmap='viridis')
        fig.colorbar(pc, ax=ax)
        ax.set_xlabel(f'impact x [{distance_unit}]')
        ax.set_ylabel(f'impact y [{distance_unit}]')
        ax.set_title(f'{title}')
        for i in range(nx_bins):
            for j in range(ny_bins):
                if np.isfinite(M[i, j]):
                    txt_color = 'black' if pc.norm(M[i, j]) > 0.5 else 'white'
                    ax.text(0.5 * (xe[i] + xe[i + 1]), 0.5 * (ye[j] + ye[j + 1]),
                            f'{M[i, j]:.1f}', ha='center', va='center', fontsize=14, color=txt_color)
    fig.tight_layout()
    plt.savefig(f'{output_folder}/efficiency_maps_{tag}.png', dpi=200)
    plt.savefig(f'{output_folder}/efficiency_maps_{tag}.pdf')
    plt.close(fig)
 
    print(f"[{tag}] map: {np.isfinite(eff).sum()}/{nx_bins * ny_bins} bins with den >= {min_den}; "
          f"count-based global eff (±{n_sigma}σ) = {num.sum() / den.sum() * 100:.1f}%")
    return eff, peak, den, num


def eff_profile(var_den, var_num_all, dtheta_num, mu, sigma, edges, scale=1.0, n_sigma=3, min_den=50):
    """Efficienza 1D vs una variabile qualunque (x, y, theta_x_in...). Errori binomiali."""
    in_win = np.abs(dtheta_num - mu) < n_sigma * sigma
    den, _ = np.histogram(var_den * scale, bins=edges)
    num, _ = np.histogram(var_num_all[in_win] * scale, bins=edges)
    with np.errstate(divide='ignore', invalid='ignore'):
        eff = np.where(den >= min_den, num / den, np.nan)
        err = np.sqrt(eff * (1 - eff) / den)
    return 0.5 * (edges[1:] + edges[:-1]), eff, err
 
def profiles_plot(x_den, y_den, th_den, fd, mu, sigma, crystal_x, crystal_y, theta_L, tag='smeared', n_bins=10):
    """fd = final_data (o final_data_s): dict con x_in, y_in, theta_x_in, dtheta_x"""
    fig, axs = plt.subplots(1, 3, figsize=(16, 4.5))
    specs = [
        (x_den, fd['x_in'], np.linspace(-crystal_x / 2, crystal_x / 2 - EDGE_CUT_RIGHT, n_bins + 1) * 1e3, 1e3, 'impact x [mm]'),
        (y_den, fd['y_in'], np.linspace(-crystal_y / 2, crystal_y / 2, n_bins + 1) * 1e3, 1e3, 'impact y [mm]'),
        (th_den, fd['theta_x_in'], np.linspace(-0.5 * theta_L, 0.5 * theta_L, n_bins + 1) * 1e6, 1e6, r'$\theta_{x,in}$ [$\mu$rad]'),
    ]
    for ax, (vd, vn, edges, sc, xl) in zip(axs, specs):
        c, e, er = eff_profile(vd, vn, fd['dtheta_x'], mu, sigma, edges, scale=sc)
        ax.errorbar(c, e * 100, yerr=er * 100, fmt='o', color="black", capsize=3)
        ax.set_xlabel(xl); ax.set_ylabel('Efficiency [%]'); ax.grid(alpha=0.3)
        if vd is th_den:
            ax.axvline(-0.5 * theta_L * 1e6, color='red', lw=1, ls='--', label=r'$\pm 1/2 \theta_L$')
            ax.axvline(0.5 * theta_L * 1e6, color='red', lw=1, ls='--')
            ax.legend()
    fig.suptitle(f'Efficiency profiles')
    fig.tight_layout()
    plt.savefig(f'{output_folder}/efficiency_profiles_{tag}.png', dpi=200)
    plt.savefig(f'{output_folder}/efficiency_profiles_{tag}.pdf')
    plt.close(fig)

def beam_distribution_plot(x_in, y_in, theta_x_in, theta_y_in, label=None, distance_unit='mm', angle_unit='urad'):
    scale = 1e6 if angle_unit == 'urad' else 1.0
    unit_str = r'\mu rad' if angle_unit == 'urad' else 'rad'

    distance_scale = 1e3 if distance_unit == 'mm' else 1.0

    # x_lim = (-3, 3) if distance_unit == 'mm' else (x_in.min() * distance_scale, x_in.max() * distance_scale)
    # y_lim = (-11, 11) if distance_unit == 'mm' else (y_in.min() * distance_scale, y_in.max() * distance_scale)
    x_lim = (x_in.min() * distance_scale, x_in.max() * distance_scale)
    y_lim = (y_in.min() * distance_scale, y_in.max() * distance_scale)
    theta_x_lim = (-150, 150) if angle_unit == 'urad' else (theta_x_in.min() * scale, theta_x_in.max() * scale)
    theta_y_lim = (-200, 200) if angle_unit == 'urad' else (theta_y_in.min() * scale, theta_y_in.max() * scale)
    # theta_x_lim = (theta_x_in.min() * scale, theta_x_in.max() * scale)
    # theta_y_lim = (theta_y_in.min() * scale, theta_y_in.max() * scale)

    fig, graph = plt.subplots(1, 4, figsize=(20, 6))
    # 1. Grafico (x, y) - Beam impact position
    h0 = graph[0].hist2d(x_in * distance_scale, y_in * distance_scale, bins=100, cmap='viridis', norm=LogNorm(), rasterized=True)
    # fig.colorbar(h0[3], ax=graph[0], label='Counts')
    graph[0].set_xlabel(f'impact x [{distance_unit}]')
    graph[0].set_ylabel(f'impact y [{distance_unit}]')
    graph[0].set_title(f'Beam Impact Position (x, y)')
    graph[0].set_xlim(x_lim)
    graph[0].set_ylim(y_lim)
    graph[0].grid(alpha=0.3)

    h00 = graph[1].hist2d(theta_x_in * scale, theta_y_in * scale, bins=100, cmap='viridis', norm=LogNorm(), rasterized=True)
    # fig.colorbar(h00[3], ax=graph[1], label='Counts')
    graph[1].set_xlabel(rf'$\theta_x$ [${unit_str}$]')
    graph[1].set_ylabel(rf'$\theta_y$ [${unit_str}$]')
    graph[1].set_title(f'Angular Distribution')
    graph[1].set_xlim(theta_x_lim)
    graph[1].set_ylim(theta_y_lim)
    graph[1].grid(alpha=0.3)

    h1 = graph[2].hist2d(x_in * distance_scale, theta_x_in * scale, bins=100, cmap='viridis', norm=LogNorm(), rasterized=True)
    # fig.colorbar(h1[3], ax=graph[2], label='Counts') 
    graph[2].set_xlabel(f'impact x [{distance_unit}]')
    graph[2].set_ylabel(rf'$\theta_x$ [${unit_str}$]')
    graph[2].set_title(f'Horizontal Phase Space (x, $θ_x$)')
    graph[2].set_xlim(x_lim)
    graph[2].set_ylim(theta_x_lim)
    graph[2].grid(alpha=0.3)

    h2 = graph[3].hist2d(y_in * distance_scale, theta_y_in * scale, bins=100, cmap='viridis', norm=LogNorm(), rasterized=True)
    # fig.colorbar(h2[3], ax=graph[3], label='Counts')
    graph[3].set_xlabel(f'impact y [{distance_unit}]')
    graph[3].set_ylabel(rf'$\theta_y$ [${unit_str}$]')
    graph[3].set_title(f'Vertical Phase Space (y, $θ_y$)')
    graph[3].set_xlim(y_lim)
    graph[3].set_ylim(theta_y_lim)
    graph[3].grid(alpha=0.3)

    fig.tight_layout()
    # label = label[2:]
    plt.savefig(f'{output_folder}/impact_distribution_{label}.png', dpi=300)
    plt.savefig(f'{output_folder}/impact_distribution_{label}.pdf')
    plt.close(fig)

def angular_scan_plot(theta_x_in, theta_x_out, theta_L, popt=None, fit_bin_width=None, label=None, angle_unit='urad'):

    scale = 1e6 if angle_unit == 'urad' else 1.0
    unit_str = r'\mu rad' if angle_unit == 'urad' else 'rad'

    dtheta_x = theta_x_out - theta_x_in

    fig, graph = plt.subplots(1, 3, figsize=(18, 5))

    x_range = (-7, 7) if label == '' else (-150, 150)
    h = graph[0].hist2d(theta_x_in * scale, dtheta_x * scale, bins=300, range=[x_range, [-750, 6200]], cmap='viridis', norm=LogNorm(), rasterized=True)
    cb = fig.colorbar(h[3], ax=graph[0])
    cb.set_label('Counts')
    graph[0].set_ylim(-750, 6200)
    # graph[0].axhline(0, color='gray', lw=0.8, ls='--')
    graph[0].axvline(-0.5 * theta_L * scale, color='red', lw=0.8, ls='--', label=r'$\pm 1/2 \theta_L$')
    graph[0].axvline(0.5 * theta_L * scale, color='red', lw=0.8, ls='--')
    graph[0].set_xlabel(rf'$\theta_{{in}}$ [${unit_str}$]')
    graph[0].set_ylabel(rf'$\Delta\theta_x = \theta_{{out}} - \theta_{{in}}$ [${unit_str}$]')
    graph[0].set_title(f'Angular scan')
    graph[0].grid(alpha=0.3)
    graph[0].legend()

    counts, bins, _ = graph[1].hist(dtheta_x * scale, bins=200, range=(-750, 6200), color=dark_color, rasterized=True)
    counts2, bins2, _ = graph[2].hist(dtheta_x * scale, bins=200, range=(5950, 6050), color=dark_color, rasterized=True)
    plot_bin_width = bins[1] - bins[0]
    plot_bin_width2 = bins2[1] - bins2[0]

    graph[1].set_xlim(-750, 6200)
    graph[1].set_yscale('log')
    graph[1].set_ylim(500, counts.max()*1.2)
    graph[1].set_xlabel(rf'$\Delta\theta_x$ [${unit_str}$]')
    graph[1].set_ylabel('Counts')
    graph[1].set_title(rf'$\Delta\theta_x$ distribution - {ROUTINE.upper()}')
    graph[1].grid(alpha=0.3)

    graph[2].set_xlim(5950, 6050)
    graph[2].set_ylim(0, counts2.max()*1.2)
    graph[2].set_xlabel(rf'$\Delta\theta_x$ [${unit_str}$]')
    graph[2].set_ylabel('Counts')
    graph[2].set_title(rf'$\Delta\theta_x$ distribution (Channelled peak zoom)')
    graph[2].grid(alpha=0.3)

    if popt is not None and fit_bin_width is not None:
        A, mu, sigma = popt

        A_scaled = A * (plot_bin_width / (fit_bin_width * scale))
        A_scaled2 = A * (plot_bin_width2 / (fit_bin_width * scale))
        # Attenzione: mu e sigma dal fit sono in radianti, mentre il grafico è in urad!
        # Dobbiamo moltiplicarli per 'scale'
        mu_scaled = mu * scale
        sigma_scaled = sigma * scale

        x_curve_wide = np.linspace(dtheta_x.min() * scale, dtheta_x.max() * scale, 2000)
        x_curve_zoom = np.linspace(5950, 6050, 2000)

        # graph[1].plot(x_curve_wide, gaussian(x_curve_wide, A_scaled, mu_scaled, sigma_scaled), color='red', lw=2)
        graph[2].plot(x_curve_zoom, gaussian(x_curve_zoom, A_scaled2, mu_scaled, sigma_scaled), color='red', lw=2)


    fig.tight_layout()
    if label == '(survived particles)': label = 'survived'
    if label == '': label = 'final'
    plt.savefig(f'{output_folder}/angular_scan_{label}.png', dpi=300)
    plt.savefig(f'{output_folder}/angular_scan_{label}.pdf')
    plt.close(fig)

def smear_comparison_plot(dtheta_x, dtheta_x_smeared, label=None, angle_unit='urad'):

    scale = 1e6 if angle_unit == 'urad' else 1.0
    unit_str = r'\mu rad' if angle_unit == 'urad' else 'rad'

    if label[-1] == 'y' or label[-1] == 'x':
        label_ax = label[-1]
    else:
        label_ax = 'x'

    fig, graph = plt.subplots(1, 2, figsize=(12, 5))

    counts_smeared, bins_smeared, _ = graph[0].hist(dtheta_x_smeared * scale, bins=200, range=(-750, 6200), color=dark_color, label='smeared', rasterized=True)
    counts, bins, _ = graph[0].hist(dtheta_x * scale, bins=200, range=(-750, 6200), color=light_color, alpha=0.7, label='original', rasterized=True)

    if label[-1] == 'y':
        graph[0].set_xlim(-250, 250)
    # graph[0].set_xlim(-750, 6200)
    graph[0].set_yscale('log')
    graph[0].set_ylim(500, max(counts.max(), counts_smeared.max())*1.2)
    graph[0].set_xlabel(rf'$\Delta\theta_{label_ax}$ [${unit_str}$]')
    graph[0].set_ylabel('Counts')
    graph[0].set_title(rf'$\Delta\theta_{label_ax}$ distribution - {ROUTINE.upper()}')
    graph[0].grid(alpha=0.3)
    graph[0].legend()

    graph[1].hist(dtheta_x_smeared * scale, bins=200, range=(5950, 6050), color=dark_color, label='smeared', rasterized=True)
    graph[1].hist(dtheta_x * scale, bins=200, range=(5950, 6050), color=light_color, alpha=0.7, label='original', rasterized=True)

    if label[-1] == 'y':
        graph[1].set_xlim(-60, 60)
    # graph[1].set_xlim(5950, 6050)
    # graph[1].set_ylim(0, max(counts.max(), counts_smeared.max())*1.2)
    graph[1].set_xlabel(rf'$\Delta\theta_{label_ax}$ [${unit_str}$]')
    graph[1].set_ylabel('Counts')
    graph[1].set_title(rf'$\Delta\theta_{label_ax}$ distribution (zoom)')
    graph[1].grid(alpha=0.3)
    graph[1].legend()

    fig.tight_layout()
    plt.savefig(f'{output_folder}/smear_comparison_{label}.png', dpi=300)
    plt.savefig(f'{output_folder}/smear_comparison_{label}.pdf')
    plt.close(fig)

def smear_comparison_plot2(x, x_smeared, y, y_smeared, angle_unit='urad'):
    scale = 1e6 if angle_unit == 'urad' else 1.0
    unit_str = r'\mu rad' if angle_unit == 'urad' else 'rad'

    fig, graph = plt.subplots(1, 2, figsize=(12, 5))

    counts_smeared, _, _ = graph[0].hist(x_smeared * scale, bins=200, range=(-150, 150), color=dark_color, label='smeared')
    counts, _, _ = graph[0].hist(x * scale, bins=200, range=(-150, 150), color=light_color, alpha=0.7, label='original')

    # graph[0].set_ylim(500, max(counts.max(), counts_smeared.max())*1.2)
    graph[0].set_xlabel(rf'$\theta_{{in,x}}$ [${unit_str}$]')
    graph[0].set_ylabel('Counts')
    graph[0].set_title(rf'Incoming $\theta_{{x}}$ distribution')
    graph[0].grid(alpha=0.3)
    graph[0].legend()

    graph[1].hist(y_smeared * scale, bins=200, range=(-150, 150), color=dark_color, label='smeared')
    graph[1].hist(y * scale, bins=200, range=(-150, 150), color=light_color, alpha=0.7, label='original')

    graph[1].set_xlabel(rf'$\theta_{{in,y}}$ [${unit_str}$]')
    graph[1].set_ylabel('Counts')
    graph[1].set_title(rf'Incoming $\theta_{{y}}$ distribution')
    graph[1].grid(alpha=0.3)
    graph[1].legend()

    fig.tight_layout()
    plt.savefig(f'{output_folder}/smear_comparison_angle.png', dpi=300)
    plt.savefig(f'{output_folder}/smear_comparison_angle.pdf')
    plt.close(fig)

def ch_footprint_plot(x_in, y_in, x_out, y_out, dtheta_x, crystal_x, crystal_y, tol_frac=0.2, bending_angle=bending_angle, label = '', distance_unit='mm'):
    distance_scale = 1e3 if distance_unit == 'mm' else 1.0

    channeled_mask = np.abs(dtheta_x - bending_angle) < (tol_frac * bending_angle)
    x_in_channeled = x_in[channeled_mask]
    y_in_channeled = y_in[channeled_mask]
    x_out_channeled = x_out[channeled_mask]
    y_out_channeled = y_out[channeled_mask]
    shift = crystal_length * bending_angle/2 # 0.222 mm

    fig = plt.figure(figsize=(20, 11), layout='constrained')
    subfigs = fig.subfigures(2, 1, hspace=0.01)

    configs = [
        {
            'subfig': subfigs[0],
            'title': 'All Particles',
            'data': [x_in, y_in, x_out, y_out],
            'color': viridis_teal,
        },
        {
            'subfig': subfigs[1],
            'title': 'Channeled Particles',
            'data': [x_in_channeled, y_in_channeled, x_out_channeled, y_out_channeled],
            'color': viridis_teal,
        }
    ]

    edge_x = crystal_x / 2 * distance_scale
    edge_y = crystal_y / 2 * distance_scale

    for cfg in configs:
        sf = cfg['subfig']
        # Titolo principale di riga (|titolo1| / |titolo2|)
        sf.suptitle(cfg['title'], fontsize=12, y=1.0)
        axes = sf.subplots(1, 4)

        # Sottotitoli raggruppati per colonne (|subtitle1| sopra plot 1-2, |subtitle2| sopra plot 3-4)
        axes[0].annotate('Upstream', xy=(1.0, 1.05), xycoords='axes fraction',
                         ha='center', va='bottom', fontsize=12)
        axes[2].annotate('Downstream', xy=(1.0, 1.05), xycoords='axes fraction',
                         ha='center', va='bottom', fontsize=12)

        labels = [
            (f'impact x [{distance_unit}]', edge_x),
            (f'impact y [{distance_unit}]', edge_y),
            (f'outgoing x [{distance_unit}]', edge_x),
            (f'outgoing y [{distance_unit}]', edge_y),
        ]

        for ax, vals, (xlabel, edge) in zip(axes, cfg['data'], labels):
            ax.hist(vals * distance_scale, bins=100, color=cfg['color'], histtype='stepfilled')
            ax.axvline(-edge, color='red', ls='--', label='crystal edges')
            ax.axvline(edge, color='red', ls='--')
            if cfg['title'] == 'Channeled Particles' and xlabel.startswith('outgoing x'):
                ax.axvline(-edge + shift * distance_scale, color='red', ls='--')
                ax.axvline(edge + shift * distance_scale, color='red', ls='--')
                ax.axvline(-edge, color='orange', ls='--', label='crystal edges no bend')
                ax.axvline(edge, color='orange', ls='--')
            ax.set_xlabel(xlabel)
            ax.set_ylabel('Counts')
            ax.legend(loc='upper right')

    plt.savefig(f'{output_folder}/ch_footprint_{label}.png', dpi=300, bbox_inches='tight')
    plt.savefig(f'{output_folder}/ch_footprint_{label}.pdf', bbox_inches='tight')
    plt.close(fig)

def mcs_edge_scan(pos_in, dtheta, crystal_edge, plane_label='x', dtheta_label='x', n_slices=50, distance_unit='mm', 
                  angle_unit='urad', fit_window=None, label='', min_per_bin=50, other_pos_in=None, other_edge=None):

    scale = 1e6 if angle_unit == 'urad' else 1.0
    unit_str = r'\mu rad' if angle_unit == 'urad' else 'rad'
    distance_scale = 1e3 if distance_unit == 'mm' else 1.0

    pos_in_scaled = pos_in * distance_scale
    dtheta_scaled = dtheta * scale

    # cut particles that are outside the crystal in the OTHER plane
    if other_pos_in is not None and other_edge is not None:
        within_other = np.abs(other_pos_in) < (other_edge / 2)
        n_before = len(pos_in_scaled)
        pos_in_scaled = pos_in_scaled[within_other]
        dtheta_scaled = dtheta_scaled[within_other]
        print(f"mcs_edge_scan ({plane_label}): kept {within_other.sum()}/{n_before} particles "
              f"within the orthogonal aperture (±{other_edge/2*distance_scale:.2f} {distance_unit})")

    # full histogram of x (or y) positions, sliced into n_slices equal-width bins
    counts_pos, bin_edges = np.histogram(pos_in_scaled, bins=n_slices)
    bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
    edge_scaled = crystal_edge / 2 * distance_scale

    sigmas, sigma_errs, valid_pos, n_used = [], [], [], []
    skipped_low_stats, skipped_fit_fail = 0, 0
    records = []  # for the example plot

    for i in range(n_slices):
        mask = (pos_in_scaled >= bin_edges[i]) & (pos_in_scaled < bin_edges[i + 1])
        dtheta_slice = dtheta_scaled[mask]

        if fit_window is not None:
            dtheta_slice = dtheta_slice[np.abs(dtheta_slice) < fit_window]

        if len(dtheta_slice) < min_per_bin:
            skipped_low_stats += 1
            continue

        counts, edges = np.histogram(dtheta_slice, bins=50)
        centers = 0.5 * (edges[:-1] + edges[1:])
        A0 = counts.max()
        mu0 = np.mean(dtheta_slice)
        sigma0 = np.std(dtheta_slice)

        popt = None
        try:
            popt, pcov = curve_fit(gaussian, centers, counts, p0=[A0, mu0, sigma0], maxfev=5000)
            sig = abs(popt[2])
            err = np.sqrt(pcov[2, 2])
        except RuntimeError:
            skipped_fit_fail += 1
            sig, err = sigma0, 0.0

        sigmas.append(sig)
        sigma_errs.append(err)
        valid_pos.append(bin_centers[i])
        n_used.append(len(dtheta_slice))
        records.append(dict(center=bin_centers[i], counts=counts, edges=edges, centers_h=centers,
                            popt=popt, n=len(dtheta_slice)))

    print(f"mcs_edge_scan ({plane_label}): {len(valid_pos)}/{n_slices} bins fit "
          f"({skipped_low_stats} skipped for low stats, {skipped_fit_fail} fell back to raw std)")

    valid_pos = np.array(valid_pos)
    sigmas = np.array(sigmas)
    sigma_errs = np.array(sigma_errs)
    n_used = np.array(n_used)

    # --- main sigma(position) plot ---
    fig, ax = plt.subplots(figsize=(8, 5))
    fell_back = sigma_errs == 0.0
    ax.errorbar(valid_pos[~fell_back], sigmas[~fell_back], yerr=sigma_errs[~fell_back],
                fmt='o', color='royalblue', markersize=4, capsize=3,
                label=rf'$\sigma$ of $\Delta\theta_{dtheta_label}$ (fit)')
    if fell_back.any():
        ax.scatter(valid_pos[fell_back], sigmas[fell_back], marker='x', color='gray',
                   s=30, label='raw std (fit failed)', zorder=3)

    ax.axvline(-edge_scaled, color='red', ls='--', label='Geometrical Crystal Edges')
    ax.axvline(edge_scaled, color='red', ls='--')
    ax.set_xlabel(f'impact {plane_label} [{distance_unit}]')
    ax.set_ylabel(rf'MCS width $\sigma(\Delta\theta_{dtheta_label})$ [${unit_str}$]')
    title = f'Material Edge Transition {label}'
    ax.set_title(title)
    ax.grid(alpha=0.3)
    ax.legend()
    fig.tight_layout()
    plt.savefig(f'{output_folder}/mcs_edge_scan_{plane_label}.png', dpi=300)
    plt.savefig(f'{output_folder}/mcs_edge_scan_{plane_label}.pdf', dpi=300)

    # --- example: one slice inside, one outside the crystal ---
    inside_recs = [r for r in records if abs(r['center']) <= edge_scaled]
    outside_recs = [r for r in records if abs(r['center']) > edge_scaled]

    if inside_recs and outside_recs:
        example_inside = min(inside_recs, key=lambda r: abs(r['center']))
        example_outside = max(outside_recs, key=lambda r: abs(r['center']))

        fig2, axes2 = plt.subplots(1, 2, figsize=(12, 5))
        for ax2, rec, subtitle in zip(axes2, [example_inside, example_outside],
                                      ['Inside crystal', 'Outside crystal']):
            width = rec['edges'][1] - rec['edges'][0]
            ax2.bar(rec['centers_h'], rec['counts'], width=width, color='darkorange', alpha=0.7)
            if rec['popt'] is not None:
                x_curve = np.linspace(rec['edges'][0], rec['edges'][-1], 300)
                ax2.plot(x_curve, gaussian(x_curve, *rec['popt']), 'r-', lw=2, label='Gaussian fit')
                ax2.legend(loc='upper center')
            ax2.set_xlabel(rf'$\Delta\theta_{dtheta_label}$ [${unit_str}$]')
            ax2.set_ylabel('Counts')
            ax2.set_title(f"{subtitle} — pos = {rec['center']:.2f} {distance_unit}, n = {rec['n']}")
            ax2.grid(alpha=0.3)
        fig2.suptitle(f'Example slice fits {label}')
        fig2.tight_layout()
        plt.savefig(f'{output_folder}/mcs_edge_scan_{plane_label}_example_fits.png', dpi=300)
        plt.savefig(f'{output_folder}/mcs_edge_scan_{plane_label}_example_fits.pdf', dpi=300)

    return valid_pos, sigmas, sigma_errs, n_used

def apply_selection(data_dict, cut):
    # apply same boolean mask to the arrrays of the dictionary
    return {k: v[cut] for k, v in data_dict.items()}

def gaussian(x, A, mu, sigma):
    return A * np.exp(-0.5 * ((x - mu) / sigma) ** 2)

def channeling_efficiency(theta_x_in, dtheta_x, bending_angle=bending_angle, tol_frac=0.4,
                          angle_unit='urad', weights=None, bins=7250, verbose=True):
    """
    Come prima, ma:
      - fit pesato Poisson (sigma = sqrt(counts), absolute_sigma=True)
      - pcov gonfiata di max(1, chi2/ndf) se il fit è povero
      - errore su N_ch propagato da pcov con la correlazione A-sigma
    weights/bins servono solo al bootstrap.
    Ritorna: n_channeled, n_channeled_err, popt, pcov, bin_width
    """
    scale = 1e6 if angle_unit == 'urad' else 1.0
 
    counts, edges = np.histogram(dtheta_x, bins=bins, weights=weights)
    centers = 0.5 * (edges[:-1] + edges[1:])
    bin_width = edges[1] - edges[0]
 
    search_window = 0.5 * abs(bending_angle)
    search_mask = np.abs(centers - bending_angle) < search_window
 
    if search_mask.sum() > 0:
        peak_x = centers[search_mask][np.argmax(counts[search_mask])]
    else:
        peak_x = bending_angle
 
    window = 0.5 * abs(bending_angle)
    mask = (centers >= (5990 / scale)) & (centers < (peak_x + 50 / scale))
    while mask.sum() < 5 and window < (dtheta_x.max() - dtheta_x.min()):
        window *= 2
        mask = np.abs(centers - bending_angle) < window
 
    x_fit, y_fit = centers[mask], counts[mask]
    A0 = y_fit.max() if len(y_fit) else 1.0
    sigma0 = 12e-6  # rad
 
    sig_y = np.sqrt(np.maximum(y_fit, 1.0))   # errori Poisson per bin
    popt, pcov = curve_fit(gaussian, x_fit, y_fit, p0=[A0, bending_angle, sigma0],
                           sigma=sig_y, absolute_sigma=True, maxfev=10000)
 
    chi2 = np.sum(((y_fit - gaussian(x_fit, *popt)) / sig_y) ** 2)
    ndf = len(x_fit) - 3
    chi2_ndf = chi2 / ndf if ndf > 0 else np.nan
    if np.isfinite(chi2_ndf) and chi2_ndf > 1:
        pcov = pcov * chi2_ndf
 
    A, mu, sigma_signed = popt
    sigma = abs(sigma_signed)
    perr = np.sqrt(np.diag(pcov))
    mu_err, sigma_err = perr[1], perr[2]
 
    n_channeled = A * sigma * np.sqrt(2 * np.pi) / bin_width
 
    # propagazione: N = A |sigma| sqrt(2pi) / bw, con covarianza A-sigma
    J = np.array([sigma * np.sqrt(2 * np.pi) / bin_width,
                  np.sign(sigma_signed) * A * np.sqrt(2 * np.pi) / bin_width])
    C = pcov[np.ix_([0, 2], [0, 2])]
    n_channeled_err = np.sqrt(J @ C @ J)
 
    if verbose:
        print('-' * 50)
        print(f"Channeling peak = {mu * scale:.2f} ± {mu_err * scale:.2f} urad | "
              f"Sigma = {sigma * scale:.2f} ± {sigma_err * scale:.2f} | "
              f"chi2/ndf = {chi2_ndf:.2f} ({len(x_fit)} bins)")
        print(f"N_ch = {n_channeled:.0f} ± {n_channeled_err:.0f} (fit)")
 
    return n_channeled, n_channeled_err, popt, pcov, bin_width

def eff_3sigma(fd, popt, n_tot, n_sigma=3):
    mu, sig = popt[1], abs(popt[2])
    n = np.sum(np.abs(fd['dtheta_x'] - mu) < n_sigma * sig)
    eff = n / n_tot
    return eff, np.sqrt(eff * (1 - eff) / n_tot)

def efficiency_bootstrap(dtheta_final, n_tot, n_boot=200, seed=1234, verbose=False):
    """
    Bootstrap Poissoniano: ogni sopravvissuta pesa Poisson(1), le perse contribuiscono
    solo al denominatore (Poisson(n_tot - n_surv)). Il fit è rifatto ad ogni replica con
    gli stessi bin del fit nominale. Cattura fit + conteggio + correlazione num/den.
    dtheta_final : dtheta_x di final_data (o final_data_s)
    n_tot        : n_tot_correct (o n_tot_correct_s)
    """
    rng = np.random.default_rng(seed)
    _, edges = np.histogram(dtheta_final, bins=7250)   # stessi bin del fit nominale
    n_surv = len(dtheta_final)
    effs = []
    for _ in range(n_boot):
        w = rng.poisson(1.0, n_surv).astype(float)
        n_tot_b = w.sum() + rng.poisson(max(n_tot - n_surv, 0))
        try:
            n_ch_b = channeling_efficiency(dtheta_final, dtheta_final, weights=w,
                                           bins=edges, verbose=False)[0]
        except RuntimeError:
            continue
        effs.append(n_ch_b / n_tot_b)
    effs = np.array(effs)
    if verbose:
        print(f"bootstrap: {len(effs)}/{n_boot} fit ok | eff = {effs.mean()*100:.2f}% | "
              f"std = {effs.std(ddof=1)*100:.3f}%")
    return effs.std(ddof=1), effs

def dechanneling_length(dtheta_x, popt_ch, mu0=None, sigma0=None, L=None, bins=200,
                        hist_range=None, angle_unit='urad', tag='', min_counts=1,
                        peak0_window=(-150e-6, 150e-6), verbose=True):
    """
    Dechannelling length da fit lineare di ln N vs dtheta_x.

    Il fit è fatto con retta (pesi Poisson sui bin) nell'intervallo
    [mu0 + 3 sigma0 , mu_ch - 3 sigma_ch], dove
      - (mu0, sigma0)   = picco a zero (particelle non canalizzate/amorfe),
      - (mu_ch, sigma_ch) = picco di channelling (popt_ch = [A, mu, sigma]).

    Modello: N(Δθ) = N0 exp(-(Δθ)/(θ_b * L_D))  ->  ln N = a + b Δθ,
    con b = -1/(θ_b L_D)  =>  L_D = -1/(b θ_b)   (θ_b = bending_angle / L, rad/m).
    Qui Δθ = θ_out - θ_in è proporzionale alla lunghezza percorsa in channelling
    prima del dechannelling: Δθ = θ_b * z.

    Parametri
    ---------
    dtheta_x   : array dtheta_x (rad) della selezione finale (final_data / final_data_s)
    popt_ch    : popt del fit gaussiano del picco di channelling [A, mu, sigma] (rad)
    mu0, sigma0: media e sigma del picco a zero (rad). Se None, vengono stimati con
                 un fit gaussiano nella finestra peak0_window.
    L          : lunghezza del cristallo [m] (default: crystal_length globale)
    """
    scale = 1e6 if angle_unit == 'urad' else 1.0
    unit_str = r'\mu rad' if angle_unit == 'urad' else 'rad'
    if L is None:
        L = crystal_length
                     # rad/m (curvatura)

    mu_ch, sigma_ch = popt_ch[1], abs(popt_ch[2])

    # --- picco a zero: stima (mu0, sigma0) se non forniti ---
    if mu0 is None or sigma0 is None:
        sel = dtheta_x[(dtheta_x > peak0_window[0]) & (dtheta_x < peak0_window[1])]
        c0, e0 = np.histogram(sel, bins=200)
        x0 = 0.5 * (e0[:-1] + e0[1:])
        p0, _ = curve_fit(gaussian, x0, c0,
                          p0=[c0.max(), np.mean(sel), np.std(sel)],
                          sigma=np.sqrt(np.maximum(c0, 1.0)), absolute_sigma=True,
                          maxfev=10000)
        mu0, sigma0 = p0[1], abs(p0[2])

    # lo = mu0 + 6 * sigma0
    lo = 1000e-6  # rad
    hi = mu_ch - 3 * sigma_ch
    if hi <= lo:
        raise ValueError(f"Intervallo di fit vuoto: lo={lo*scale:.1f}, hi={hi*scale:.1f} {angle_unit}")

    # --- istogramma ---
    if hist_range is None:
        hist_range = (dtheta_x.min(), dtheta_x.max())
    else:
        hist_range = (hist_range[0] / scale, hist_range[1] / scale)
    counts, edges = np.histogram(dtheta_x, bins=bins, range=hist_range)
    centers = 0.5 * (edges[:-1] + edges[1:])

    fit_mask = (centers >= lo) & (centers <= hi) & (counts >= min_counts)
    if fit_mask.sum() < 3:
        raise ValueError(f"Troppi pochi bin nel fit ({fit_mask.sum()}): cambia bins/min_counts")

    xf = centers[fit_mask]
    yf = np.log(counts[fit_mask])
    sig_lnN = 1.0 / np.sqrt(counts[fit_mask])        # Poisson: σ(ln N) = 1/√N

    # --- fit lineare pesato: ln N = a + b x ---
    line = lambda x, a, b: a + b * x
    popt_l, pcov_l = curve_fit(line, xf, yf, sigma=sig_lnN, absolute_sigma=True)
    a, b = popt_l
    a_err, b_err = np.sqrt(np.diag(pcov_l))

    chi2 = np.sum(((yf - line(xf, *popt_l)) / sig_lnN) ** 2)
    ndf = len(xf) - 2
    chi2_ndf = chi2 / ndf if ndf > 0 else np.nan

    # --- dechannelling length: L_D = -1/(b θ_b) ---
    L_D = - L / (b * bending_angle)
    L_D_err = L_D * (b_err / abs(b))

    if verbose:
        print('-' * 50)
        print(f"Dechannelling fit [{tag}]: range = [{lo*scale:.1f}, {hi*scale:.1f}] {angle_unit} "
              f"({fit_mask.sum()} bins)")
        print(f"  peak0: mu0 = {mu0*scale:.2f}, sigma0 = {sigma0*scale:.2f} {angle_unit} | "
              f"peak_ch: mu = {mu_ch*scale:.2f}, sigma = {sigma_ch*scale:.2f} {angle_unit}")
        print(f"  slope b = {b/scale:.4e} ± {b_err/scale:.1e} per {angle_unit} | chi2/ndf = {chi2_ndf:.2f}")
        print(f"  L_D = {L_D*1e2:.2f} ± {L_D_err*1e2:.2f} cm")

    # --- grafico ---
    fig, ax = plt.subplots(figsize=(8, 4.5))
    good = counts > 0
    ax.errorbar(centers[good] * scale, np.log(counts[good]),
                yerr=1.0 / np.sqrt(counts[good]), fmt='o', ms=3, color='gray',
                alpha=0.5, label='data (all bins)')
    ax.errorbar(xf * scale, yf, yerr=sig_lnN, fmt='o', ms=3, color=dark_color,
                capsize=2, label='bins used in fit')
    xl = np.linspace(lo, hi, 200)
    ax.plot(xl * scale, line(xl, *popt_l), color='red', lw=2, ls='--', zorder=10,
            label=rf'linear fit: $L_D$ = {L_D*1e2:.1f} ± {L_D_err*1e2:.1f} cm')
    for v in (lo, hi):
        ax.axvline(v * scale, color='k', ls=':', lw=1)
    ax.axvline(mu0 * scale, color='blue', ls='--', lw=0.8, label=r'$\mu_0$')
    ax.axvline(mu_ch * scale, color='green', ls='--', lw=0.8, label=r'$\mu_{ch}$')
    ax.set_xlabel(rf'$\Delta\theta_x$ [${unit_str}$]')
    ax.set_ylabel(r'$\ln N$')
    ax.set_title(rf'Dechannelling length fit - {ROUTINE.upper()}')
    ax.grid(alpha=0.3)
    ax.legend(fontsize=11)
    fig.tight_layout()
    plt.savefig(f'{output_folder}/dechanneling_fit_{tag}.png', dpi=300)
    plt.savefig(f'{output_folder}/dechanneling_fit_{tag}.pdf')
    plt.close(fig)

    return L_D, L_D_err, popt_l, pcov_l, chi2_ndf


# ANALYSIS ========================================================================================

beam_distribution_plot(x_in, y_in, theta_x_in, theta_y_in, label='all')
# beam_distribution_plot(x_in_survived, y_in_survived, theta_x_in_survived, theta_y_in_survived, label='- survived particles')

ch_footprint_plot(x_in_survived, y_in_survived, x_out_survived, y_out_survived, dtheta_x, crystal_x, crystal_y, 0.2, bending_angle, label='survived')
# ch_footprint_plot(x_in, y_in, x_out, y_out, dtheta_x, crystal_x, crystal_y, 0.2, bending_angle, label='all')

# pos_x, sigmas_x, errs_x, n_x = mcs_edge_scan(
#     pos_in=survived_data['x_in'], dtheta=survived_data['dtheta_y'], crystal_edge=crystal_x,
#     other_pos_in=survived_data['y_in'], other_edge=crystal_y,
#     plane_label='x', dtheta_label='x', n_slices=50, fit_window=30, label='- x plane')

# pos_y, sigmas_y, errs_y, n_y = mcs_edge_scan(
#     pos_in=survived_data['y_in'], dtheta=survived_data['dtheta_y'], crystal_edge=crystal_y,
#     other_pos_in=survived_data['x_in'], other_edge=crystal_x,
#     plane_label='y', dtheta_label='y', n_slices=50, fit_window=None, label='- y plane')


# geometric cut
EDGE_CUT_RIGHT = 0.22e-3 if ROUTINE == 'bdsim' else 0.0   # m
def geometric_mask(x, y):
    return ((x > -crystal_x/2) &
            (x <  crystal_x/2 - EDGE_CUT_RIGHT) &
            (np.abs(y) < crystal_y/2))

geometric_hit     = geometric_mask(survived_data['x_in'], survived_data['y_in'])
geometric_hit_all = geometric_mask(x_in, y_in)
geom_data = apply_selection(survived_data, geometric_hit)
geom_data_s = apply_selection(survived_data_smeared, geometric_hit)

print('-'*50)
print(f"Sopravvissute che hanno colpito geometricamente il cristallo: {geometric_hit.sum()} = {geometric_hit.sum()/len(survived_data['x_in'])*100:.1f}%")
print(f"Sopravvissute che hanno mancato il cristallo: {(~geometric_hit).sum()} = {(~geometric_hit).sum()/len(survived_data['x_in'])*100:.1f}%")
# beam_distribution_plot(geom_data['x_in'], geom_data['y_in'], geom_data['theta_x_in'], geom_data['theta_y_in'], label='- entered particles')

vr_mask = (dtheta_x < -0.3*theta_L) & (dtheta_x > -3*theta_L)
n_vr = vr_mask.sum()
# print(f"Volume reflection: {n_vr}/{len(survived_data['x_in'])} = {n_vr/len(survived_data['x_in'])*100:.1f}%")

# lindhard selection
lindhard_cut = np.abs(geom_data['theta_x_in']) < (0.5 * theta_L); lindhard_all = np.abs(theta_x_in) < (0.5 * theta_L)
lindhard_cut_s = np.abs(geom_data_s['theta_x_in']) < (0.5 * theta_L); lindhard_all_s = np.abs(theta_x_in_meas_all) < (0.5 * theta_L)

final_data = apply_selection(geom_data, lindhard_cut)
final_data_s = apply_selection(geom_data_s, lindhard_cut_s)

# also keep the particles with the potential to be channeled (N_tot)
attempted_mask = geometric_hit_all & lindhard_all
attempted_ids = particle_id_in[attempted_mask]
n_tot_correct = attempted_mask.sum()

lost_mask = ~np.isin(attempted_ids, idx_out)
n_lost_attempted = lost_mask.sum()
# print(f"Persi tra i tagli geometrici + Lindhard: {n_lost_attempted}/{n_tot_correct} = {n_lost_attempted/n_tot_correct*100:.1f}%")

attempted_mask_s = geometric_hit_all & lindhard_all_s
n_tot_correct_s  = attempted_mask_s.sum()

n_geom = geometric_hit_all.sum()

# f: frazione di fascio nella finestra di Lindhard (angoli veri, non smeared)
f = attempted_mask.sum() / n_geom
print(f"f = {attempted_mask.sum()}/{n_geom} = {f*100:.1f}%")

# P_inel: frazione persa, sulle geometriche
out_win = geometric_hit_all & ~lindhard_all
P_inel_all = lost_in[geometric_hit_all].sum() / n_geom
P_inel_out = lost_in[out_win].sum() / out_win.sum()
P_inel_in  = lost_in[attempted_mask].sum() / attempted_mask.sum()

print(f"P_inel (tutte le geometriche)   = {P_inel_all*100:.1f}%")
print(f"P_inel (fuori finestra Lindhard) = {P_inel_out*100:.1f}%")
print(f"P_inel (dentro finestra)         = {P_inel_in*100:.1f}%")

# EFFICIENCY ========================================================================================
n_channeled,   n_ch_err,   popt,   pcov,   fit_bin_width   = channeling_efficiency(final_data['theta_x_in'],   final_data['dtheta_x'])
n_channeled_s, n_ch_err_s, popt_s, pcov_s, fit_bin_width_s = channeling_efficiency(final_data_s['theta_x_in'], final_data_s['dtheta_x'])

eff_ns,   eff_s     = n_channeled / n_tot_correct,     n_channeled_s / n_tot_correct_s
stat_ns,  stat_s    = n_ch_err / n_tot_correct,        n_ch_err_s / n_tot_correct_s       # errore statistico (fit)
binom_ns, binom_s   = np.sqrt(eff_ns*(1-eff_ns)/n_tot_correct), np.sqrt(eff_s*(1-eff_s)/n_tot_correct_s)  # solo conteggio
print(f"[not smeared] eff = ({eff_ns*100:.2f} ± {stat_ns*100:.2f} [fit]) %   (binomiale: ±{binom_ns*100:.2f})")
print(f"[smeared]     eff = ({eff_s*100:.2f} ± {stat_s*100:.2f} [fit]) %   (binomiale: ±{binom_s*100:.2f})")

eff3_ns, err3_ns = eff_3sigma(final_data,   popt,   n_tot_correct)
eff3_s,  err3_s  = eff_3sigma(final_data_s, popt_s, n_tot_correct_s)
print(f"[not smeared] eff(3σ) = ({eff3_ns*100:.2f} ± {err3_ns*100:.2f}) %")
print(f"[smeared]     eff(3σ) = ({eff3_s*100:.2f} ± {err3_s*100:.2f}) %")


# controllo indipendente (opzionale, ~secondi)
boot_ns, _ = efficiency_bootstrap(final_data['dtheta_x'],   n_tot_correct)
boot_s,  _ = efficiency_bootstrap(final_data_s['dtheta_x'], n_tot_correct_s)

L_D_s, L_D_err_s, *_ = dechanneling_length(final_data_s['dtheta_x'], popt_s, hist_range=(-750,6200), tag='smeared')
# L_D_ns, L_D_err_ns, *_ = dechanneling_length(final_data['dtheta_x'], popt, hist_range=(-750,6200), tag='original')
 


# # --- GUARD RAIL: verifica affidabilità del fit ---
# sigma_err = np.sqrt(pcov_s[2, 2]) if pcov_s is not None else np.inf
# sigma = abs(popt_s[2]) if popt_s is not None else 0
# fit_reliable = (
#     pcov_s is not None
#     and np.all(np.isfinite(pcov_s))
#     and sigma_err < sigma          # errore sul sigma minore del sigma stesso
#     and n_channeled_s <= n_tot_correct_s  # non puoi channelare più particelle di quelle disponibili
# )
# if not fit_reliable:
#     print("ATTENZIONE: fit inaffidabile (statistica insufficiente) — efficienza non significativa")
#     efficiency_correct = np.nan
# else:
#     efficiency_correct = n_channeled_s / n_tot_correct_s
#     eff_err_binom = np.sqrt(efficiency_correct*(1-efficiency_correct)/n_tot_correct_s)
#     print(f"[smeared] Channeling efficiency: {n_channeled_s:.0f}/{n_tot_correct_s} = ({efficiency_correct*100:.1f}± {eff_err_binom*100:.1f})%")
# # --- fine guard rail ---


# angular_scan_plot(final_data['theta_x_in'], final_data['theta_x_out'], theta_L, label='', popt=popt, fit_bin_width=fit_bin_width)
# angular_scan_plot(survived_data['theta_x_in'], survived_data['theta_x_out'], theta_L, label='(survived particles)', popt=popt, fit_bin_width=fit_bin_width)
angular_scan_plot(geom_data_s['theta_x_in'], geom_data_s['theta_x_out'], theta_L, label='entered', popt=popt, fit_bin_width=fit_bin_width_s)
angular_scan_plot(final_data_s['theta_x_in'], final_data_s['theta_x_out'], theta_L, label='smeared', popt=popt_s, fit_bin_width=fit_bin_width_s)

# smear_comparison_plot(geom_data['dtheta_x'], geom_data_s['dtheta_x'], popt=popt, popt_smeared=popt_s, fit_bin_width=fit_bin_width, fit_bin_width_smeared=fit_bin_width_s, label='geometric')
smear_comparison_plot(final_data['dtheta_x'], final_data_s['dtheta_x'], label='final_x')
smear_comparison_plot2(survived_data['theta_x_in'], survived_data_smeared['theta_x_in'], survived_data['theta_y_in'], survived_data_smeared['theta_y_in'])

for ax, th, th_m in [('x', theta_x_in, theta_x_in_meas_all), ('y', theta_y_in, theta_y_in_meas_all)]:
    print(f"theta_{ax}: sigma originale = {th.std()*1e6:.2f} urad, smeared = {th_m.std()*1e6:.2f} urad, "
          f"atteso = {np.hypot(th.std(), s_in)*1e6:.2f} urad")

# --- non smeared ---
mu_ns, sig_ns = popt[1], abs(popt[2])
efficiency_maps(x_in[attempted_mask], y_in[attempted_mask],
                final_data['x_in'], final_data['y_in'], final_data['dtheta_x'],
                mu_ns, sig_ns, crystal_x, crystal_y, nx_bins=5, ny_bins=7, tag='original')
profiles_plot(x_in[attempted_mask], y_in[attempted_mask], theta_x_in[attempted_mask],
              final_data, mu_ns, sig_ns, crystal_x, crystal_y, theta_L, tag='original')

# --- smeared (denominatore con lindhard_all_s, già in attempted_mask_s) ---
mu_s, sig_s = popt_s[1], abs(popt_s[2])
efficiency_maps(x_in[attempted_mask_s], y_in[attempted_mask_s],
                final_data_s['x_in'], final_data_s['y_in'], final_data_s['dtheta_x'],
                mu_s, sig_s, crystal_x, crystal_y, nx_bins=5, ny_bins=7, tag='smeared')
profiles_plot(x_in[attempted_mask_s], y_in[attempted_mask_s], theta_x_in_meas_all[attempted_mask_s],
              final_data_s, mu_s, sig_s, crystal_x, crystal_y, theta_L, tag='smeared')


np.save(f'{output_folder}/dtheta_x_final_smeared_{ROUTINE}.npy', final_data_s['dtheta_x'])
print(f"Salvato array dtheta_x per la routine {ROUTINE}")

print('='*50)
print(f"Plot saved in {output_folder}")
