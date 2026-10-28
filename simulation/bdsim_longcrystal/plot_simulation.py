import sys

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from scipy.optimize import curve_fit
import os

viridis_yellow = "#fde725"
viridis_green = "#5ec962"
viridis_teal = "#21918c"
viridis_blue = "#3b528b"
viridis_purple = "#440154"

n = sys.argv[1] if len(sys.argv) > 1 else '6'
filename = f'./runs/cry2_merged{n}_slim.npz'
data = np.load(filename)

output_folder = f'./plots{n}'
os.makedirs(output_folder, exist_ok=True)

    # cov_avg=cov_avg,
    # mu_avg=mu_avg,

# reconstruct the arrays from the saved data
x_in = data['x_in']; theta_x_in = data['px_in']
y_in = data['y_in']; theta_y_in = data['py_in']
# zeta_in = data['zeta_in']; delta_in = data['delta_in']; ptau_in = data['ptau_in']

x_out=data['x_out']; theta_x_out=data['px_out']
y_out=data['y_out']; theta_y_out=data['py_out']
# zeta_out=data['zeta_out']; delta_out=data['delta_out']; ptau_out=data['ptau_out']

particle_id_in = data['particle_id_in']
# state_in = data['state_in']; pdg_id_in = data['pdg_id_in']
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
survived = state_out == 1
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
print(f"Survived particles: {np.sum(survived)} = {np.sum(survived)/len(x_in)*100:.2f}%")

theta_L = 12.992e-6  # rad Lindhard critical angle for Si(110) at 180 GeV

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

s_in  = np.hypot(sigma_track, highland(p_GeV, L_up   / X0_air))
s_out = np.hypot(sigma_track, highland(p_GeV, L_down / X0_air))

# 1) θ_in misurato per TUTTE le particelle valide (un'estrazione sola)
theta_x_in_meas_all = theta_x_in + rng_smear.standard_normal(theta_x_in.size, dtype=np.float32) * s_in
theta_y_in_meas_all = theta_y_in + rng_smear.standard_normal(theta_y_in.size, dtype=np.float32) * s_in

# 2) stesso valore propagato alle sopravvissute (stessa mappatura già usata sopra)
theta_x_in_meas_surv = theta_x_in_meas_all[src]
theta_y_in_meas_surv = theta_y_in_meas_all[src]

# 3) θ_out misurato (estrazioni indipendenti, solo sulle sopravvissute)
theta_x_out_meas_surv = theta_x_out_survived + rng_smear.standard_normal(theta_x_out_survived.size, dtype=np.float32) * s_out
theta_y_out_meas_surv = theta_y_out_survived + rng_smear.standard_normal(theta_y_out_survived.size, dtype=np.float32) * s_out

survived_data_smeared = dict(survived_data)
survived_data_smeared.update({
    'theta_x_in':  theta_x_in_meas_surv,
    'theta_y_in':  theta_y_in_meas_surv,
    'theta_x_out': theta_x_out_meas_surv,
    'theta_y_out': theta_y_out_meas_surv,
    'dtheta_x':    theta_x_out_meas_surv - theta_x_in_meas_surv,
    'dtheta_y':    theta_y_out_meas_surv - theta_y_in_meas_surv,
})
print(f"Smear: s_in = {s_in*1e6:.2f} urad, s_out = {s_out*1e6:.2f} urad")

def beam_distribution_plot(x_in, y_in, theta_x_in, theta_y_in, label=None, distance_unit='mm', angle_unit='urad'):
    scale = 1e6 if angle_unit == 'urad' else 1.0
    unit_str = r'\mu rad' if angle_unit == 'urad' else 'rad'

    distance_scale = 1e3 if distance_unit == 'mm' else 1.0

    # x_lim = (-3, 3) if distance_unit == 'mm' else (x_in.min() * distance_scale, x_in.max() * distance_scale)
    # y_lim = (-11, 11) if distance_unit == 'mm' else (y_in.min() * distance_scale, y_in.max() * distance_scale)
    x_lim = (x_in.min() * distance_scale, x_in.max() * distance_scale)
    y_lim = (y_in.min() * distance_scale, y_in.max() * distance_scale)
    theta_x_lim = (-150, 150) if angle_unit == 'urad' else (theta_x_in.min() * scale, theta_x_in.max() * scale)
    theta_y_lim = (-150, 150) if angle_unit == 'urad' else (theta_y_in.min() * scale, theta_y_in.max() * scale)
    # theta_x_lim = (theta_x_in.min() * scale, theta_x_in.max() * scale)
    # theta_y_lim = (theta_y_in.min() * scale, theta_y_in.max() * scale)

    fig, graph = plt.subplots(2, 2, figsize=(10, 10))
    # 1. Grafico (x, y) - Beam impact position
    h0 = graph[0, 0].hist2d(x_in * distance_scale, y_in * distance_scale, bins=100, cmap='viridis', norm=LogNorm(), rasterized=True)
    fig.colorbar(h0[3], ax=graph[0, 0], label='Counts')
    graph[0, 0].set_xlabel(f'impact x [{distance_unit}]')
    graph[0, 0].set_ylabel(f'impact y [{distance_unit}]')
    graph[0, 0].set_title(f'Beam Impact Position (x, y) {label}')
    graph[0, 0].set_xlim(x_lim)
    graph[0, 0].set_ylim(y_lim)
    graph[0, 0].grid(alpha=0.3)

    h00 = graph[0, 1].hist2d(theta_x_in * scale, theta_y_in * scale, bins=100, cmap='viridis', norm=LogNorm(), rasterized=True)
    fig.colorbar(h00[3], ax=graph[0, 1], label='Counts')
    graph[0, 1].set_xlabel(rf'$\theta_x$ [${unit_str}$]')
    graph[0, 1].set_ylabel(rf'$\theta_y$ [${unit_str}$]')
    graph[0, 1].set_title(f'Angular Distribution {label}')
    graph[0, 1].set_xlim(theta_x_lim)
    graph[0, 1].set_ylim(theta_y_lim)
    graph[0, 1].grid(alpha=0.3)

    h1 = graph[1, 0].hist2d(x_in * distance_scale, theta_x_in * scale, bins=100, cmap='viridis', norm=LogNorm(), rasterized=True)
    fig.colorbar(h1[3], ax=graph[1, 0], label='Counts')
    graph[1, 0].set_xlabel(f'impact x [{distance_unit}]')
    graph[1, 0].set_ylabel(rf'$\theta_x$ [${unit_str}$]')
    graph[1, 0].set_title(f'Horizontal Phase Space (x, $θ_x$) {label}')
    graph[1, 0].set_xlim(x_lim)
    graph[1, 0].set_ylim(theta_x_lim)
    graph[1, 0].grid(alpha=0.3)

    h2 = graph[1, 1].hist2d(y_in * distance_scale, theta_y_in * scale, bins=100, cmap='viridis', norm=LogNorm(), rasterized=True)
    fig.colorbar(h2[3], ax=graph[1, 1], label='Counts')
    graph[1, 1].set_xlabel(f'impact y [{distance_unit}]')
    graph[1, 1].set_ylabel(rf'$\theta_y$ [${unit_str}$]')
    graph[1, 1].set_title(f'Vertical Phase Space (y, $θ_y$) {label}')
    graph[1, 1].set_xlim(y_lim)
    graph[1, 1].set_ylim(theta_y_lim)
    graph[1, 1].grid(alpha=0.3)

    fig.tight_layout()
    label = label[2:]
    plt.savefig(f'{output_folder}/impact_distribution_{label}.png', dpi=300)
    plt.savefig(f'{output_folder}/impact_distribution_{label}.pdf')
    plt.close(fig)

def angular_scan_plot(theta_x_in, theta_x_out, theta_L, popt=None, fit_bin_width=None, label=None, angle_unit='urad'):

    scale = 1e6 if angle_unit == 'urad' else 1.0
    unit_str = r'\mu rad' if angle_unit == 'urad' else 'rad'

    dtheta_x = theta_x_out - theta_x_in

    fig, graph = plt.subplots(1, 3, figsize=(18, 5))

    h = graph[0].hist2d(theta_x_in * scale, dtheta_x * scale, bins=300, cmap='viridis', norm=LogNorm(), rasterized=True)
    cb = fig.colorbar(h[3], ax=graph[0])
    cb.set_label('Counts')
    graph[0].set_ylim(-750, 6200)
    if label=='':
        graph[0].set_xlim(-7, 7)
    # graph[0].axhline(0, color='gray', lw=0.8, ls='--')
    graph[0].axvline(-0.5 * theta_L * scale, color='red', lw=0.8, ls='--', label=r'$\pm 1/2 \theta_L$')
    graph[0].axvline(0.5 * theta_L * scale, color='red', lw=0.8, ls='--')
    graph[0].set_xlabel(rf'$\theta_{{in}}$ [${unit_str}$]')
    graph[0].set_ylabel(rf'$\Delta\theta_x = \theta_{{out}} - \theta_{{in}}$ [${unit_str}$]')
    graph[0].set_title(f'Angular scan')
    graph[0].grid(alpha=0.3)
    graph[0].legend()

    counts, bins, _ = graph[1].hist(dtheta_x * scale, bins=200, color=viridis_green, rasterized=True)
    counts2, bins2, _ = graph[2].hist(dtheta_x * scale, bins=7250, color=viridis_green, rasterized=True)
    plot_bin_width = bins[1] - bins[0]
    plot_bin_width2 = bins2[1] - bins2[0]

    graph[1].set_xlim(-750, 6200)
    graph[1].set_yscale('log')
    graph[1].set_ylim(500, counts.max()*1.2)
    graph[1].set_xlabel(rf'$\Delta\theta_x$ [${unit_str}$]')
    graph[1].set_ylabel('Counts')
    graph[1].set_title(rf'$\Delta\theta_x$ distribution')
    graph[1].grid(alpha=0.3)

    graph[2].set_xlim(5950, 6050)
    graph[2].set_ylim(0, counts2.max()*1.2)
    graph[2].set_xlabel(rf'$\Delta\theta_x$ [${unit_str}$]')
    graph[2].set_ylabel('Counts')
    graph[2].set_title(rf'$\Delta\theta_x$ distribution (zoom)')
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

    counts_smeared, bins_smeared, _ = graph[0].hist(dtheta_x_smeared * scale, bins=200, range=(-750, 6200), color=viridis_blue, label='smeared', rasterized=True)
    counts, bins, _ = graph[0].hist(dtheta_x * scale, bins=200, range=(-750, 6200), color=viridis_green, alpha=0.5, label='original', rasterized=True)

    if label[-1] == 'y':
        graph[0].set_xlim(-250, 250)
    # graph[0].set_xlim(-750, 6200)
    graph[0].set_yscale('log')
    graph[0].set_ylim(500, max(counts.max(), counts_smeared.max())*1.2)
    graph[0].set_xlabel(rf'$\Delta\theta_{label_ax}$ [${unit_str}$]')
    graph[0].set_ylabel('Counts')
    graph[0].set_title(rf'$\Delta\theta_{label_ax}$ distribution {label}')
    graph[0].grid(alpha=0.3)
    graph[0].legend()

    graph[1].hist(dtheta_x * scale, bins=200, range=(5950, 6050), color=viridis_green, label='original', rasterized=True)
    graph[1].hist(dtheta_x_smeared * scale, bins=200, range=(5950, 6050), color=viridis_blue, alpha=0.5, label='smeared', rasterized=True)

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

    counts_smeared, _, _ = graph[0].hist(x_smeared * scale, bins=200, range=(-150, 150), color=viridis_blue, label='smeared')
    counts, _, _ = graph[0].hist(x * scale, bins=200, range=(-150, 150), color=viridis_green, alpha=0.7, label='original')

    # graph[0].set_ylim(500, max(counts.max(), counts_smeared.max())*1.2)
    graph[0].set_xlabel(rf'$\theta_{{in,x}}$ [${unit_str}$]')
    graph[0].set_ylabel('Counts')
    graph[0].set_title(rf'Incoming $\theta_{{x}}$ distribution')
    graph[0].grid(alpha=0.3)
    graph[0].legend()

    graph[1].hist(y_smeared * scale, bins=200, range=(-150, 150), color=viridis_blue, label='smeared')
    graph[1].hist(y * scale, bins=200, range=(-150, 150), color=viridis_green, alpha=0.7, label='original')

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

def channeling_efficiency(theta_x_in, dtheta_x, bending_angle=bending_angle, tol_frac=0.4, angle_unit='urad'):
    scale = 1e6 if angle_unit == 'urad' else 1.0
    unit_str = r'\mu rad' if angle_unit == 'urad' else 'rad'

    n_tot = len(theta_x_in)

    counts, edges = np.histogram(dtheta_x, bins=7250)
    centers = 0.5 * (edges[:-1] + edges[1:])
    bin_width = edges[1] - edges[0]

    search_window = 0.5 * abs(bending_angle)
    search_mask = np.abs(centers - bending_angle) < search_window

    if search_mask.sum() > 0:
        peak_x = centers[search_mask][np.argmax(counts[search_mask])]
    else:
        peak_x = bending_angle

    window = 0.5 * abs(bending_angle)
    # mask = (centers >= (peak_x - 150/scale)) & (centers < (peak_x + window))
    mask = (centers >= (5990/scale)) & (centers < (peak_x + 50/scale))
    while mask.sum() < 5 and window < (dtheta_x.max() - dtheta_x.min()):
        window *= 2
        mask = np.abs(centers - bending_angle) < window
 
    x_fit, y_fit = centers[mask], counts[mask]
    A0 = y_fit.max() if len(y_fit) else 1.0
    sigma0 = 12e-6  # initial guess for sigma in radians
 
    popt, pcov = curve_fit(gaussian, x_fit, y_fit, p0=[A0, bending_angle, sigma0], maxfev=10000)
 
    A, mu, sigma = popt
    sigma = abs(sigma)
    perr = np.sqrt(np.diag(pcov))
    mu_err = perr[1]
    sigma_err = perr[2]
 
    n_channeled = A * sigma * np.sqrt(2 * np.pi) / bin_width
    efficiency = n_channeled / n_tot if n_tot > 0 else np.nan

    n_3sigma = np.sum(np.abs(dtheta_x - mu) < 3 * sigma)
    n_2sigma = np.sum(np.abs(dtheta_x - mu) < 2 * sigma)
    n_above  = np.sum(dtheta_x > (mu - 2 * sigma))

    print('-'*50)
    # print(f"Channeling efficiency: {n_channeled:.0f}/{n_tot} = {efficiency*100:.1f}%")
    # print(f"Channeling efficiency (within 3 sigma): {n_3sigma:.0f}/{n_tot} = {n_3sigma/n_tot*100:.1f}%")
    # print(f"Channeling efficiency (within 2 sigma): {n_2sigma:.0f}/{n_tot} = {n_2sigma/n_tot*100:.1f}%")
    # print(f"Channeling efficiency (above 2sigma): {n_above:.0f}/{n_tot} = {n_above/n_tot*100:.1f}%")
    print(f"Channeling peak = {mu * scale:.2f} ± {mu_err * scale:.2f} urad | Sigma = {sigma * scale:.2f} ± {sigma_err * scale:.2f}")

    return n_channeled, popt, pcov, bin_width

# ANALYSIS

beam_distribution_plot(x_in, y_in, theta_x_in, theta_y_in, label='- all particles')
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
geometric_hit = (np.abs(survived_data['x_in']) < crystal_x/2) & (np.abs(survived_data['y_in']) < crystal_y/2)
geometric_hit_all = (np.abs(x_in) < crystal_x/2) & (np.abs(y_in) < crystal_y/2)
geom_data = apply_selection(survived_data, geometric_hit)
geom_data_s = apply_selection(survived_data_smeared, geometric_hit)

print('-'*50)
print(f"Sopravvissute che hanno colpito geometricamente il cristallo: {geometric_hit.sum()} = {geometric_hit.sum()/len(survived_data['x_in'])*100:.1f}%")
print(f"Sopravvissute che hanno mancato il cristallo: {(~geometric_hit).sum()} = {(~geometric_hit).sum()/len(survived_data['x_in'])*100:.1f}%")
# beam_distribution_plot(geom_data['x_in'], geom_data['y_in'], geom_data['theta_x_in'], geom_data['theta_y_in'], label='- entered particles')

vr_mask = (dtheta_x < -0.3*theta_L) & (dtheta_x > -3*theta_L)
n_vr = vr_mask.sum()
print(f"Volume reflection: {n_vr}/{len(survived_data['x_in'])} = {n_vr/len(survived_data['x_in'])*100:.1f}%")

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
print(f"Persi tra i tagli geometrici + Lindhard: {n_lost_attempted}/{n_tot_correct} = {n_lost_attempted/n_tot_correct*100:.1f}%")

attempted_mask_s = geometric_hit_all & lindhard_all_s
n_tot_correct_s  = attempted_mask_s.sum()

# efficiency
n_channeled, popt, pcov, fit_bin_width = channeling_efficiency(final_data['theta_x_in'], final_data['dtheta_x'])
n_channeled_s, popt_s, pcov_s, fit_bin_width_s = channeling_efficiency(final_data_s['theta_x_in'], final_data_s['dtheta_x'])

efficiency_correct_s = n_channeled_s / n_tot_correct_s
print(f"[smeared] Channeling efficiency: {n_channeled_s:.0f}/{n_tot_correct_s} = {efficiency_correct_s*100:.1f}%")


# --- GUARD RAIL: verifica affidabilità del fit ---
sigma_err = np.sqrt(pcov_s[2, 2]) if pcov_s is not None else np.inf
sigma = abs(popt_s[2]) if popt_s is not None else 0
fit_reliable = (
    pcov_s is not None
    and np.all(np.isfinite(pcov_s))
    and sigma_err < sigma          # errore sul sigma minore del sigma stesso
    and n_channeled_s <= n_tot_correct_s  # non puoi channelare più particelle di quelle disponibili
)
if not fit_reliable:
    print("ATTENZIONE: fit inaffidabile (statistica insufficiente) — efficienza non significativa")
    efficiency_correct = np.nan
else:
    efficiency_correct = n_channeled_s / n_tot_correct_s
    eff_err_binom = np.sqrt(efficiency_correct*(1-efficiency_correct)/n_tot_correct_s)
    print(f"[smeared] Channeling efficiency: {n_channeled_s:.0f}/{n_tot_correct_s} = ({efficiency_correct*100:.1f}± {eff_err_binom*100:.1f})%")
# --- fine guard rail ---

efficiency_notsmeared = n_channeled / n_tot_correct
eff_err_binom = np.sqrt(efficiency_notsmeared*(1-efficiency_notsmeared)/n_tot_correct)
print(f"[not smeared] Channeling efficiency: {n_channeled:.0f}/{n_tot_correct} = ({efficiency_notsmeared*100:.1f}± {eff_err_binom*100:.1f})%")


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

print('='*50)
print(f"Plot saved in {output_folder}")
