import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from scipy.optimize import curve_fit

# filename = 'cry1_20260916_11.npz'
filename = 'cry1_20260917_16.npz'
data = np.load(filename)
output_folder = 'plots_cry1'

    # cov_avg=cov_avg,
    # mu_avg=mu_avg,

# reconstruct the arrays from the saved data
x_in = data['x_in']; theta_x_in = data['px_in']
y_in = data['y_in']; theta_y_in = data['py_in']
zeta_in = data['zeta_in']; delta_in = data['delta_in']; ptau_in = data['ptau_in']

x_out=data['x_out']; theta_x_out=data['px_out']
y_out=data['y_out']; theta_y_out=data['py_out']
zeta_out=data['zeta_out']; delta_out=data['delta_out']; ptau_out=data['ptau_out']

particle_id_in = data['particle_id_in']; state_in = data['state_in']; pdg_id_in = data['pdg_id_in']
particle_id_out = data['particle_id_out']; state_out = data['state_out']; pdg_id_out = data['pdg_id_out']; at_element_out = data['at_element_out']

# single values
npart = int(data['npart'])
crystal_length = float(data['crystal_length'])
crystal_x = float(data['crystal_x'])
crystal_y = float(data['crystal_y'])
bending_angle = float(data['bending_angle'])
crystal_material = str(data['crystal_material'])
p0c = float(data['p0c'])
mass0 = float(data['mass0'])
label = str(data['label'])


# LOST PARTICLES
# where lost particles got lost?
print('='*50)
unique_states, counts = np.unique(state_out, return_counts=True)
for s, c in zip(unique_states, counts):
    print(f"state={s}: {c} particles")
lost_333 = state_out == -333
elements, counts = np.unique(at_element_out[lost_333], return_counts=True)
for e, c in zip(elements, counts):
    print(f"at_element={e}: {c} lost particles at -333")

# where they would have been?
particle_id_in_valid = data['particle_id_in'][data['state_in'] == 1]
x_in_valid = data['x_in'][data['state_in'] == 1]
y_in_valid = data['y_in'][data['state_in'] == 1]

lost_mask = data['state_out'] == -333
idx_lost = data['particle_id_out'][lost_mask]

order = np.argsort(particle_id_in_valid)
sorted_ids = particle_id_in_valid[order]
pos = np.searchsorted(sorted_ids, idx_lost)
x_in_lost = x_in_valid[order][pos]
y_in_lost = y_in_valid[order][pos]

geometric_hit_lost = (np.abs(x_in_lost) < crystal_x/2) & (np.abs(y_in_lost) < crystal_y/2)
print(f"Frazione delle perse che avrebbe comunque colpito il cristallo: {geometric_hit_lost.mean()*100:.1f}%")
print(f"Numero di particelle perse che avrebbero dovuto MANCARE il cristallo: {(~geometric_hit_lost).sum()}")

# SURVIVED PARTICLES
valid_in = state_in == 1
survived = state_out == 1

x_in, y_in, theta_x_in, theta_y_in = x_in[valid_in], y_in[valid_in], theta_x_in[valid_in], theta_y_in[valid_in]
zeta_in, delta_in, ptau_in = zeta_in[valid_in], delta_in[valid_in], ptau_in[valid_in]
particle_id_in, pdg_id_in = particle_id_in[valid_in], pdg_id_in[valid_in]

# filtering survived particles  
idx_out = particle_id_out[survived]
theta_x_out_survived = theta_x_out[survived]
pdg_id_out_survived = pdg_id_out[survived]

order_in = np.argsort(particle_id_in)
sorted_idx_in = particle_id_in[order_in]
idx = np.searchsorted(sorted_idx_in, idx_out)

x_in_survived = x_in[order_in][idx]
y_in_survived = y_in[order_in][idx]
theta_x_in_survived = theta_x_in[order_in][idx]
theta_y_in_survived = theta_y_in[order_in][idx]
zeta_in_survived = zeta_in[order_in][idx]; delta_in_survived = delta_in[order_in][idx]; ptau_in_survived = ptau_in[order_in][idx]

x_out_survived = x_out[survived]
y_out_survived = y_out[survived]
theta_y_out_survived = theta_y_out[survived]
zeta_out_survived = zeta_out[survived]; delta_out_survived = delta_out[survived]; ptau_out_survived = ptau_out[survived]

assert np.all(sorted_idx_in[idx] == idx_out), "Mismatch nel matching particle_id in/out!"

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
}

print(f"Total particles: {len(x_in)}")
print(f"Survived particles: {np.sum(survived)}")

theta_L1 = 13.228e-6  # rad Lindhard critical angle for Si(110) at 180 GeV
theta_L2 = 12.992e-6  # rad Lindhard critical angle for Si(110) at 180 GeV

def beam_distribution_plot(x_in, y_in, theta_x_in, theta_y_in, label=None, distance_unit='mm', angle_unit='urad'):
    scale = 1e6 if angle_unit == 'urad' else 1.0
    unit_str = r'\mu rad' if angle_unit == 'urad' else 'rad'

    distance_scale = 1e3 if distance_unit == 'mm' else 1.0

    x_lim = (x_in.min() * distance_scale, x_in.max() * distance_scale)
    y_lim = (y_in.min() * distance_scale, y_in.max() * distance_scale)
    theta_x_lim = (-150, 150) if angle_unit == 'urad' else (theta_x_in.min() * scale, theta_x_in.max() * scale)
    theta_y_lim = (-150, 150) if angle_unit == 'urad' else (theta_y_in.min() * scale, theta_y_in.max() * scale)
    # theta_x_lim = (theta_x_in.min() * scale, theta_x_in.max() * scale)
    # theta_y_lim = (theta_y_in.min() * scale, theta_y_in.max() * scale)

    fig, graph = plt.subplots(2, 2, figsize=(10, 10))
    # 1. Grafico (x, y) - Beam impact position
    h0 = graph[0, 0].hist2d(x_in * distance_scale, y_in * distance_scale, bins=100, cmap='viridis', norm=LogNorm())
    fig.colorbar(h0[3], ax=graph[0, 0], label='Counts')
    graph[0, 0].set_xlabel(f'impact x [{distance_unit}]')
    graph[0, 0].set_ylabel(f'impact y [{distance_unit}]')
    graph[0, 0].set_title(f'Beam Impact Position (x, y) {label}')
    graph[0, 0].set_xlim(x_lim)
    graph[0, 0].set_ylim(y_lim)
    graph[0, 0].grid(alpha=0.3)

    h00 = graph[0, 1].hist2d(theta_x_in * scale, theta_y_in * scale, bins=100, cmap='viridis', norm=LogNorm())
    fig.colorbar(h00[3], ax=graph[0, 1], label='Counts')
    graph[0, 1].set_xlabel(rf'$\theta_x$ [${unit_str}$]')
    graph[0, 1].set_ylabel(rf'$\theta_y$ [${unit_str}$]')
    graph[0, 1].set_title(f'Angular Distribution {label}')
    graph[0, 1].set_xlim(theta_x_lim)
    graph[0, 1].set_ylim(theta_y_lim)
    graph[0, 1].grid(alpha=0.3)

    h1 = graph[1, 0].hist2d(x_in * distance_scale, theta_x_in * scale, bins=100, cmap='viridis', norm=LogNorm())
    fig.colorbar(h1[3], ax=graph[1, 0], label='Counts')
    graph[1, 0].set_xlabel(f'impact x [{distance_unit}]')
    graph[1, 0].set_ylabel(rf'$\theta_x$ [${unit_str}$]')
    graph[1, 0].set_title(f'Horizontal Phase Space (x, $θ_x$) {label}')
    graph[1, 0].set_xlim(x_lim)
    graph[1, 0].set_ylim(theta_x_lim)
    graph[1, 0].grid(alpha=0.3)

    h2 = graph[1, 1].hist2d(y_in * distance_scale, theta_y_in * scale, bins=100, cmap='viridis', norm=LogNorm())
    fig.colorbar(h2[3], ax=graph[1, 1], label='Counts')
    graph[1, 1].set_xlabel(f'impact y [{distance_unit}]')
    graph[1, 1].set_ylabel(rf'$\theta_y$ [${unit_str}$]')
    graph[1, 1].set_title(f'Vertical Phase Space (y, $θ_y$) {label}')
    graph[1, 1].set_xlim(y_lim)
    graph[1, 1].set_ylim(theta_y_lim)
    graph[1, 1].grid(alpha=0.3)

    fig.tight_layout()
    label = label[2:]
    plt.savefig(f'{output_folder}/gaus_impact_distribution_{label}.png')
    plt.savefig(f'{output_folder}/gaus_impact_distribution_{label}.pdf')
    plt.show(block=False)

def angular_scan_plot(theta_x_in, theta_x_out, theta_L, popt=None, fit_bin_width=None, label=None, angle_unit='urad'):
    scale = 1e6 if angle_unit == 'urad' else 1.0
    unit_str = r'\mu rad' if angle_unit == 'urad' else 'rad'

    dtheta_x = theta_x_out - theta_x_in

    fig, graph = plt.subplots(1, 2, figsize=(12, 5))

    h = graph[0].hist2d(theta_x_in * scale, dtheta_x * scale, bins=100, cmap='viridis', norm=LogNorm())
    cb = fig.colorbar(h[3], ax=graph[0])
    cb.set_label('Counts')
    graph[0].set_ylim(-150, 150)
    # graph[0].axhline(0, color='gray', lw=0.8, ls='--')
    graph[0].axvline(-0.5 * theta_L * scale, color='red', lw=0.8, ls='--', label=r'$\pm 1/2 \theta_L$')
    graph[0].axvline(0.5 * theta_L * scale, color='red', lw=0.8, ls='--')
    graph[0].set_xlabel(rf'$\theta_{{in}}$ [${unit_str}$]')
    graph[0].set_ylabel(rf'$\Delta\theta = \theta_{{out}} - \theta_{{in}}$ [${unit_str}$]')
    graph[0].set_title(f'Angular scan – {label}')
    graph[0].grid(alpha=0.3)
    graph[0].legend()

    counts, bins, _ = graph[1].hist(dtheta_x * scale, bins=200, color='darkorange')
    plot_bin_width = bins[1] - bins[0]

    graph[1].set_xlim(-75, 100)
    graph[1].set_xlabel(rf'$\Delta\theta$ [${unit_str}$]')
    graph[1].set_ylabel('Counts')
    graph[1].set_title(rf'$\Delta\theta$ distribution – {label}')
    graph[1].grid(alpha=0.3)

    if popt is not None and fit_bin_width is not None:
        A, mu, sigma = popt

        A_scaled = A * (plot_bin_width / (fit_bin_width * scale))
        # Attenzione: mu e sigma dal fit sono in radianti, mentre il grafico è in urad!
        # Dobbiamo moltiplicarli per 'scale'
        mu_scaled = mu * scale
        sigma_scaled = sigma * scale

        x_curve = np.linspace(dtheta_x.min() * scale, dtheta_x.max() * scale, 1000)
        y_curve = gaussian(x_curve, A_scaled, mu_scaled, sigma_scaled)
        graph[1].plot(x_curve, y_curve, color='red', lw=2)

    fig.tight_layout()
    if label == '(survived particles)': label = 'survived'
    if label == '($θ_b =$ 50 µrad)': label = None
    plt.savefig(f'{output_folder}/angular_scan_{label}.png')
    plt.savefig(f'{output_folder}/angular_scan_{label}.pdf')
    plt.show(block=False)

def ch_footprint_plot(x_in, y_in, x_out, y_out, dtheta_x, crystal_x, crystal_y, tol_frac=0.2, bending_angle=bending_angle, distance_unit='mm'):
    distance_scale = 1e3 if distance_unit == 'mm' else 1.0

    channeled_mask = np.abs(dtheta_x - bending_angle) < (tol_frac * bending_angle)
    x_in_channeled = x_in[channeled_mask]
    y_in_channeled = y_in[channeled_mask]
    x_out_channeled = x_out[channeled_mask]
    y_out_channeled = y_out[channeled_mask]

    fig, graph = plt.subplots(2, 2, figsize=(10, 10))

    h0 = graph[0, 0].hist(x_in_channeled * distance_scale, bins=100, color='royalblue')
    graph[0, 0].axvline(-crystal_x/2 * distance_scale, color='red', ls='--', label='crystal edge')
    graph[0, 0].axvline(crystal_x/2 * distance_scale, color='red', ls='--')
    graph[0, 0].set_xlabel(f'impact x [{distance_unit}]')
    graph[0, 0].set_ylabel('Counts')
    graph[0, 0].set_title('Upstream - channeled particles')
    graph[0, 0].legend()

    h1 = graph[0, 1].hist(y_in_channeled * distance_scale, bins=100, color='royalblue')
    graph[0, 1].axvline(-crystal_y/2 * distance_scale, color='red', ls='--', label='crystal edge')
    graph[0, 1].axvline(crystal_y/2 * distance_scale, color='red', ls='--')
    graph[0, 1].set_xlabel(f'impact y [{distance_unit}]')
    graph[0, 1].set_ylabel('Counts')
    graph[0, 1].set_title('Upstream - channeled particles')
    graph[0, 1].legend()

    h2 = graph[1, 0].hist(x_out_channeled * distance_scale, bins=100, color='royalblue')
    graph[1, 0].axvline(-crystal_x/2 * distance_scale, color='red', ls='--', label='crystal edge')
    graph[1, 0].axvline(crystal_x/2 * distance_scale, color='red', ls='--')
    graph[1, 0].set_xlabel(f'outgoing x [{distance_unit}]')
    graph[1, 0].set_ylabel('Counts')
    graph[1, 0].set_title('Downstream - channeled particles')
    graph[1, 0].legend()

    h3 = graph[1, 1].hist(y_out_channeled * distance_scale, bins=100, color='royalblue')
    graph[1, 1].axvline(-crystal_y/2 * distance_scale, color='red', ls='--', label='crystal edge')
    graph[1, 1].axvline(crystal_y/2 * distance_scale, color='red', ls='--')
    graph[1, 1].set_xlabel(f'outgoing y [{distance_unit}]')
    graph[1, 1].set_ylabel('Counts')
    graph[1, 1].set_title('Downstream - channeled particles')
    graph[1, 1].legend()
    
    fig.tight_layout()
    plt.savefig(f'{output_folder}/ch_footprint.png')
    plt.savefig(f'{output_folder}/ch_footprint.pdf')
    plt.show(block=False)

def mcs_edge_scan(pos_in, dtheta, crystal_edge, plane_label='x', dtheta_label='x', bins=100, 
                  distance_unit='mm', angle_unit='urad', fit_window=None, label=''):

    scale = 1e6 if angle_unit == 'urad' else 1.0
    unit_str = r'\mu rad' if angle_unit == 'urad' else 'rad'
    distance_scale = 1e3 if distance_unit == 'mm' else 1.0
    
    pos_in_scaled = pos_in * distance_scale
    dtheta_scaled = dtheta * scale
    
    bin_edges = np.linspace(pos_in_scaled.min(), pos_in_scaled.max(), bins + 1)
    bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
    
    sigmas = []
    sigma_errs = []
    valid_pos = []
    
    for i in range(bins):
        mask = (pos_in_scaled >= bin_edges[i]) & (pos_in_scaled < bin_edges[i+1])
        dtheta_slice = dtheta_scaled[mask]
        
        if fit_window is not None:
            dtheta_slice = dtheta_slice[np.abs(dtheta_slice) < fit_window]
            
        if len(dtheta_slice) < 20:
            continue
            
        counts, edges = np.histogram(dtheta_slice, bins=50)
        centers = 0.5 * (edges[:-1] + edges[1:])
        
        A0 = counts.max()
        mu0 = np.mean(dtheta_slice)
        sigma0 = np.std(dtheta_slice)
        
        try:
            popt, pcov = curve_fit(gaussian, centers, counts, p0=[A0, mu0, sigma0], maxfev=2000)
            sig = abs(popt[2])
            err = np.sqrt(pcov[2, 2])
            
            # Filter out wildly failed fits
            if sig < 1000 and err < sig: 
                sigmas.append(sig)
                sigma_errs.append(err)
                valid_pos.append(bin_centers[i])
        except RuntimeError:
            # Fallback to standard deviation if fit fails
            sigmas.append(sigma0)
            sigma_errs.append(0.0)
            valid_pos.append(bin_centers[i])

    fig, ax = plt.subplots(figsize=(8, 5))
    
    ax.errorbar(valid_pos, sigmas, yerr=sigma_errs, fmt='o', color='royalblue', 
                markersize=4, capsize=3, label=rf'$\sigma$ of $\Delta\theta_{dtheta_label}$')
    
    edge_scaled = crystal_edge / 2 * distance_scale
    ax.axvline(-edge_scaled, color='red', ls='--', label='Geometrical Crystal Edges')
    ax.axvline(edge_scaled, color='red', ls='--')
    
    ax.set_xlabel(f'Impact {plane_label} [{distance_unit}]')
    ax.set_ylabel(rf'MCS width $\sigma(\Delta\theta_{dtheta_label})$ [${unit_str}$]')
    
    title = f'Material Edge Transition {label}'
    if fit_window is not None:
        title += f' (Fit Window $\pm${fit_window} {unit_str})'
    ax.set_title(title)
    
    ax.grid(alpha=0.3)
    ax.legend()
    fig.tight_layout()
    
    label_str = f"_{label}" if label else ""
    plt.savefig(f'{output_folder}/mcs_edge_scan_{plane_label}{label_str}.png')
    plt.savefig(f'{output_folder}/mcs_edge_scan_{plane_label}{label_str}.pdf')
    plt.show(block=False)
    
    return np.array(valid_pos), np.array(sigmas), np.array(sigma_errs)

def apply_selection(data_dict, cut):
    # apply same boolean mask to the arrrays of the dictionary
    return {k: v[cut] for k, v in data_dict.items()}

def gaussian(x, A, mu, sigma):
    return A * np.exp(-0.5 * ((x - mu) / sigma) ** 2)

def channeling_efficiency(theta_x_in, dtheta_x, bending_angle=bending_angle, tol_frac=0.4, angle_unit='urad'):
    scale = 1e6 if angle_unit == 'urad' else 1.0
    unit_str = r'\mu rad' if angle_unit == 'urad' else 'rad'

    n_tot = len(theta_x_in)

    counts, edges = np.histogram(dtheta_x, bins=200)
    centers = 0.5 * (edges[:-1] + edges[1:])
    bin_width = edges[1] - edges[0]

    search_window = 0.5 * abs(bending_angle)
    search_mask = np.abs(centers - bending_angle) < search_window

    if search_mask.sum() > 0:
        peak_x = centers[search_mask][np.argmax(counts[search_mask])]
    else:
        peak_x = bending_angle

    window = 0.5 * abs(bending_angle)
    mask = (centers >= (peak_x - 10/scale)) & (centers < (peak_x + window))
    while mask.sum() < 5 and window < (dtheta_x.max() - dtheta_x.min()):
        window *= 2
        mask = np.abs(centers - bending_angle) < window
 
    x_fit, y_fit = centers[mask], counts[mask]
    A0 = y_fit.max() if len(y_fit) else 1.0
    sigma0 = window / 2
 
    popt, pcov = curve_fit(gaussian, x_fit, y_fit, p0=[A0, bending_angle, sigma0], maxfev=10000)
 
    A, mu, sigma = popt
    sigma = abs(sigma)
    perr = np.sqrt(np.diag(pcov))
    mu_err = perr[1]
 
    n_channeled = A * sigma * np.sqrt(2 * np.pi) / bin_width
    efficiency = n_channeled / n_tot if n_tot > 0 else np.nan

    n_2sigma = np.sum(np.abs(dtheta_x - mu) < 2 * sigma)
    n_above  = np.sum(dtheta_x > (mu - 2 * sigma))

    print('-'*50)
    print(f"Channeling efficiency: {n_channeled:.0f}/{n_tot} = {efficiency*100:.1f}%")
    print(f"Channeling efficiency (within 2sigma): {n_2sigma:.0f}/{n_tot} = {n_2sigma/n_tot*100:.1f}%")
    print(f"Channeling efficiency (above 2sigma): {n_above:.0f}/{n_tot} = {n_above/n_tot*100:.1f}%")
    print(f"Channeling peak = {mu * scale:.2f} ± {mu_err * scale:.2f} urad - sigma = ±")

    return efficiency, popt, pcov, bin_width

# ANALYSIS

beam_distribution_plot(x_in, y_in, theta_x_in, theta_y_in, label='- all particles')
beam_distribution_plot(x_in_survived, y_in_survived, theta_x_in_survived, theta_y_in_survived, label='- survived particles')

ch_footprint_plot(x_in_survived, y_in_survived, x_out_survived, y_out_survived, dtheta_x, crystal_x, crystal_y, 0.2, bending_angle)

pos_x, sigmas_x, errs_x = mcs_edge_scan(
    pos_in=survived_data['x_in'], 
    dtheta=survived_data['dtheta_x'], 
    crystal_edge=crystal_x, 
    plane_label='x', 
    dtheta_label='x',
    bins=120, 
    fit_window=20, # urad window around 0
    label='- X Plane'
)
pos_y, sigmas_y, errs_y = mcs_edge_scan(
    pos_in=survived_data['y_in'], 
    dtheta=dtheta_y, 
    crystal_edge=crystal_y, 
    plane_label='y', 
    dtheta_label='y',
    bins=120, 
    fit_window=None, 
    label='- Y Plane'
)

# geometric cut
geometric_hit = (np.abs(survived_data['x_in']) < crystal_x/2) & (np.abs(survived_data['y_in']) < crystal_y/2)
geom_data = apply_selection(survived_data, geometric_hit)
print(f"Sopravvissute che hanno colpito geometricamente il cristallo: {geometric_hit.sum()}")
print(f"Sopravvissute che hanno mancato il cristallo: {(~geometric_hit).sum()}")
beam_distribution_plot(geom_data['x_in'], geom_data['y_in'], geom_data['theta_x_in'], geom_data['theta_y_in'], label='- entered particles')

# lindhard selection
lindhard_cut = np.abs(geom_data['theta_x_in']) < (0.5 * theta_L1)
final_data = apply_selection(geom_data, lindhard_cut)

# efficiency
eff, popt, pcov, fit_bin_width = channeling_efficiency(final_data['theta_x_in'], final_data['dtheta_x'])

angular_scan_plot(final_data['theta_x_in'], final_data['theta_x_out'], theta_L1, label='($θ_b =$ 50 µrad)', popt=popt, fit_bin_width=fit_bin_width)
angular_scan_plot(survived_data['theta_x_in'], survived_data['theta_x_out'], theta_L1, label='(survived particles)', popt=popt, fit_bin_width=fit_bin_width)

