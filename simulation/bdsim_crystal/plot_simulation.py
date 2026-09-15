import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

filename = 'cry1_20260915_16.npz'
data = np.load(filename)

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
    h0 = graph[0, 0].hist2d(x_in * distance_scale, y_in * distance_scale, bins=200, cmap='viridis', norm=LogNorm())
    fig.colorbar(h0[3], ax=graph[0, 0], label='Counts')
    graph[0, 0].set_xlabel(f'impact x [{distance_unit}]')
    graph[0, 0].set_ylabel(f'impact y [{distance_unit}]')
    graph[0, 0].set_title(f'Beam Impact Position (x, y) {label}')
    graph[0, 0].set_xlim(x_lim)
    graph[0, 0].set_ylim(y_lim)
    graph[0, 0].grid(alpha=0.3)

    h00 = graph[0, 1].hist2d(theta_x_in * scale, theta_y_in * scale, bins=200, cmap='viridis', norm=LogNorm())
    fig.colorbar(h00[3], ax=graph[0, 1], label='Counts')
    graph[0, 1].set_xlabel(rf'$\theta_x$ [${unit_str}$]')
    graph[0, 1].set_ylabel(rf'$\theta_y$ [${unit_str}$]')
    graph[0, 1].set_title(f'Angular Distribution {label}')
    graph[0, 1].set_xlim(theta_x_lim)
    graph[0, 1].set_ylim(theta_y_lim)
    graph[0, 1].grid(alpha=0.3)

    h1 = graph[1, 0].hist2d(x_in * distance_scale, theta_x_in * scale, bins=200, cmap='viridis', norm=LogNorm())
    fig.colorbar(h1[3], ax=graph[1, 0], label='Counts')
    graph[1, 0].set_xlabel(f'impact x [{distance_unit}]')
    graph[1, 0].set_ylabel(rf'$\theta_x$ [${unit_str}$]')
    graph[1, 0].set_title(f'Horizontal Phase Space (x, $\theta_x$) {label}')
    graph[1, 0].set_xlim(x_lim)
    graph[1, 0].set_ylim(theta_x_lim)
    graph[1, 0].grid(alpha=0.3)

    h2 = graph[1, 1].hist2d(y_in * distance_scale, theta_y_in * scale, bins=200, cmap='viridis', norm=LogNorm())
    fig.colorbar(h2[3], ax=graph[1, 1], label='Counts')
    graph[1, 1].set_xlabel(f'impact y [{distance_unit}]')
    graph[1, 1].set_ylabel(rf'$\theta_y$ [${unit_str}$]')
    graph[1, 1].set_title(f'Vertical Phase Space (y, $\theta_y$) {label}')
    graph[1, 1].set_xlim(y_lim)
    graph[1, 1].set_ylim(theta_y_lim)
    graph[1, 1].grid(alpha=0.3)

    fig.tight_layout()
    plt.show(block=False)
    plt.pause(0.1)

def angular_scan_plot(theta_x_in, theta_x_out, theta_L, label, angle_unit='urad'):
    scale = 1e6 if angle_unit == 'urad' else 1.0
    unit_str = r'\mu rad' if angle_unit == 'urad' else 'rad'

    dtheta_x = theta_x_out - theta_x_in

    fig, graph = plt.subplots(1, 2, figsize=(12, 5))

    h = graph[0].hist2d(theta_x_in * scale, dtheta_x * scale, bins=200, cmap='viridis', norm=LogNorm())
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

    graph[1].hist(dtheta_x * scale, bins=500, color='#1f77b4', edgecolor='white', linewidth=0.3)
    graph[1].set_xlim(-75, 100)
    graph[1].set_xlabel(rf'$\Delta\theta$ [${unit_str}$]')
    graph[1].set_ylabel('Counts')
    graph[1].set_title(rf'$\Delta\theta$ distribution – {label}')
    graph[1].grid(alpha=0.3)

    fig.tight_layout()
    plt.show(block=False)

def apply_selection(data_dict, cut):
    # apply same boolean mask to the arrrays of the dictionary
    return {k: v[cut] for k, v in data_dict.items()}

def channeling_efficiency(theta_x_in, dtheta_x, bending_angle, tol_frac=0.2, angle_unit='urad'):
    scale = 1e6 if angle_unit == 'urad' else 1.0
    unit_str = r'\mu rad' if angle_unit == 'urad' else 'rad'

    n_tot = len(theta_x_in)

    tol = tol_frac * bending_angle
    channeled = (np.abs(dtheta_x - bending_angle) < tol)
    n_channeled = np.sum(channeled)

    efficiency = n_channeled / n_tot if n_tot > 0 else np.nan
    print(f"Particles within critical angle: {n_tot}")
    print(f"Channeling particles: {n_channeled}")
    print(f"Channeling efficiency: {efficiency*100:.1f}%")
    return efficiency

# ANALYSIS

beam_distribution_plot(x_in, y_in, theta_x_in, theta_y_in, label='- all particles')
beam_distribution_plot(x_in_survived, y_in_survived, theta_x_in_survived, theta_y_in_survived, label='- survived particles')

angular_scan_plot(theta_x_in_survived, theta_x_out_survived, theta_L1, label='($θ_b =$ 50 µrad)')

# geometric cut
geometric_hit = (np.abs(survived_data['x_in']) < crystal_x/2) & (np.abs(survived_data['y_in']) < crystal_y/2)
geom_data = apply_selection(survived_data, geometric_hit)
print(f"Sopravvissute che hanno colpito geometricamente il cristallo: {geometric_hit.sum()}")
print(f"Sopravvissute che hanno mancato il cristallo: {(~geometric_hit).sum()}")

# lindhard selection
lindhard_cut = np.abs(geom_data['theta_x_in']) < (0.5 * theta_L1)
final_data = apply_selection(geom_data, lindhard_cut)

# efficiency
eff = channeling_efficiency(final_data['theta_x_in'], final_data['dtheta_x'], bending_angle=50e-6)
