import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

filename = 'cry1_20260827_18.npz'
data = np.load(filename)

# reconstruct the arrays from the saved data
x_in = data['x_in']
y_in = data['y_in']
py_in = data['py_in']
px_in = data['px_in']
theta_out_all = data['px_out']
particle_id_in = data['particle_id_in']
particle_id_out = data['particle_id_out']
state_out = data['state_out']
bending_angle = float(data['bending_angle'])
label = str(data['label'])

survived = state_out == 1
ids_out = particle_id_out[survived]
theta_out_x = theta_out_all[survived]

order_in = np.argsort(particle_id_in)
sorted_ids_in = particle_id_in[order_in]
idx = np.searchsorted(sorted_ids_in, ids_out)

x_in_survived = data['x_in'][order_in][idx]
y_in_survived = data['y_in'][order_in][idx]
theta_in_x = px_in[order_in][idx]
px_in_survived = theta_in_x
py_in_survived = data['py_in'][order_in][idx]

dtheta = theta_out_x - theta_in_x

theta_L1 = 13.228e-6  # rad Lindhard critical angle for Si(110) at 180 GeV
theta_L2 = 12.992e-6  # rad Lindhard critical angle for Si(110) at 180 GeV

def beam_distribution_plot(x_in, y_in, px_in, py_in, label=None, distance_unit='mm', angle_unit='urad'):
    scale = 1e6 if angle_unit == 'urad' else 1.0
    unit_str = r'\mu rad' if angle_unit == 'urad' else 'rad'

    distance_scale = 1e3 if distance_unit == 'mm' else 1.0

    x_lim = (x_in.min() * distance_scale, x_in.max() * distance_scale)
    y_lim = (y_in.min() * distance_scale, y_in.max() * distance_scale)
    px_lim = (-150, 150) if angle_unit == 'urad' else (px_in.min() * scale, px_in.max() * scale)
    py_lim = (-150, 150) if angle_unit == 'urad' else (py_in.min() * scale, py_in.max() * scale)
    # px_lim = (px_in.min() * scale, px_in.max() * scale)
    # py_lim = (py_in.min() * scale, py_in.max() * scale)

    fig, graph = plt.subplots(1, 3, figsize=(15, 5))
    # 1. Grafico (x, y) - Beam impact position
    h0 = graph[0].hist2d(x_in * distance_scale, y_in * distance_scale, bins=200, cmap='viridis', norm=LogNorm())
    fig.colorbar(h0[3], ax=graph[0], label='Counts')
    graph[0].set_xlabel(f'impact x [{distance_unit}]')
    graph[0].set_ylabel(f'impact y [{distance_unit}]')
    graph[0].set_title(f'Beam Impact Position (x, y) {label}')
    graph[0].set_xlim(x_lim)
    graph[0].set_ylim(y_lim)
    graph[0].grid(alpha=0.3)
    # 2. Grafico (x, px) - Phase space orizzontale (scala px in urad per coerenza)
    h1 = graph[1].hist2d(x_in * distance_scale, px_in * scale, bins=200, cmap='viridis', norm=LogNorm())
    fig.colorbar(h1[3], ax=graph[1], label='Counts')
    graph[1].set_xlabel(f'impact x [{distance_unit}]')
    graph[1].set_ylabel(rf'$px_{{in}}$ [${unit_str}$]')
    graph[1].set_title(f'Horizontal Phase Space (x, px) {label}')
    graph[1].set_xlim(x_lim)
    graph[1].set_ylim(px_lim)
    graph[1].grid(alpha=0.3)
    # 3. Grafico (y, py) - Phase space verticale (scala py in urad)
    h2 = graph[2].hist2d(y_in * distance_scale, py_in * scale, bins=200, cmap='viridis', norm=LogNorm())
    fig.colorbar(h2[3], ax=graph[2], label='Counts')
    graph[2].set_xlabel(f'impact y [{distance_unit}]')
    graph[2].set_ylabel(rf'$py_{{in}}$ [${unit_str}$]')
    graph[2].set_title(f'Vertical Phase Space (y, py) {label}')
    graph[2].set_xlim(y_lim)
    graph[2].set_ylim(py_lim)
    graph[2].grid(alpha=0.3)

    fig.tight_layout()
    plt.show(block=False)
    plt.pause(0.1)

def angular_scan_plot(px_in, particle_id_in, px_out, particle_id_out, state_out,
                       theta_L, bending_angle, label, angle_unit='urad'):
    scale = 1e6 if angle_unit == 'urad' else 1.0
    unit_str = r'\mu rad' if angle_unit == 'urad' else 'rad'

    survived = state_out == 1
    ids_out = particle_id_out[survived]
    theta_out_x = px_out[survived]

    # match via particle_id
    order_in = np.argsort(particle_id_in)
    sorted_ids_in = particle_id_in[order_in]
    idx = np.searchsorted(sorted_ids_in, ids_out)
    theta_in_x = px_in[order_in][idx]
    dtheta = theta_out_x - theta_in_x

    # cutting on the acceptance of the crystal
    cut = np.abs(theta_in_x) < theta_L
    theta_in_cut = theta_in_x[cut]
    dtheta_cut = dtheta[cut]

    fig, graph = plt.subplots(1, 2, figsize=(12, 5))

    h = graph[0].hist2d(theta_in_x * scale, dtheta * scale, bins=200, cmap='viridis', norm=LogNorm())
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

    graph[1].hist(dtheta * scale, bins=500, color='#1f77b4', edgecolor='white', linewidth=0.3)
    graph[1].set_xlim(-75, 100)
    graph[1].set_xlabel(rf'$\Delta\theta$ [${unit_str}$]')
    graph[1].set_ylabel('Counts')
    graph[1].set_title(rf'$\Delta\theta$ distribution – {label}')
    graph[1].grid(alpha=0.3)

    fig.tight_layout()
    plt.show(block=False)
    plt.pause(0.1)

    return theta_in_x, theta_in_cut, dtheta, dtheta_cut

def channeling_efficiency(theta_in_x, dtheta, bending_angle, theta_L, tol_frac=0.2):

    # cutting within the critical angle
    in_acceptance = np.abs(theta_in_x) < theta_L
    n_accepted = np.sum(in_acceptance)

    tol = tol_frac * bending_angle
    channeled = in_acceptance & (np.abs(dtheta - bending_angle) < tol)
    n_channeled = np.sum(channeled)

    efficiency = n_channeled / n_accepted if n_accepted > 0 else np.nan
    print(f"Particles within critical angle: {n_accepted}")
    print(f"Channeling particles: {n_channeled}")
    print(f"Channeling efficiency: {efficiency*100:.1f}%")
    return efficiency

beam_distribution_plot(x_in, y_in, px_in, py_in, label='- all particles')
beam_distribution_plot(x_in_survived, y_in_survived, px_in_survived, py_in_survived, label='- survived particles')

theta_in_1, theta_in_1_cut, dtheta_1, dtheta_1_cut = angular_scan_plot(px_in=data['px_in'], particle_id_in=data['particle_id_in'],
    px_out=data['px_out'], particle_id_out=data['particle_id_out'], state_out=data['state_out'], theta_L=theta_L1,
    bending_angle=float(data['bending_angle']), label='($θ_b = 50 µrad)')
# theta_in_2, theta_in_2_cut, dtheta_2, dtheta_2_cut = angular_scan_plot(
#     px_in=data['px_in'], particle_id_in=data['particle_id_in'],
#     px_out=data['px_out'], particle_id_out=data['particle_id_out'],
#     state_out=data['state_out'],
#     theta_L=theta_L2, bending_angle=float(data['bending_angle']),
#     label='cry2 (bending 6 mrad)')

eff1 = channeling_efficiency(theta_in_1, dtheta_1, bending_angle=50e-6, theta_L=theta_L1)
# eff2 = channeling_efficiency(theta_in_2, dtheta_2, bending_angle=6e-3, theta_L=theta_L2)

# exact_zero = dtheta == 0.0
# near_zero = np.abs(dtheta) < 1e-9  # soglia numerica minima

# print(f"Particles with dtheta exactly 0: {exact_zero.sum()}")
# print(f"Particles with dtheta < 1e-9: {near_zero.sum()}")

# suspect = dtheta == 0.0   # o near_zero, a seconda di cosa trovi sopra

# print("=== Particelle sospette (dtheta=0) ===")
# print(f"N: {suspect.sum()} su {len(dtheta)} ({100*suspect.sum()/len(dtheta):.1f}%)")
# print(f"x_in range: [{theta_in_x[suspect].min():.2e}, {theta_in_x[suspect].max():.2e}]")  # occhio: qui serve x, non theta_in_x
# print(f"state_out valori unici: {np.unique(state_out[suspect])}")

# # posizione d'impatto x,y per le sospette vs il resto
# fig, graph = plt.subplots(1, 2, figsize=(12,5))
# graph[0].hist(x_impact[suspect]*1e3, bins=50, alpha=0.6, label='dtheta=0', density=True)
# graph[0].hist(x_impact[~suspect]*1e3, bins=50, alpha=0.6, label='resto', density=True)
# graph[0].set_xlabel('x impact [mm]')
# graph[0].legend()

# graph[1].hist(y_impact[suspect]*1e3, bins=50, alpha=0.6, label='dtheta=0', density=True)
# graph[1].hist(y_impact[~suspect]*1e3, bins=50, alpha=0.6, label='resto', density=True)
# graph[1].set_xlabel('y impact [mm]')
# graph[1].legend()
# plt.show()