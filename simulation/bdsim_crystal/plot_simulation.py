import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

filename = 'cry1_20260827_18.npz'
data = np.load(filename)

# reconstruct the arrays from the saved data
theta_in_all = data['px_in']
theta_out_all = data['px_out']
particle_id_in = data['particle_id_in']
particle_id_out = data['particle_id_out']
state_out = data['state_out']
bending_angle = float(data['bending_angle'])
label = str(data['label'])

survived = state_out == 1
ids_out = particle_id_out[survived]
theta_out = theta_out_all[survived]

order_in = np.argsort(particle_id_in)
sorted_ids_in = particle_id_in[order_in]
idx = np.searchsorted(sorted_ids_in, ids_out)
theta_in = theta_in_all[order_in][idx]

dtheta = theta_out - theta_in

theta_L1 = 13.228e-6  # rad Lindhard critical angle for Si(110) at 180 GeV
theta_L2 = 12.992e-6  # rad Lindhard critical angle for Si(110) at 180 GeV

def angular_scan_plot(px_in, particle_id_in, px_out, particle_id_out, state_out,
                       theta_L, bending_angle, label, angle_unit='urad'):
    scale = 1e6 if angle_unit == 'urad' else 1.0
    unit_str = r'\mu rad' if angle_unit == 'urad' else 'rad'

    survived = state_out == 1
    ids_out = particle_id_out[survived]
    theta_out = px_out[survived]

    # match via particle_id
    order_in = np.argsort(particle_id_in)
    sorted_ids_in = particle_id_in[order_in]
    idx = np.searchsorted(sorted_ids_in, ids_out)
    theta_in = px_in[order_in][idx]
    dtheta = theta_out - theta_in

    # cutting on the acceptance of the crystal
    cut = np.abs(theta_in) < theta_L
    theta_in_cut = theta_in[cut]
    dtheta_cut = dtheta[cut]

    fig, graph = plt.subplots(1, 2, figsize=(12, 5))

    h = graph[0].hist2d(theta_in * scale, dtheta * scale, bins=200, cmap='viridis', norm=LogNorm())
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

    return theta_in, theta_in_cut, dtheta, dtheta_cut

def channeling_efficiency(theta_in, dtheta, bending_angle, theta_L, tol_frac=0.2):

    # cutting within the critical angle
    in_acceptance = np.abs(theta_in) < theta_L
    n_accepted = np.sum(in_acceptance)

    tol = tol_frac * bending_angle
    channeled = in_acceptance & (np.abs(dtheta - bending_angle) < tol)
    n_channeled = np.sum(channeled)

    efficiency = n_channeled / n_accepted if n_accepted > 0 else np.nan
    print(f"Particelle in accettanza: {n_accepted}")
    print(f"Particelle channelate: {n_channeled}")
    print(f"Efficienza di channeling: {efficiency*100:.1f}%")
    return efficiency

theta_in_1, theta_in_1_cut, dtheta_1, dtheta_1_cut = angular_scan_plot(
    px_in=data['px_in'], particle_id_in=data['particle_id_in'],
    px_out=data['px_out'], particle_id_out=data['particle_id_out'],
    state_out=data['state_out'],
    theta_L=theta_L1, bending_angle=float(data['bending_angle']),
    label='cry1 (bending 50 µrad)')
# theta_in_2, theta_in_2_cut, dtheta_2, dtheta_2_cut = angular_scan_plot(
#     px_in=data['px_in'], particle_id_in=data['particle_id_in'],
#     px_out=data['px_out'], particle_id_out=data['particle_id_out'],
#     state_out=data['state_out'],
#     theta_L=theta_L2, bending_angle=float(data['bending_angle']),
#     label='cry2 (bending 6 mrad)')

eff1 = channeling_efficiency(theta_in_1, dtheta_1, bending_angle=50e-6, theta_L=theta_L1)
# eff2 = channeling_efficiency(theta_in_2, dtheta_2, bending_angle=6e-3, theta_L=theta_L2)