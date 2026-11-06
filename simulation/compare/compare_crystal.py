import sys
import numpy as np
import matplotlib.pyplot as plt

viridis_blue = "#3b528b"
viridis_green = "#5ec962"

CRYSTAL = sys.argv[1] if len(sys.argv) > 1 else 'long'

if CRYSTAL == 'long':
    file1 = '../bdsim_longcrystal/plots6/dtheta_x_final_smeared_bdsim.npy'
    file2 = '../xcoll_longcrystal/plots/dtheta_x_final_smeared_xcoll.npy' 
    range = (-750, 6200)
    range_zoom = (5950, 6050)
    bin_num = 6950
elif CRYSTAL == 'short':
    file1 = '../bdsim_shortcrystal/plots/dtheta_x_final_smeared_bdsim.npy'
    file2 = '../xcoll_shortcrystal/plots/dtheta_x_final_smeared_xcoll.npy' 
    range = (-100, 120)
    range_zoom = (20, 80)
    bin_num = 440

def comparison_plot(dtheta_1, dtheta_2, range, range_zoom, label1='Long 1', label2='Long 2', angle_unit='urad'):
    scale = 1e6 if angle_unit == 'urad' else 1.0
    unit_str = r'\mu rad' if angle_unit == 'urad' else 'rad'

    fig, graph = plt.subplots(1, 2, figsize=(12, 5))

    graph[0].hist(dtheta_1 * scale, bins=200, range=range, 
                  color=viridis_blue, alpha=1, label=label1, density=True, rasterized=True)
    graph[0].hist(dtheta_2 * scale, bins=200, range=range, 
                  color=viridis_green, alpha=0.7, label=label2, density=True, rasterized=True)

    if CRYSTAL == 'long':
        graph[0].set_yscale('log')
    graph[0].set_xlabel(rf'$\Delta\theta_x$ [${unit_str}$]')
    graph[0].set_ylabel('Normalized Counts')
    graph[0].set_title(rf'$\Delta\theta_x$ distribution comparison')
    graph[0].grid(alpha=0.3)
    graph[0].legend()

    graph[1].hist(dtheta_1 * scale, bins=bin_num, range=range, 
                  color=viridis_blue, alpha=1, label=label1, density=True, rasterized=True)
    graph[1].hist(dtheta_2 * scale, bins=bin_num, range=range, 
                  color=viridis_green, alpha=0.7, label=label2, density=True, rasterized=True)

    graph[1].set_xlabel(rf'$\Delta\theta_x$ [${unit_str}$]')
    graph[1].set_ylabel('Normalized Counts')
    graph[1].set_xlim(range_zoom[0], range_zoom[1])
    graph[1].set_title(rf'$\Delta\theta_x$ distribution (zoom)')
    graph[1].grid(alpha=0.3)
    graph[1].legend()

    fig.tight_layout()
    plt.savefig(f'{CRYSTAL}_comparison_dtheta_x.png', dpi=300)
    plt.savefig(f'{CRYSTAL}_comparison_dtheta_x.pdf')
    plt.close(fig)
    print(f"Plot saved as '{CRYSTAL}_comparison_dtheta_x.png'")


dtheta_x_1 = np.load(file1)
dtheta_x_2 = np.load(file2)

comparison_plot(dtheta_x_1, dtheta_x_2, range, range_zoom, label1='Geant4/BDSIM', label2='Xsuite/XColl')



