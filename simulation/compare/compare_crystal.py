import sys
import numpy as np
import matplotlib.pyplot as plt

from scipy.optimize import curve_fit

def gauss(x, A, mu, s):
    return A * np.exp(-0.5 * ((x - mu) / s) ** 2)

def fit_peak(arr, lo=5950, hi=6050, bins=100):
    # arr in urad; ritorna (mu, sigma, errore su mu)
    h, e = np.histogram(arr, bins=bins, range=(lo, hi))
    c = 0.5 * (e[:-1] + e[1:])
    p, cov = curve_fit(gauss, c, h, p0=[h.max(), c[np.argmax(h)], 5])
    return p[1], abs(p[2]), np.sqrt(cov[1, 1])

plt.rcParams.update({
    'font.size': 16,          # Dimensione generale del font
    'axes.titlesize': 18,     # Dimensione dei titoli dei singoli assi/subplot
    'axes.labelsize': 16,     # Dimensione delle etichette (asse x, y)
    'xtick.labelsize': 14,    # Dimensione dei numeri sull'asse x
    'ytick.labelsize': 14,    # Dimensione dei numeri sull'asse y
    'legend.fontsize': 14,    # Dimensione del testo nella legenda
    'figure.titlesize': 20    # Dimensione del suptitle (titolo globale della figura)
})


CRYSTAL = sys.argv[1] if len(sys.argv) > 1 else 'long'

if CRYSTAL == 'long':
    file1 = '../bdsim_longcrystal/plots6/dtheta_x_final_smeared_bdsim.npy'
    file2 = '../xcoll_longcrystal/plots/dtheta_x_final_smeared_xcoll.npy' 
    xrange = (-750, 6200)
    xrange_zoom = (5950, 6050)
    bin_num = 6950
    eff1 = (0.224, 0.001)
    eff2 = (0.276, 0.001)
    color1 = "#982916"
    color2 = "#68238B"
elif CRYSTAL == 'short':
    file1 = '../bdsim_shortcrystal/plots/dtheta_x_final_smeared_bdsim.npy'
    file2 = '../xcoll_shortcrystal/plots/dtheta_x_final_smeared_xcoll.npy' 
    xrange = (-100, 120)
    xrange_zoom = (20, 80)
    bin_num = 440
    eff1 = (0.633, 0.005)
    eff2 = (0.635, 0.002)
    color1 = "#0072B2"
    color2 = "#228833"

def comparison_plot(dtheta_1, dtheta_2, xrange, xrange_zoom, label1='Long 1', label2='Long 2', angle_unit='urad', eff1=None, eff2=None, add_data=False):
    if CRYSTAL == 'short': add_data = False

    scale = 1e6 if angle_unit == 'urad' else 1.0
    unit_str = r'\mu rad' if angle_unit == 'urad' else 'rad'

    fig, graph = plt.subplots(1, 2, figsize=(12, 5))

    graph[0].hist(dtheta_1 * scale, bins=200, range=xrange, 
                  color=color1, alpha=1, label=label1, density=True, rasterized=True)
    graph[0].hist(dtheta_2 * scale, bins=200, range=xrange, 
                  color=color2, alpha=0.7, label=label2, density=True, rasterized=True)

    if add_data:
        data_file = f'../../recoDataSimple/plots/plots_8430_py/dtheta_data.npy'
        dtheta_data = np.load(data_file)
        y, edges = np.histogram(dtheta_data, bins=200, range=xrange, density=True)
        graph[0].plot(0.5 * (edges[:-1] + edges[1:]), y, 'ok', ms=3, alpha=1, label='Data')


    if CRYSTAL == 'long':
        graph[0].set_yscale('log')
        graph[0].set_ylim(3e-6, 1e-2)
    graph[0].set_xlabel(rf'$\Delta\theta_x$ [${unit_str}$]')
    graph[0].set_ylabel('Normalized Counts')
    graph[0].set_title(rf'$\Delta\theta_x$ distribution comparison')
    graph[0].grid(alpha=0.3)
    graph[0].legend()

    graph[1].hist(dtheta_1 * scale, bins=bin_num, range=xrange, 
                  color=color1, alpha=1, label=label1, density=True, rasterized=True)
    graph[1].hist(dtheta_2 * scale, bins=bin_num, range=xrange, 
                  color=color2, alpha=0.7, label=label2, density=True, rasterized=True)
    if add_data:
        y, edges = np.histogram(dtheta_data, bins=bin_num, range=xrange, density=True)
        graph[1].plot(0.5 * (edges[:-1] + edges[1:]), y, 'ok', ms=3, alpha=0.8, label='Data')

    graph[1].set_xlabel(rf'$\Delta\theta_x$ [${unit_str}$]')
    graph[1].set_ylabel('Normalized Counts')
    graph[1].set_xlim(xrange_zoom[0], xrange_zoom[1])
    graph[1].set_title(rf'$\Delta\theta_x$ distribution (zoom)')
    graph[1].grid(alpha=0.3)
    if eff1 is not None and eff2 is not None:
            txt = (rf'$\epsilon$ ({label1[-5:]}) = ({eff1[0]*100:.1f} ± {eff1[1]*100:.1f}) %' '\n'
                   rf'$\epsilon$ ({label2[-5:]}) = ({eff2[0]*100:.1f} ± {eff2[1]*100:.1f}) %')
            graph[1].text(0.97, 0.95, txt, transform=graph[1].transAxes,
                          ha='right', va='top', fontsize=12,
                          bbox=dict(boxstyle='round', facecolor='white', edgecolor='gray', alpha=0.85))

    fig.tight_layout()
    plt.savefig(f'{CRYSTAL}_comparison_dtheta_x.png', dpi=300)
    plt.savefig(f'{CRYSTAL}_comparison_dtheta_x.pdf')
    plt.close(fig)
    print(f"Plot saved as '{CRYSTAL}_comparison_dtheta_x.png'")

    def peak_stats(arr, mu, half=40):
        # arr in urad; frazione di eventi nel picco rispetto al totale nel range
        tot = np.sum((arr > xrange[0]) & (arr < xrange[1]))
        in_peak = np.sum(np.abs(arr - mu) < half)
        return in_peak / tot

    for name, arr, sc in [('BDSIM', dtheta_1, 1e6), ('XColl', dtheta_2, 1e6), ('Data', dtheta_data, 1)]:
        a = arr * sc
        mu, s, err = fit_peak(a)
        print(f'{name}: sigma = {s:.2f} urad, frazione nel picco = {peak_stats(a, mu):.3f}')


dtheta_x_1 = np.load(file1)
dtheta_x_2 = np.load(file2)


comparison_plot(dtheta_x_1, dtheta_x_2, xrange, xrange_zoom, label1='Geant4/BDSIM', label2='Xsuite/XColl', eff1=eff1, eff2=eff2, add_data=True)





