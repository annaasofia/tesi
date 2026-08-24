import numpy as np
import xtrack as xt
import bdsim
import xpart as xp
import matplotlib.pyplot as plt

drift = xt.Drift(length=1)

line = xt.Line(elements=[drift])
line.particle_ref = xt.Particles(mass0=xt.PROTON_MASS_EV, q0=1, energy0=180e9)

bds_link = bdsim.BDSLinkTrackerInterface.GetInstance("../bdsim_crystal/trackerInterface.gmad",
                                                     referenceParticlePDG=2212,
                                                     referenceKineticEnergy=180e3, # MeV
                                                    relativeEnergyCut=0.1,
                                                    batchMode = True,
                                                    seed = 82)
l = bds_link.GetBDSIMLink()

# short crystal 2 mm
cry1 = bdsim.Element()
cry1.type = bdsim.elementtype.ElementType.USERCOMPONENT
cry1.userTypeName = "crystaldeflector"
cry1.name = "crystaldeflector"
cry1.material="G4_Si"
cry1.xsize = 0.05
cry1.ysize = 0.2
cry1.materialThickness = 2e-3
cry1.l = 2e-3
cry1.horizontalWidth = 0.1
cry1.userParameters="crystalRegion:crystaldeflector crystalLattice:(110) crystalBendingAngle:50e-6;"
cry1.apertureType = "rectangular" # aperture of box around crystal
cry1.aper1 = 10e1 # m
l.AddLinkElement(cry1)

# long crystal 74 mm
# cry2 = bdsim.Element()
# cry2.type = bdsim.elementtype.ElementType.USERCOMPONENT
# cry2.userTypeName = "crystaldeflector"
# cry2.name = "crystaldeflector2"
# cry2.material="G4_Si"
# cry2.xsize = 0.05 # check if this is correct, in the GMAD file it is 0.05
# cry2.ysize = 0.2 # check
# cry2.materialThickness = 74e-3
# cry2.l = 74e-3
# cry2.horizontalWidth = 0.1
# cry2.userParameters="crystalRegion:crystaldeflector crystalLattice:(110) crystalBendingAngle:6e-3;"
# cry2.apertureType = "rectangular" # aperture of box around crystal
# cry2.aper1 = 10e1 # m
# l.AddLinkElement(cry2)

npart = 200000
particles1 = line.build_particles(
    nemitt_x=2.5e-6, nemitt_y=1e-6,
    x=np.zeros(npart), y=np.zeros(npart),
    # px=(np.arange(-2,2,(2-(-2))/npart)+0)*1e-6
    px=np.linspace(-40e-6, 100e-6, npart), # distribuzione piatta
    py=np.zeros(npart),
    zeta=np.zeros(npart), delta=np.zeros(npart), _capacity = int(npart*2))
particles1.pdg_id[:npart] = np.ones(npart)*2212
particles10 = particles1.copy()
bds_link.TrackXSuite(0,'crystaldeflector',particles1,180e3)

# npart = 20000
# particles2 = line.build_particles(
#     nemitt_x=2.5e-6, nemitt_y=1e-6,
#     x=np.zeros(npart), y=np.zeros(npart),
#     px=(np.arange(-2,2,(2-(-2))/npart)+0)*1e-6, py=np.zeros(npart),
#     zeta=np.zeros(npart), delta=np.zeros(npart), _capacity = int(npart*2))
# particles2.pdg_id[:npart] = np.ones(npart)*2212
# particles20 = particles2.copy()
# bds_link.TrackXSuite(0,'crystaldeflector2',particles2,180e3)

def angular_scan_plot(p_in, p_out, label, angle_unit='urad'):
    scale = 1e6 if angle_unit == 'urad' else 1.0
    unit_str = r'\mu rad' if angle_unit == 'urad' else 'rad'

    survived = p_out.state == 1
    ids_out = p_out.particle_id[survived]
    theta_out = p_out.px[survived]

    # match via particle_id
    order_in = np.argsort(p_in.particle_id)
    sorted_ids_in = p_in.particle_id[order_in]
    idx = np.searchsorted(sorted_ids_in, ids_out)
    theta_in = p_in.px[order_in][idx]

    dtheta = theta_out - theta_in

    fig, graph = plt.subplots(1, 2, figsize=(12, 5))

    # angular scan
    graph[0].scatter(theta_in * scale, dtheta * scale, s=6, alpha=0.6)
    graph[0].axhline(0, color='gray', lw=0.8, ls='--')
    graph[0].set_xlabel(rf'$\theta_{{in}}$ [${unit_str}$]')
    graph[0].set_ylabel(rf'$\Delta\theta = \theta_{{out}} - \theta_{{in}}$ [${unit_str}$]')
    graph[0].set_title(f'Angular scan – {label}')
    graph[0].grid(alpha=0.3)

    # histo 1d delta theta
    graph[1].hist(dtheta * scale, bins=100, color='#1f77b4', edgecolor='white', linewidth=0.3)
    graph[1].set_xlabel(rf'$\Delta\theta$ [${unit_str}$]')
    graph[1].set_ylabel('Counts')
    graph[1].set_title(f'$\Delta\\theta$ distribution – {label}')
    graph[1].grid(alpha=0.3)

    fig.tight_layout()
    plt.show()

    return theta_in, dtheta

theta_in_1, dtheta_1 = angular_scan_plot(particles10, particles1, 'cry1 (bending 50 µrad)')
theta_in_2, dtheta_2 = angular_scan_plot(particles20, particles2, 'cry2 (bending 6 mrad)')

def channeling_efficiency(theta_in, dtheta, bending_angle, theta_L, tol_frac=0.2):

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

theta_L1 = 13.228e-6  # rad Lindhard critical angle for Si(110) at 180 GeV
eff1 = channeling_efficiency(theta_in_1, dtheta_1, bending_angle=50e-6, theta_L=theta_L1)

theta_L2 = 12.992e-6  # rad Lindhard critical angle for Si(110) at 180 GeV
eff2 = channeling_efficiency(theta_in_2, dtheta_2, bending_angle=6e-3, theta_L=theta_L2)