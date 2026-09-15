import numpy as np
import xtrack as xt
import bdsim
import xpart as xp
from datetime import datetime
from scipy.stats import norm
import simulation_functions as simfun


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
cry1.xsize = 2e-3 # meters
cry1.ysize = 35e-3 # meters
cry1.materialThickness = 4e-3 # meters
cry1.l = 4e-3 # meters
cry1.horizontalWidth = 0.1 # meters
cry1.userParameters="crystalRegion:crystaldeflector crystalLattice:(110) crystalBendingAngle:50e-6;"
cry1.apertureType = "rectangular" # aperture of box around crystal
cry1.aper1 = 0.05 # 50 mm
cry1.aper2 = 0.05 # 50 mm
l.AddLinkElement(cry1)

# x_half_range = cry1.xsize / 2
# y_half_range = cry1.ysize / 2

#correlated gaussian beam
cov_avg, eigvals, n_list = simfun.cov_matrix()
assert eigvals.min() > 0, f"cov_avg non è definita positiva! min eigenvalue = {eigvals.min()}"

mu_avg = simfun.mean_vector()  # [mu_x, mu_y, mu_px, mu_py]
rng = np.random.default_rng(seed=82)  # same seed of BDSLinkTrackerInterface, per coerenza/riproducibilità

npart = 200000
samples = rng.multivariate_normal(mean=mu_avg, cov=cov_avg, size=npart)
x_impact  = samples[:, 0]
y_impact  = samples[:, 1]
px_impact = samples[:, 2]
py_impact = samples[:, 3]
# x_impact = np.random.normal(-0.12e-3, 2.20e-3, npart)
# y_impact = np.random.normal(0.64e-3, 2.34e-3, npart)
# px_impact = np.random.normal(-0.93e-6, 26.92e-6, npart)
# py_impact = np.random.normal(6.79e-6, 40.22e-6, npart)

particles1 = line.build_particles(
    nemitt_x=2.5e-6, nemitt_y=1e-6,
    # x=np.zeros(npart), y=np.zeros(npart),
    x=x_impact, y=y_impact,
    px=px_impact, py=py_impact,
    zeta=np.zeros(npart), delta=np.zeros(npart), _capacity = int(npart*2))
particles1.pdg_id[:npart] = np.ones(npart)*2212
particles10 = particles1.copy()
bds_link.TrackXSuite(0,'crystaldeflector',particles1,180e3)

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

# npart = 20000
# particles2 = line.build_particles(
#     nemitt_x=2.5e-6, nemitt_y=1e-6,
#     x=np.zeros(npart), y=np.zeros(npart),
#     px=(np.arange(-2,2,(2-(-2))/npart)+0)*1e-6, py=np.zeros(npart),
#     zeta=np.zeros(npart), delta=np.zeros(npart), _capacity = int(npart*2))
# particles2.pdg_id[:npart] = np.ones(npart)*2212
# particles20 = particles2.copy()
# bds_link.TrackXSuite(0,'crystaldeflector2',particles2,180e3)

timestamp = datetime.now().strftime('%Y%m%d_%H')
outfile = f'cry1_{timestamp}.npz'

np.savez(outfile,
    # input
    x_in=particles10.x, px_in=particles10.px,
    y_in=particles10.y, py_in=particles10.py,
    zeta_in=particles10.zeta, delta_in=particles10.delta,
    ptau_in=particles10.ptau,

    # output
    x_out=particles1.x, px_out=particles1.px,
    y_out=particles1.y, py_out=particles1.py,
    zeta_out=particles1.zeta, delta_out=particles1.delta,
    ptau_out=particles1.ptau,

    # identifiers
    particle_id_in=particles10.particle_id, particle_id_out=particles1.particle_id,
    state_in=particles10.state, state_out=particles1.state,   
    at_element_out=particles1.at_element, # should tell you where they were lost
    pdg_id_in=particles10.pdg_id, pdg_id_out=particles1.pdg_id,

    # metadata beam
    bending_angle=50e-6,
    crystal_length=cry1.l,
    crystal_x=cry1.xsize,
    crystal_y=cry1.ysize,
    crystal_material=cry1.material,
    p0c=particles10.p0c[0], 
    mass0=particles10.mass0,

    # metadata simulation
    npart=npart,
    label='cry1',
    cov_avg=cov_avg,
    mu_avg=mu_avg,
)

print(f"Data saved in: {outfile}")

