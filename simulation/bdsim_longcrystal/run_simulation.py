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

# long crystal 74 mm
cry2 = bdsim.Element()
cry2.type = bdsim.elementtype.ElementType.USERCOMPONENT
cry2.userTypeName = "crystaldeflector"
cry2.name = "crystaldeflector2"
cry2.material="G4_Si"
cry2.xsize = 5e-3 # m
cry2.ysize = 12.8e-3 # m
cry2.materialThickness = 74e-3 # m
cry2.l = 74e-3
cry2.horizontalWidth = 0.1
cry2.userParameters="crystalRegion:crystaldeflector crystalLattice:(110) crystalBendingAngle:6e-3;"
cry2.apertureType = "rectangular" # aperture of box around crystal
cry2.aper1 = 10e1 # m
cry2.aper2 = 10e1 # m
l.AddLinkElement(cry2)

#correlated gaussian beam
cov_avg, eigvals, n_list = simfun.cov_matrix()
assert eigvals.min() > 0, f"cov_avg non è definita positiva! min eigenvalue = {eigvals.min()}"

mu_avg = simfun.mean_vector()  # [mu_x, mu_y, mu_px, mu_py]
rng = np.random.default_rng(seed=82)  # same seed of BDSLinkTrackerInterface, per coerenza/riproducibilità

npart = 20000
samples = rng.multivariate_normal(mean=mu_avg, cov=cov_avg, size=npart)
x_impact  = samples[:, 0]
y_impact  = samples[:, 1]
px_impact = samples[:, 2]
py_impact = samples[:, 3]

particles1 = line.build_particles(
    nemitt_x=2.5e-6, nemitt_y=1e-6,
    x=x_impact, y=y_impact,
    px=px_impact, py=py_impact,
    zeta=np.zeros(npart), delta=np.zeros(npart), _capacity = int(npart*2))
particles1.pdg_id[:npart] = np.ones(npart)*2212
particles10 = particles1.copy()
bds_link.TrackXSuite(0,'crystaldeflector',particles1,180e3)

timestamp = datetime.now().strftime('%Y%m%d_%H')
outfile = f'cry2_{timestamp}.npz'

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
    # bending_angle=50e-6,
    # crystal_length=cry1.l,
    # crystal_x=cry1.xsize,
    # crystal_y=cry1.ysize,
    # crystal_material=cry1.material,
    bending_angle=6e-3,
    crystal_length=cry2.l,
    crystal_x=cry2.xsize,
    crystal_y=cry2.ysize,
    crystal_material=cry2.material,
    p0c=particles10.p0c[0], 
    mass0=particles10.mass0,

    # metadata simulation
    npart=npart,
    label='cry1',
    cov_avg=cov_avg,
    mu_avg=mu_avg,
)

print(f"Data saved in: {outfile}")

