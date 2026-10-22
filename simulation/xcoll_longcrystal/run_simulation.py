import time
from datetime import datetime
import numpy as np
# import matplotlib as mpl
# import matplotlib.pyplot as plt
import simulation_functions as simfun

import xobjects as xo
import xtrack as xt
import xpart as xp
import xcoll as xc


# Make a context and get a buffer
context = xo.ContextCpu()         # For CPU
# context = xo.ContextCupy()      # For CUDA GPUs
# context = xo.ContextPyopencl()  # For OpenCL GPUs

num_part = int(20e6)
particle_ref = xt.Particles('proton', p0c=180e9)

coll = xc.EverestCrystal(length=74e-3, material=xc.materials.Silicon, bending_angle=6e-3,
                         width=0.002, height=0.0128, side='+', miscut=0., lattice='strip', jaw=-0.001,
                         _context=context)

SIGMA_DET = 9e-6
cov_avg, eigvals, n_list = simfun.cov_matrix(sigma_det_x=SIGMA_DET, sigma_det_y=SIGMA_DET)
mu_avg = simfun.mean_vector()  # [mu_x, mu_y, mu_px, mu_py]
rng = np.random.default_rng(seed=82)
samples = rng.multivariate_normal(mean=mu_avg, cov=cov_avg, size=num_part)
x_init  = samples[:, 0]
y_init  = samples[:, 1]
px_init = samples[:, 2]
py_init = samples[:, 3]

# x_init   = np.random.normal(loc=1.5e-3, scale=75.e-6, size=num_part)
# px_init  = np.random.uniform(low=-50.e-6, high=250.e-6, size=num_part)
# y_init   = np.random.normal(loc=0., scale=1e-3, size=num_part)
# py_init  = np.random.normal(loc=0., scale=5.e-6, size=num_part)
part = xp.build_particles(x=x_init, px=px_init, y=y_init, py=py_init,
                          particle_ref=particle_ref, _context=context)
part_init = part.copy()

print(f"Tracking {num_part} particles (Everest)...     ", end='', flush=True)
start = time.time()
coll.track(part)
print(f"Done in {time.time() - start:.2f} seconds.")

# Sort particles to be able to compare to part_init
part.sort(interleave_lost_particles=True)

# # Select only surviving particles and only within the window of interest
# mask = (part.state > 0 ) & ( part.px - part_init.px < 250.e-6) & ( part.px - part_init.px > -50.e-6)

# plt.figure(figsize=(12,8))
# plt.hist2d(part_init.px[mask]*1.e6, part.px[mask]*1.e6 - part_init.px[mask]*1.e6, 500, norm=mpl.colors.LogNorm())
# plt.xlim(-30, 180)
# plt.ylim(-55, 205)
# plt.ylabel(r'$\Delta\theta$ [$\mu$rad]')
# plt.xlabel(r'$\theta_{in}$ [$\mu$rad]')
# plt.tight_layout()
# plt.show()


timestamp = datetime.now().strftime('%Y%m%d_%H')
outfile = f'runs/cry2_{timestamp}.npz'

np.savez(outfile,
    # input
    x_in=part_init.x, px_in=part_init.px,
    y_in=part_init.y, py_in=part_init.py,
    zeta_in=part_init.zeta, delta_in=part_init.delta,
    ptau_in=part_init.ptau,

    # output
    x_out=part.x, px_out=part.px,
    y_out=part.y, py_out=part.py,
    zeta_out=part.zeta, delta_out=part.delta,
    ptau_out=part.ptau,

    # identifiers
    particle_id_in=part_init.particle_id, particle_id_out=part.particle_id,
    state_in=part_init.state, state_out=part.state,   
    at_element_out=part.at_element, # should tell you where they were lost
    pdg_id_in=part_init.pdg_id, pdg_id_out=part.pdg_id,

    # metadata beam
    bending_angle=6e-3,
    crystal_length=coll.length,
    crystal_x=coll.width,
    crystal_y=coll.height,
    crystal_material='Si',
    p0c=part_init.p0c[0], 
    mass0=part_init.mass0,

    # metadata simulation
    npart=num_part,
    label='cry2',
    cov_avg=cov_avg,
    mu_avg=mu_avg,
)

print(f"Data saved in: {outfile}")