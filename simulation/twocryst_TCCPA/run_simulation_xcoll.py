import time
from datetime import datetime
import numpy as np
# import matplotlib as mpl
# import matplotlib.pyplot as plt
# import simulation_functions as simfun

import xobjects as xo
import xtrack as xt
import xpart as xp
import xcoll as xc


# Make a context and get a buffer
context = xo.ContextCpu()         # For CPU
# context = xo.ContextCupy()      # For CUDA GPUs
# context = xo.ContextPyopencl()  # For OpenCL GPUs

num_part = int(1e7)
particle_ref = xt.Particles('proton', p0c=180e9)

coll = xc.EverestCrystal(length=70.5e-3, material=xc.materials.Silicon, bending_angle=13.3e-3,
                         width=0.002, height=0.0225, side='+', miscut=0., lattice='strip', jaw=-0.001,
                         _context=context)

sigma_x = 2e-3
sigma_y = 2e-3
sigma_px = 60e-6
sigma_py = 60e-6
rng = np.random.default_rng(seed=17)
x_impact  = rng.normal(0.0, sigma_x, num_part)
y_impact  = rng.normal(0.0, sigma_y, num_part)
px_impact = rng.normal(0.0, sigma_px, num_part)
py_impact = rng.normal(0.0, sigma_py, num_part)

part = xp.build_particles(x=x_impact, px=px_impact, y=y_impact, py=py_impact,
                          particle_ref=particle_ref, _context=context)
part_init = part.copy()

print(f"Tracking {num_part} particles (Everest)...     ", end='', flush=True)
start = time.time()
coll.track(part)
print(f"Done in {time.time() - start:.2f} seconds.")

# Sort particles to be able to compare to part_init
part.sort(interleave_lost_particles=True)

outfile = f'long_xcoll.npz'

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
    bending_angle=13.3e-3,
    crystal_length=coll.length,
    crystal_x=coll.width,
    crystal_y=coll.height,
    crystal_material='Si',
    p0c=part_init.p0c[0], 
    mass0=part_init.mass0,

    # metadata simulation
    npart=num_part,
    label='cry2',
    # cov_avg=cov_avg,
    # mu_avg=mu_avg,
    sigma_x = sigma_x,
    sigma_y = sigma_y,
    sigma_px = sigma_px,
    sigma_py = sigma_py,
)

print(f"Data saved in: {outfile}")