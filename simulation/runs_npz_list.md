# General info from simulation

1. [Available features](#available-features-of-the-particles)
2. [SHORT CRYSTAL RUNS (50x200x2 mm3)](#short-crystal-1)
    - [run 25 aug](#specifics-for-cry1_20260825_16npz)
    - [run 26 aug](#specifics-for-cry1_20260826_15npz)
    - [run 26 aug](#specifics-for-cry1_20260826_16npz-and-cry1_20260826_18npz)
2. [SHORT CRYSTAL RUNS (2x35x4 mm3)](#short-crystal-2)
    - [run 1 sep](#specifics-for-cry1_20260908_13npz)
    - [run 2/3 sep](#specifics-for-cry1_20260914_18npz-and-cry1_20260916_11npz)
    - [run 4/5 sep](#back-to-aper1100-cry1_20260917_16npz-and-cry1_20260918_10npz)
    - [long run - oct](#specifics-for-cry2_20261007_20npz)
3. [LONG CRYSTAL](#long-crystal---lxplus)

$\to$ [run_simulation.py](run_simulation.py)
| parameter | value |
| --- | --- |
| pdg_id (particle type) | `2212` (proton) |


## Available features of the particles

$\to$ full documentation here:  
[https://xsuite.readthedocs.io/en/latest/particlesmanip.html](https://xsuite.readthedocs.io/en/latest/particlesmanip.html)

or to see them directly:
```python
import xtrack as xt

particles_test = xt.Particles(mass0=xt.PROTON_MASS_EV, q0=1, energy0=180e9,
                                x=[0], px=[0], y=[0], py=[0],
                                zeta=[0], delta=[0])

# lista tutti gli attributi pubblici disponibili
print([a for a in dir(particles_test) if not a.startswith('_')])

# oppure, se disponibile, per una panoramica più leggibile 
# (tabella coordinate/reference)
particles1.show()   
```

`['XoStruct', 'add_particles', 'add_to_energy', 'anomalous_magnetic_moment', 'at_element', 'at_turn', 'ax', 'ay', 'beta0', 'charge', 'charge_ratio', 'chi', 'compile_kernels', 'copy', 'delta', 'energy', 'energy0', 'extra_sources', 'filter', 'from_dict', 'from_pandas', 'gamma0', 'get_active_particle_id_range', 'get_classical_particle_radius0', 'get_table', 'hide_first_n_particles', 'hide_lost_particles', 'init_independent_per_part_vars', 'init_pipeline', 'kin_ps', 'kin_px', 'kin_py', 'kin_xp', 'kin_xprime', 'kin_yp', 'kin_yprime', 'kinetic_energy0', 'lost_particles_are_hidden', 'mass', 'mass0', 'mass_ratio', 'merge', 'move', 'p0c', 'parent_particle_id', 'part_energy_varnames', 'part_energy_vars', 'particle_id', 'pdg_id', 'per_particle_vars', 'ptau', 'px', 'py', 'pzeta', 'q0', 'reference_from_pdg_id', 'remove_unused_space', 'reorganize', 'rigidity0', 'rpp', 'rvv', 's', 'scalar_vars', 'set_particle', 'show', 'size_vars', 'sort', 'spin_x', 'spin_y', 'spin_z', 'start_tracking_at_element', 'state', 't_sim', 'to_dict', 'to_json', 'to_pandas', 'to_table', 'unhide_first_n_particles', 'unhide_lost_particles', 'update_beta0', 'update_delta', 'update_gamma0', 'update_p0c', 'update_p0c_and_energy_deviations', 'update_ptau', 'weight', 'x', 'xoinitialize', 'y', 'zeta']`

attribute|what is it?|useful?
--|--|--
`zeta`, `delta` | relative longitudinal position, relative moment deviation | to see correlations energy-channeling
`ptau`, `pzeta`	| normalized energy momentum | redundant with delta for your purpose
`s` |	absolute longitudinal position along the line | not very useful here — you only have a drift + crystal
`at_element` |	index of the element where the particle is located/stops | tells you where a particle was lost with `state != 1`
`at_turn` |	number of turn | irrelevant for me - i don't do multi turns
`weight` |	statistical weight of the macroparticle |  doesn't matter if you don't use weighed macroparticles
`charge`, `mass`, `pdg_id`	| properties of the particle (charge, mass, type)	| static, the same for all (protons) - there's no need to save them for each particle - they're already in the metadata
`spin_x`, `spin_y`, `spin_z` |	 spin components | irrelevant except for polarization studies
`ax`, `ay` |	range of motion (action angle)	| useful only for emission/optical analysis, not for direct channeling

instead data like `beta0`, `gamma0`, `energy0`, `maa0`, `q0`, `p0c`. `kinetic_energy0`, `rigidity0` are properties of the reference particle and so always the same 

## SHORT CRYSTAL 1

| parameter | value |
| --- | --- |
| material | "G4_Si" |
| $x$  | 50 mm |
| $y$  | 200 mm |
| material thickness  | 2 mm |
| $l$ (box around crystal) - keep same as thickness | 2 mm |
| horizontal width (box around crystal) - keep large  | 100 mm |
| bending angle | 50 urad |
| aperture of box around crystal  | rectangular |
| aper1  | 10e1 m |
| x emittance | 2.5e-6 |
| y emittance | 1e-6 |
| zeta coordinate | `zeros(npart)` |
| capacity | `int(npart * 2)` |

### specifics for `cry1_20260826_15.npz`
| parameter | value |
| --- | --- |
| number of particles | 200k|
| x position | `np.zeros(npart)` |
| y position | `np.zeros(npart)` |
| theta_in x divergence | `linspace(-150e-6, 150e-6, npart)`|
| theta_in y divergence | `zeros(npart)` |

![cry1_thetain150_x0.png](plots_cry1/thetain150_x0.png)

### specifics for `cry1_20260826_18.npz`
| parameter | value |
| --- | --- |
| number of particles | 200k |
| x position | `np.random.uniform(-x_half_range, x_half_range, npart)` |
| y position | `np.random.uniform(-y_half_range, y_half_range, npart)` |
| theta_in x divergence | `linspace(-150e-6, 150e-6, npart)`|
| theta_in y divergence | `zeros(npart)` |

![cry1_thetain150_xUnif.png](plots_cry1/thetain150_xUnif.png)

The great peak in $\Delta\theta=0$ is probably given by all of that particles that survived the aperture but did not pass through the crystal (we were not selecting what entered in the crystal and what did not): the `aper` value was 100 meters $\to$ definitely too big!  
Variables `aper1` and `aper2` were changed to 50 mm.  
!!! $\to$ actually 50 mm were giving me problems: too little particles surviving, so it was changed back to 100

## SHORT CRYSTAL 2
**= same dimensions as TCCS =**
| parameter | value |
| --- | --- |
| material | "G4_Si" |
| $x$  | 2 mm |
| $y$  | 35 mm |
| material thickness  | 4 mm |
| $l$ (box around crystal) - keep same as thickness | 4 mm |
| horizontal width (box around crystal) - keep large  | 100 mm |
| bending angle | 50 urad |
| aperture of box around crystal  | rectangular |
| aper1  | 50 mm |
| aper2  | 50 mm |
| x emittance | 2.5e-6 |
| y emittance | 1e-6 |
| zeta coordinate | `zeros(npart)` |
| capacity | `int(npart * 2)` |

### specifics for `cry1_20260908_13.npz`
| parameter | value |
| --- | --- |
| number of particles | 200k |
| x position | `np.random.uniform(-x_half_range, x_half_range, npart)` |
| y position | `np.random.uniform(-y_half_range, y_half_range, npart)` |
| theta_in x divergence | `np.random.uniform(-150e-6, 150e-6, npart)`|
| theta_in y divergence | `np.random.uniform(-150e-6, 150e-6, npart)` |

![uniform distribution](plots_cry1/random_impact_distribution.png)

### specifics for `cry1_20260914_18.npz` and `cry1_20260916_11.npz`

Now instead we are trying to reproduce the gaussian feature of the beam. From real data (8430 run) i fitted the impact positions $(x,y)$ and the entrance angles $(\theta_x,\theta_y)$.

|*Results (8430 run):* ||$\mu$|$\sigma$|
|--|--|--|--|
|all         | $\theta_x$ | -0.93 $\pm$ 0.01 | 26.92 $\pm$ 0.01|
|            | $\theta_y$ | 6.79 $\pm$ 0.01 | 40.22 $\pm$ 0.02|
|all         | $x$ | -0.21 $\pm$ 0.00 | 2.30 $\pm$ 0.00|
|            | $y$ | 0.60 $\pm$ 0.00 | 2.40 $\pm$ 0.00|

| parameter | value |
| --- | --- |
| number of particles | 200k |
|samples | `rng.multivariate_normal(mean=mu_avg, cov=cov_avg, size=npart)`|
| x position | `x_impact = samples[:, 0]` |
| y position | `y_impact = samples[:, 1]` |
| theta_in x divergence | `py_impact = samples[:, 2]`|
| theta_in y divergence | `py_impact = samples[:, 3]` |

### back to `aper1=100`: `cry1_20260917_16.npz` and `cry1_20260918_10.npz`
respectively: 20k and 200k particles  
parameters as before 

### specifics for `cry1_20261008_.npz` and `cry1_20261008_.npz`

#### PROTONS
$10^6$ protons $p$  
```text
Design particle properties: 
Particle:       "proton"
Mass:            0.938272013 GeV
Charge:          1 e
Total Energy:    180.002445412 GeV
Kinetic Energy:  179.064173399 GeV
Momentum:        180 GeV
Gamma:           191.844628123
Beta:            0.999986414561
FFact:           1
Rigidity (Brho): 600.415371357 T*m
```
$\to$ plots in [plots_cry1](./bdsim_shortcrystal/plots_cry1/)  

#### PIONS
$10^6$ pions $\pi^+$  
```text
Design particle properties: 
Particle:       "pi+"
Mass:            0.1395701 GeV
Charge:          1 e
Total Energy:    180.000054111 GeV
Kinetic Energy:  179.860484011 GeV
Momentum:        180 GeV
Gamma:           1289.67489534
Beta:            0.999999699386
FFact:           1
Rigidity (Brho): 600.415371357 T*m
```
changing also [trackerInterface_pi.gmad](./bdsim_shortcrystal/trackerInterface_pi.gmad)
$\to$ plots in [plots_cry1_π](./bdsim_shortcrystal/plots_cry1_π/)



## LONG CRYSTAL - LXPLUS

To count how many files in a folder:
``` bash
ls -1q <folder>/cry2_* | wc -l
```
To see the structure of the files:
``` bash
tree -a -I '.git|node_modules'
```

```text
.
├── condor
│   ├── backup
|   |   ├── logs
│   │   |   └── ...
|   |   ├── output
│   │   |   └── ...
│   │   └── plots_cry2
│   ├── logs
|   |   ├── job_*.log
|   |   ├── job_*.out
│   │   └── job_*.err
│   ├── output
│   │   └── cry2_*.npz
│   ├── htcondor.sub
│   ├── job.py
│   ├── run.sh
│   ├── simulation_functions.py
│   └── trackerInterface.gmad
├── condor2
│   ├── logs
│   │   └── ...
│   ├── output
│   │   └── cry2_*.npz
│   ├── htcondor.sub
│   └── run.sh
├── condor3
│   ├── logs
│   │   └── ...
│   ├── output
│   │   └── cry2_*.npz
│   ├── htcondor.sub
│   └── run.sh
├── plot
│   ├── merge_npz.py
│   ├── plots_cry2
│   │   └── *.png
│   └── plot_simulation.py
├── test
│   ├── cry2_17.npz
│   ├── job.py
│   ├── __pycache__
│   │   └── simulation_functions.cpython-312.pyc
│   ├── simulation_functions.py
│   └── trackerInterface.gmad
├── simulation_functions.py
└── trackerInterface.gmad

17 directories, 26999 files
```

- condor/backup: 2
- condor: 242 $\to$ 2'400'400 particles
- condor2: 1927 $\to$  particles
- condor3: 1600 $\to$  particles
