#
# pip install xsuite, matplotlib, xplt, pint 
#
import bdsim
import xpart as xp
import matplotlib.pyplot as plt
import xtrack as xt
import numpy as np
env = xt.Environment()
env.set_particle_ref('proton', kinetic_energy0=200e6)
env.vars.default_to_zero = True # Undefined variables in env are automatically 
                                # set to zero.

plot = False

# Element geometry
n_bends = 16
env['ang_mb'] = 2*np.pi/n_bends
env['l_mb'] = 1.65
env['l_mq'] = 0.35

env.new('mb', xt.RBend, length_straight='l_mb', angle='ang_mb')
env.new('mq', xt.Quadrupole, length='l_mq');

# Quadrupole families with different strengths
env.new('qfa', 'mq', k1= 'kqfa')
env.new('qfb', 'mq', k1= 'kqfb')
env.new('qd',  'mq', k1= 'kqd');

cell_a = env.new_line(length=7.405, components=[
    env.place('qfa', at=0.3875),
    env.place('mb', at=1.8125),
    env.place('qd', at=3.2925),
    env.place('mb', at=5.0475),
    env.place('qfa', at=6.3275),
])

if plot:
    cell_a.survey().plot()
    plt.show()

cell_b = env.new_line(name='cell_b', length=8.405, components=[
    env.place('qfb', at=1.2725),
    env.place('mb', at= 2.7275),
    env.place('qd', at=4.8575),
    env.place('mb', at=6.5125),
    env.place('qfb', at=7.7925),
])

if plot:
    plt.figure(2)
    cell_b.survey().plot()
    plt.show()

# concatenate the two cells
arc = cell_a + cell_b
if plot: 
    arc.survey().plot()
    plt.show()

long_straight = env.new_line(length=2., components=[
    env.new('bgv.s', xt.Marker, at=0.),
    env.new('bgv.e', xt.Marker, at=0.1),
    env.new('mid.lss', xt.Marker, at=1.),
    env.new('coll.s', xt.Marker, at=1.1),
    env.new('coll.e', xt.Marker, at=1.5),
])
short_straight = env.new_line(length=1., components=[
    env.new('mid.sss', xt.Marker, at=1.)
])

half_ring = (long_straight
             + arc
             + short_straight
             - arc # mirror symmetric lattice
            )
if plot:
    half_ring.survey().plot()
    plt.show()

ring = 2 * half_ring
if plot:
    ring.survey().plot()
    plt.show()
	
ring.replace_all_repeated_elements() # give all elements unique names
tt = ring.get_table()
tt.cols['element_type s_start s_center s_end']

# Inspect all quadrupoles
tt_quad = tt.rows[tt.element_type=='Quadrupole']
tt_quad.cols['s_start s_center s_end']

# Tag all quadrupoles in survey plot
sv = ring.survey()
if plot:
    sv.plot(labels=tt_quad.name);
    plt.show()

# Magnet type
env.new('ms', xt.Sextupole, length=0.2)

# Magnet instances
env.new('msf.1', 'ms', k2='ksf')
env.new('msf.2', 'ms', k2='ksf')
env.new('msd.1', 'ms', k2='ksd')
env.new('msd.2', 'ms', k2='ksd')
env.new('mse',   'ms', k2='kse');

ring.insert([
    env.place('msf.1', at=-0.2, from_='qfb.0@start'),
    env.place('msf.2', at=-0.2, from_='qfb.4@start'),
    env.place('msd.1', at=0.3,  from_='qd.2@end'),
    env.place('msd.2', at=0.3,  from_='qd.6@end'),
    env.place('mse',   at=-0.3, from_='qfa.4@start')
])
	
# Inspect sextupoles in the survey
sv = ring.survey()
if plot:
    sv.plot(labels=['msf.1', 'msf.2', 'msd.1', 'msd.2', 'mse'])
    plt.show()

env.new('rf1', xt.Cavity, voltage='vrf', frequency='frf')
ring.insert('rf1', at=0.5, from_='qfa.3@start')

env['kqfa'] = 0.33772844308134037
env['kqfb'] = 0.5469407043373885
env['kqd']  = -0.5904728358276611
env['ksf'] = -2.8457172165156828
env['ksd'] = 2.5169216191564843
env['vrf'] = 1e6            # V
tw4d = ring.twiss(method='4d') 
env['frf'] = 20 / tw4d.t_rev0 # h=1

tw = ring.twiss4d()
if plot:
    pl = tw.plot()
    pl.ylim(left_hi=40, right_lo=-20, right_hi=20,
            lattice_hi=1.5, lattice_lo=-7)
    plt.show()

### define the BDSIM elements:

bds_link = bdsim.BDSLinkTrackerInterface.GetInstance("./trackerInterface.gmad", \
                                                    referenceParticlePDG=2212, \
                                                    referenceKineticEnergy=200e0, \
                                                    relativeEnergyCut=0.1, \
                                                    batchMode=True)
l = bds_link.GetBDSIMLink()

coll1 = bdsim.Element()
coll1.type = bdsim.elementtype.ElementType.TARGET
coll1.set_value('material','He')
coll1.name = 'coll1'
coll1.set_value('l',0.001)
coll1.set_value('materialThickness',0.001)
coll1.set_value('ysize',2.0)
coll1.set_value('horizontalWidth',2.0)
l.AddLinkElement(coll1)     

coll2 = bdsim.Element()
coll2.type = bdsim.elementtype.ElementType.JCOL
coll2.set_value('material','W')
coll2.name = 'coll2'
coll2.set_value('l',0.4)
coll2.set_value('xsizeLeft',10e-3)
coll2.set_value('xsizeRight',10e-3)
coll2.set_value('ysize',0.2)
coll2.set_value('horizontalWidth',2.0)
l.AddLinkElement(coll2)    

def calc_emit(part):
    m = part.state == 1
    n = np.sum(m)
    xs = part.x[m]
    pxs = part.px[m]
    ys = part.y[m]
    pys = part.py[m]
    xsq = np.sum((xs-np.mean(xs))**2)/n
    pxsq = np.sum((pxs-np.mean(pxs))**2)/n
    xpx = np.sum((xs-np.mean(xs))*(pxs-np.mean(pxs)))/n
    ex = np.sqrt(xsq*pxsq - xpx**2)
    ys = part.y[m]
    pys = part.py[m]
    ysq = np.sum((ys-np.mean(ys))**2)/n
    pysq = np.sum((pys-np.mean(pys))**2)/n
    ypy = np.sum((ys-np.mean(ys))*(pys-np.mean(pys)))/n
    ey = np.sqrt(ysq*pysq - ypy**2)  
    return ex,ey



bunch = xp.generate_matched_gaussian_bunch(line=ring,num_particles=10000,
    total_intensity_particles=1e11, nemitt_x=2e-6, nemitt_y=2e-6, sigma_z=10)

part = xp.Particles(
    x=bunch.x, px=bunch.px,
    y=bunch.y, py=bunch.py,
    zeta=bunch.zeta, ptau=bunch.ptau,
    p0c=bunch.p0c, pdg_id=bunch.pdg_id, _capacity=20000)  

ring.track(part,num_turns=1,ele_stop='bgv.s.1')
m = part.state == 1
px0 = np.max(np.abs(part.px[m]))
py0 = np.max(np.abs(part.py[m]))

#for turn in np.arange(100):
#    bds_link.TrackXSuite(0,'coll1',part,env.particle_ref.p0c[0]*1e-6)
#    ex,ey = calc_emit(part)
#    print(f'ex: {ex}, ey: {ey}')
#    m = part.state == 1
#    if m.sum() == 0:
#        break    
                                    
coords_at_coll1 = []
coords_at_coll2 = []
emittance_x = []
emittance_y = []
n_parts = []
for turn in np.arange(100):
    print(turn)
    bds_link.TrackXSuite(0,'coll1',part,env.particle_ref.p0c[0]*1e-6)
    part.beta0 = np.ones(len(part.beta0))*part.beta0[0]
    part.update_delta(part.delta)
    m = part.state == 1
    print(f'n surviving particles: {m.sum()}')
    if m.sum() == 0:
        break
    ex,ey = calc_emit(part)
    emittance_x.append(ex)
    emittance_y.append(ey)
    coords_at_coll1.append([part.x[m],part.y[m]])           
    ring.track(part,num_turns=1,ele_start='bgv.s.1',ele_stop='coll.s.1')
    bds_link.TrackXSuite(1,"coll2",part,env.particle_ref.p0c[0]*1e-6)
    part.beta0 = np.ones(len(part.beta0))*part.beta0[0]
    part.update_delta(part.delta) 
    ring.track(part,num_turns=1,ele_start='coll.s.1',ele_stop='bgv.s.1')
    m = part.state > -10000
    print(f'px0: {px0}, py0: {py0}')
    m = part.state == 1
    px1 = np.max(np.abs(part.px[m]))
    py1 = np.max(np.abs(part.py[m]))
    print(f'px1: {px1}, py1: {py1}')
    n_parts.append(m.sum())
    print(f'n surviving particles: {m.sum()}')
    if m.sum() == 0:
        break
    coords_at_coll2.append([part.x[m],part.y[m]])           
emittance_x = np.array(emittance_x)
emittance_y = np.array(emittance_y)
n_parts = np.array(n_parts)
coords_at_coll1 = np.array(coords_at_coll1[-1])
coords_at_coll2 = np.array(coords_at_coll2[-1])

FONT  = 14
TITLE = 16
LW    = 2.0
turns = np.arange(len(emittance_x))
 
fig, ax1 = plt.subplots(figsize=(10, 5))
 
# ── Left axis: emittances ─────────────────────────────────────────────────────
color_x = "#1f77b4"   # blue
color_y = "#ff7f0e"   # orange
 
l1, = ax1.plot(turns, emittance_x * 1e6, color=color_x, lw=LW, label=r"$\varepsilon_x$")
l2, = ax1.plot(turns, emittance_y * 1e6, color=color_y, lw=LW, label=r"$\varepsilon_y$", ls="--")
 
ax1.set_xlabel("Turn", fontsize=FONT)
ax1.set_ylabel(r"Emittance  [$\mu$m·rad]", fontsize=FONT, color="black")
ax1.tick_params(axis="both", labelsize=FONT - 1)
ax1.tick_params(axis="y", colors="black")
 
# ── Right axis: number of particles ──────────────────────────────────────────
color_n = "#2ca02c"   # green
 
ax2 = ax1.twinx()
l3, = ax2.plot(turns, n_parts, color=color_n, lw=LW, label=r"$N_\mathrm{particles}$", ls=":")
 
ax2.set_ylabel("Number of particles", fontsize=FONT, color=color_n)
ax2.tick_params(axis="y", labelsize=FONT - 1, colors=color_n)
ax2.spines["right"].set_edgecolor(color_n)
 
# ── Shared legend ─────────────────────────────────────────────────────────────
lines  = [l1, l2, l3]
labels = [l.get_label() for l in lines]
ax1.legend(lines, labels, fontsize=FONT - 1, loc="upper right", framealpha=0.9)
 
# ── Cosmetics ─────────────────────────────────────────────────────────────────
ax1.set_title("Emittance & Beam Intensity vs. Turn", fontsize=TITLE, pad=10)
ax1.grid(True, linestyle="--", alpha=0.4)
fig.tight_layout()
 
plt.show()


np.unique(part.at_element[part.state==-333], return_counts=True)

data = part.x[part.at_element == 68]
data2 = data[np.abs(data)<20e-3]

fig, ax = plt.subplots(figsize=(8, 5))

ax.hist(data2, bins=50, color="#1f77b4", edgecolor="white", linewidth=0.5)

ax.set_xlabel("x  [m]", fontsize=14)
ax.set_ylabel("Counts", fontsize=14)
ax.set_title("Horizontal position distribution at element 68", fontsize=15, pad=10)
ax.tick_params(axis="both", labelsize=13)
ax.grid(True, linestyle="--", alpha=0.4)

fig.tight_layout()
plt.show()

