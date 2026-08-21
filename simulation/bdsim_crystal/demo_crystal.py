import numpy as np
import xtrack as xt
import bdsim
import xpart as xp

drift = xt.Drift(length=1)

line = xt.Line(elements=[drift])
line.particle_ref = xt.Particles(mass0=xt.PROTON_MASS_EV, q0=1, energy0=180e9)

bds_link = bdsim.BDSLinkTrackerInterface.GetInstance("../bdsim_crystal/trackerInterface.gmad",
                                                     referenceParticlePDG=2212,
                                                     referenceKineticEnergy=180e3, # MeV!!!! !!!
                                                    relativeEnergyCut=0.1,
                                                    batchMode = True,
                                                    seed = 82)
l = bds_link.GetBDSIMLink()


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

cry2 = bdsim.Element()
cry2.type = bdsim.elementtype.ElementType.USERCOMPONENT
cry2.userTypeName = "crystaldeflector"
cry2.name = "crystaldeflector2"
cry2.material="G4_Si"
cry2.xsize = 0.05
cry2.ysize = 0.2
cry2.materialThickness = 74e-3
cry2.l = 74e-3
cry2.horizontalWidth = 0.1
cry2.userParameters="crystalRegion:crystaldeflector crystalLattice:(110) crystalBendingAngle:6e-3;"
cry2.apertureType = "rectangular" # aperture of box around crystal
cry2.aper1 = 10e1 # m
l.AddLinkElement(cry2)

npart = 200
particles1 = line.build_particles(
    nemitt_x=2.5e-6, nemitt_y=1e-6,
    x=np.zeros(npart), y=np.zeros(npart),
    px=(np.arange(-2,2,(2-(-2))/npart)+0)*1e-6, py=np.zeros(npart),
    zeta=np.zeros(npart), delta=np.zeros(npart), _capacity = int(npart*2))
particles1.pdg_id[:npart] = np.ones(npart)*2212
particles10 = particles1.copy()
bds_link.TrackXSuite(0,'crystaldeflector',particles1,180e3)

npart = 200
particles2 = line.build_particles(
    nemitt_x=2.5e-6, nemitt_y=1e-6,
    x=np.zeros(npart), y=np.zeros(npart),
    px=(np.arange(-2,2,(2-(-2))/npart)+0)*1e-6, py=np.zeros(npart),
    zeta=np.zeros(npart), delta=np.zeros(npart), _capacity = int(npart*2))
particles2.pdg_id[:npart] = np.ones(npart)*2212
particles20 = particles2.copy()
bds_link.TrackXSuite(0,'crystaldeflector2',particles2,180e3)