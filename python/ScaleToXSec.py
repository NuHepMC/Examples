import pyhepmc
import typing
import sys
import matplotlib.pyplot as plt
import numpy as np

#event projection function
def GetQ2(event):
  #get the hard scatter vertex
  prim_vertex = [v for v in event.vertices if (v.status == 1)][0]

  in_beam = [p for p in prim_vertex.particles_in if (p.status == 4)][0]

  #determine expeced outgoing lepton pid
  cc_lep_pid = in_beam.pid + (1 if in_beam.pid < 0 else -1)

  #look for an outgoing particle from the hard scatter with the right pid
  #might subsequently undergo FSI, so may not be status == 1, but we want the momentum transfer at the vertex
  out_charged_leps = [p for p in prim_vertex.particles_out if (p.pid == cc_lep_pid)]

  #if no expected outgoing leptons then it isn't a valid CC event
  if len(out_charged_leps) == 0:
    return None

  return -(in_beam.momentum - out_charged_leps[0].momentum).interval()

with pyhepmc.open(sys.argv[1]) as file:

  GeV2Scale = 1
  projections = []
  weights = []

  have_GC5 = False
  have_EC2 = False
  have_EC4 = False

  GC5_weight = 1
  EC4_best_estimate = 1

  for i, event in enumerate(file):

    if (i == 0):
      # on the first event verify that we know the units scaling of the input file
      GeV2Scale = 1 if (event.momentum_unit == pyhepmc.Units.MomentumUnit.GEV) else 1E-6

      conventions = []
      if "NuHepMC.Conventions" in event.run_info.attributes:
        conventions = event.run_info.attributes["NuHepMC.Conventions"].astype(typing.List[str])

      if "G.C.5" in conventions and "G.C.2" in conventions:
        have_GC5 = True
        print("Have information to scale to XSec via G.C.5+G.C.2")
        FATX = event.run_info.attributes["NuHepMC.FluxAveragedTotalCrossSection"].astype(float)
        NEvents = event.run_info.attributes["NuHepMC.Exposure.NEvents"].astype(float)
        GC5_weight = FATX / NEvents

      if "E.C.2" in conventions and "G.C.7" in conventions:
        have_EC2 = True
        print("Have information to scale to XSec via E.C.2+G.C.7")

      if "E.C.4" in conventions:
        have_EC4 = True
        print("Have information to scale to XSec via E.C.4")

      if  (not have_GC5) and \
          (not have_EC2) and \
          (not have_EC4):
          print("We can not find the required information to correctly normalize the input vector %s" % sys.argv[1])

    Q2 = GetQ2(event)

    if have_EC4:
      EC4_best_estimate = event.cross_section.astype(pyhepmc.GenCrossSection).xsec()
    
    if Q2 is not None:
      projections.append(Q2 * GeV2Scale)
      weights.append(event.weights[0])

  #use a numpy array for ease of multiplying
  weights_arr = np.array(weights)
  if have_GC5:
    weights_arr = weights_arr * GC5_weight
    print("G.C.5 Scale: %s" % GC5_weight)
  elif have_EC4:
    weights_arr = weights_arr * EC4_best_estimate
    print("E.C.4 Scale: %s" % EC4_best_estimate)
  elif have_EC2:
    pass

  #plot the weighted histogram dividing by the bin width
  plt.hist(projections, bins=100, range=(0,3), weights=(weights_arr * (100.0/3.0)))
  plt.show()
