import pyhepmc
import sys
import matplotlib.pyplot as plt

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

  for i, event in enumerate(file):

    # on the first event verify that we know the units scaling of the input file
    if (i == 0):
      GeV2Scale = 1 if (event.momentum_unit == pyhepmc.Units.MomentumUnit.GEV) else 1E-6

    Q2 = GetQ2(event)
    
    if Q2 is not None:
      projections.append(Q2 * GeV2Scale)
      weights.append(event.weights[0])

  plt.hist(projections, bins=100, range=(0,3), weights=weights)
  plt.show()
