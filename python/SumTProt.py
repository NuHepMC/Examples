import pyhepmc
import sys
import functools
import matplotlib.pyplot as plt

#helper function for getting the KE from a HepMC3 GenParticle
def GetKE(part):
  return part.momentum.e - part.momentum.m()

#event projection function
def GetSumTProt(event):
  #look for all physical undecayed final state protons in the event
  all_fs_prot = [p for p in event.particles if ((p.status == 1) and (p.pid == 2212))]

  #sum up their KE
  return functools.reduce(lambda tot, p : tot + GetKE(p), all_fs_prot, 0)

with pyhepmc.open(sys.argv[1]) as file:

  GeVScale = 1
  projections = []
  weights = []

  for i, event in enumerate(file):

    # on the first event verify that we know the units scaling of the input file
    if (i == 0):
      GeVScale = 1 if (event.momentum_unit == pyhepmc.Units.MomentumUnit.GEV) else 1E-3

    sumtprot = GetSumTProt(event)

    if sumtprot > 0:
      projections.append(sumtprot * GeVScale)
      weights.append(event.weights[0])

  plt.hist(projections, bins=100, range=(0,1.5), weights=weights)
  plt.show()