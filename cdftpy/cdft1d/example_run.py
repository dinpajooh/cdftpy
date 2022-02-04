#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Command line interface to 1D CDFT calculations
"""

from cdftpy import IonSolvation
from cdftpy.cdft1d.rdf import write_rdf_sim

# Au using UFF potential.
sim = IonSolvation(charge=0.0,sigma=2.934,eps=0.1627576, solvent="s2", max_iter=2000)
free_energy = sim.cdft()

print("free_energy, kj/mol", free_energy)

write_rdf_sim(sim)
