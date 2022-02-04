#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Command line interface to 1D CDFT calculations
"""

import os
from scipy.signal import argrelmax
from scipy.signal import argrelmin
from cdftpy.cdft1d.io_utils import print_banner
from cdftpy import IonSolvation

def write_rdf_sim(sim, dirpath="./"):
    write_rdf(sim.ifft, sim.name, sim.solvent.aname, sim.h_r, dirpath=dirpath)

def write_rdf(ifft, solute_name, solvent_name, h_r, dirpath="./"):
    rgrid = ifft.rgrid
    print("\nGenerating solvent-solute RDF's")
    for i in range(len(solvent_name)):
        filename = os.path.join(dirpath, f"rdf_{solute_name}{solvent_name[i]}.dat")
        print(f"Generating {os.path.abspath(filename)}")
        with open(filename, "w") as fp:
            fp.write(f"{'# r':<10} {'g(r)':<10}\n")

            for j in range(len(rgrid)):
                r = rgrid[j]
                rdf = h_r[i, j] - h_r[i, 0]
                fp.write(f"{r:<10.4} {rdf:<10.4}\n")

# Au using UFF potential.
sim = IonSolvation(charge=0.0,sigma=2.934,eps=0.1627576, solvent="s2", max_iter=2000)
free_energy = sim.cdft()

print("free_energy, kj/mol", free_energy)

write_rdf_sim(sim)
