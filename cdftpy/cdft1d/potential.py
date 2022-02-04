#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
This module provides collection
of routines for Lennard-Jones potential calculations
"""
import numpy as np


def compute_lj_potential(sig_solute, eps_solute, sig_solvent, eps_solvent, rgrid):

    eps_solute = float(eps_solute)
    sig_solute = float(sig_solute)

    eps = eps_solute * eps_solvent
    eps = np.sqrt(eps)

    sig = 0.5 * (sig_solute + sig_solvent)

    sig6 = sig ** 6

    r_m6 = 1.0 / rgrid ** 6
    r_m12 = 1.0 / rgrid ** 12
    sig12 = sig ** 12

    lj_repuls = 4.0 * np.multiply.outer(eps * sig12, r_m12)
    lj_attrac = -4.0 * np.multiply.outer(eps * sig6, r_m6)

    lj = lj_repuls + lj_attrac
    return lj


def compute_lj_potential_mod(sig_solute, eps_solute, sig_solvent, eps_solvent, rgrid):

    eps_solute = float(eps_solute)
    sig_solute = float(sig_solute)

    eps = eps_solute * eps_solvent
    eps = np.sqrt(eps)

    sig = 0.5 * (sig_solute + sig_solvent)

    sig6 = sig ** 6
    sig6[1] = 0

    r_m6 = 1.0 / rgrid ** 6
    r_m12 = 1.0 / rgrid ** 12
    sig12 = sig ** 12

    lj_repuls = 4.0 * np.multiply.outer(eps * sig12, r_m12)
    lj_attrac = -4.0 * np.multiply.outer(eps * sig6, r_m6)

    lj = lj_repuls + lj_attrac
    return lj


def compute_lj_nano_potential(sig_solute, eps_solute, sig_solvent, eps_solvent, rgrid):

    eps_solute = float(eps_solute)
    sig_solute = float(sig_solute)

    eps = eps_solute * eps_solvent
    eps = np.sqrt(eps)

    sig = 0.5 * (sig_solute + sig_solvent)
    sig[1] = 0

    nano_radius = 20 * sig_solute

    #rho_nano = 0.0589
    #v_nano = (4/3)* np.pi * nano_radius**3
    #number_beads = v_nano * rho_nano

    number_beads = 49851
    v_nano = (4/3)* np.pi * nano_radius**3
    rho_nano = number_beads/v_nano

    A_ns = 24 * np.pi * eps * rho_nano * sig ** 3

    print("sigma_ns (beads-solvent) = ", sig, " Ang")
    print("eps_ns (beads-solvent) = ", eps, " kj/mol")
    print("nanoparticle radius = ", nano_radius, " Ang")
    print("nanoparticle diameter = ", 2 * nano_radius, " Ang")
    print("number of beads in nanoparticle", number_beads )
    print("Hamaker constant (A_ns) = ", A_ns, " kj/mol")

    term_top = np.multiply.outer(sig ** 6,  ( (5 * nano_radius ** 6) + (45 * (nano_radius ** 4) * (rgrid ** 2)) + (
                63 * (nano_radius ** 2) * (rgrid ** 4)) + (15 * (rgrid ** 6))) )
    term_bottom = 15 * (nano_radius ** 2 - rgrid ** 2) ** 6
    term_ratio = 1 - (term_top / term_bottom)

    term_factor = (2 / 9)  * ((nano_radius ** 3) * (sig ** 3) * A_ns )
    term_bottom0 = 1 / ((nano_radius ** 2 - rgrid ** 2) ** 3 )
    lj_nano = np.multiply.outer(term_factor, term_bottom0) * term_ratio

    print("info about potential at grid[5900]", lj_nano[0,5900], lj_nano[1,5900] )
    for i in range(5900):
        lj_nano[:,i] = lj_nano[:,5900]
    print("info about potential at grid[590]",  lj_nano[0,590], lj_nano[1,590] )


    #print("term_top",term_top)
    #print("term_bottom",term_bottom)
    #print("term_ratio",term_ratio)
    #print("term_factor",term_factor)
    #print("term_bottom0",term_bottom0)
    #print("rgrid",rgrid)
    #print("lj_nano",lj_nano)

    dataout = np.column_stack((rgrid, lj_nano[0]))
    np.savetxt('lj_nano.dat', dataout)

    return lj_nano