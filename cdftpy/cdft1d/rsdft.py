#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Main RSDFT module
"""
import io
import logging
from contextlib import redirect_stdout

logging.getLogger('matplotlib').setLevel(logging.WARNING)
import os
import sys
from types import SimpleNamespace

import numpy as np

from cdftpy.cdft1d.coulomb import compute_long_range_coul_pot_kspace, compute_coulomb_potential
from cdftpy.cdft1d.coulomb import compute_long_range_coul_pot_rspace
from cdftpy.cdft1d.coulomb import compute_short_range_coul_pot_rspace
from cdftpy.cdft1d.diis import diis_session
from cdftpy.cdft1d.io_utils import print_banner, print_parameters, print_solute
from cdftpy.cdft1d.io_utils import print_simulation
from cdftpy.cdft1d.potential import compute_lj_potential_mod
from cdftpy.cdft1d.solvent import Solvent, solvent_model_locate, Molecule1
from cdftpy.utils.units import R
from cdftpy import __version__

HEADER = F"""
==================================
1D RSDFT PROGRAM

version {__version__}
Marat Valiev and Gennady Chuev
==================================
"""

kb = R
PI = np.pi

DEFAULT_PARAMS = dict(diis_iterations=2, tol=1.0e-9, output_rate=10, max_iter=200, rcoul=1.25, method="rsdft")

# DG_STRING = "|\u0394\u03B3|"
DG_STRING = "d_g"

DEFAULT_PRINT_LEVEL=frozenset(['header','parameters','solvent','solute'])

def rsdft_1d(solute, solvent, params=None, quiet=False, gr_guess=None,
             print_level = None):

    if print_level is None:
        print_level = DEFAULT_PRINT_LEVEL

    if solvent.nv > 2:
        print("RSDFT method is currently restricted to diatomic liquids")
        print("As an alternative you can try RISM")
        sys.exit(1)
    if quiet:
        sys.stdout = open(os.devnull, 'w')

    qs = solute["charge"]
    sig_s = solute["sigma"]
    eps_s = solute["eps"]

    rho_0 = solvent.density
    sig_v = solvent.sigma
    eps_v = solvent.eps
    qv = solvent.charge

    if params is None:
        params = DEFAULT_PARAMS
    params = {**DEFAULT_PARAMS, **params}

    ndiis = int(params["diis_iterations"])
    rcoul = float(params["rcoul"])
    tol = float(params["tol"])
    max_iter = int(params["max_iter"])
    output_rate = int(params["output_rate"])

    if "temp" not in params:
        params["temp"] = solvent.temp

    if "rmax" in params:
        rmax = params["rmax"]
        solvent.extend(rmax)
    else:
        rmax = solvent.ifft.rgrid[-1]
        params["rmax"] = rmax

    temp = float(params["temp"])
    beta = 1.0 / (kb * temp)

    if 'header' in print_level:
        print(HEADER)
    if 'solute' in print_level:
        print_solute(solute)
    if 'parameters' in print_level:
        print_parameters(params)
    if 'solvent' in print_level:
        solvent.report()

    zeta = solvent.zeta(beta)

    # initialize fft
    ifft = solvent.ifft
    rgrid = ifft.rgrid
    kgrid = ifft.kgrid

    # calculate lj potential
    vlj_r = beta * compute_lj_potential_mod(sig_s, eps_s, sig_v, eps_v, rgrid)

    # coulomb long and short range
    vl_k = beta * compute_long_range_coul_pot_kspace(qs, qv, kgrid, r_s=rcoul)
    vl_r = beta * compute_long_range_coul_pot_rspace(qs, qv, rgrid, r_s=rcoul)
    vcs_r = beta * compute_short_range_coul_pot_rspace(qs, qv, rgrid, r_s=rcoul)

    # total short range potential
    vs_r = vcs_r + vlj_r

    # structure factors
    s_k = solvent.s_k
    sm = solvent.sm_k(kgrid)

    inv_sm = np.linalg.inv(sm.transpose(2, 0, 1))
    delta_s_k = s_k - sm
    hbar = np.einsum('imk,kmj->ijk', np.einsum('kij,jmk->imk', inv_sm, delta_s_k), inv_sm)

    gl_r, gl_k = compute_gamma_long_range(ifft, hbar, sm, vl_k, vl_r, zeta)

    if gr_guess is None:
        g_r = np.zeros(hbar[0].shape)
    else:
        g_r = gr_guess

    print("")
    print_banner("   Self-consistent cycle     ")

    converged = False

    diis_update = diis_session()

    for it in range(max_iter):

        f_r, f_k = compute_mayer_function(ifft, vs_r, gl_r, gl_k, g_r)

        delta_h_r = compute_delta_h(ifft, sm, f_r, f_k)

        hm_r = compute_h_mol(f_r)

        h_r = hm_r + delta_h_r

        cs_k = compute_c_k(ifft, sm, h_r, g_r)

        gn_r = compute_g_r(ifft, hbar, cs_k, gl_r)

        # compute error
        dg_r = gn_r - g_r
        err = np.sum(dg_r ** 2)
        err = np.sqrt(err / gn_r.size)

        converged = err < tol

        if it % output_rate == 0 or converged:
            fe_tot, _ = compute_free_energy(beta, rho_0, ifft, vl_r, g_r, h_r, cs_k)
            if it == 0:
                print(f"{'iter':<5} {DG_STRING:<12}{'Free Energy':<10} ")
            print(f"{it:<5} {err:>.2e}   {fe_tot:<.7f}")

        if converged:
            print(f"\nReached convergence, {DG_STRING} < {tol}")
            g_r = gn_r
            break

        g_r = diis_update(ndiis, gn_r, dg_r)

    if not converged:
        print(
            f"Could not reach specified convergence criteria after {max_iter} iterations"
        )
        # sys.exit(1)

    xi_r = compute_correlation_hole(ifft, sm, f_r, f_k)

    print("\n")
    print(f"{'Total Free Energy ':<30} {fe_tot:>12.6f}")

    cs_r = np.apply_along_axis(ifft.to_rspace, 1, cs_k)

    fe_tot, fe_extra = compute_free_energy(beta, rho_0, ifft, vl_r, g_r, h_r, cs_k)

    fe_gas, fe_inter = fe_extra

    h_k = np.apply_along_axis(ifft.to_kspace, 1, h_r)
    epot_r, epot_k = compute_coulomb_potential(rho_0, qv, ifft, h_k)

    sys.stdout = sys.__stdout__
    return SimpleNamespace(
        method=F"rsdft v{__version__}",
        solute=SimpleNamespace(**solute),
        solvent=solvent,
        ifft=ifft,
        sm=sm,
        f_k=f_k,
        f_r=f_r,
        gl_r=gl_r,
        gl_k=gl_k,
        vl_r=vl_r,
        vl_k=vl_k,
        vs_r=vs_r,
        vcs_r=vcs_r,
        h_r=h_r,
        xi_r=xi_r,
        g_r=g_r,
        fe_tot=fe_tot,
        fe_gas=fe_gas,
        fe_inter=fe_inter,
        phi_r=-(g_r + vl_r) / beta,
        epot_r=epot_r,
        epot_k=epot_k,
        cs_k=cs_k,
        cs_r=cs_r,
        beta=beta,
        rmax=rmax
    )


def compute_mayer_function(ifft, vs_r, gl_r, gl_k, g_r):
    f_r = np.exp(-vs_r + g_r) - 1.0

    delta_f_r = f_r - gl_r

    delta_f_k = np.apply_along_axis(ifft.to_kspace, 1, delta_f_r)

    f_k = delta_f_k + gl_k

    return f_r, f_k


def compute_h_mol(f_r):
    f1 = f_r + 1.0
    hm_r = np.prod(f1, axis=0) - 1

    return hm_r


def compute_delta_h(ifft, sm, f_r, f_k):
    dd = sm - 1.0

    dd_fk = np.einsum("abn,bn->an", dd, f_k)

    dd_fr = np.apply_along_axis(ifft.to_rspace, 1, dd_fk)

    dh_r = dd_fr * (1.0 + f_r)

    return dh_r


def compute_correlation_hole(ifft, sm, f_r, f_k):
    dd = sm - 1.0

    dd_fk = np.einsum("abn,bn->an", dd, f_k)

    dd_fr = np.apply_along_axis(ifft.to_rspace, 1, dd_fk)

    xi_r = np.sum(f_r, axis=0) + dd_fr - f_r

    return xi_r


def compute_c_k(ifft, sm, h_r, g_r):
    h_k = np.apply_along_axis(ifft.to_kspace, 1, h_r)
    g_k = np.apply_along_axis(ifft.to_kspace, 1, g_r)

    c_k = h_k - np.einsum("abn,bn->an", sm, g_k)

    return c_k


def compute_g_r(ifft, hb, c_k, gl_r):
    hb_c_k = np.einsum("abn,bn->an", hb, c_k)

    hb_c_r = np.apply_along_axis(ifft.to_rspace, 1, hb_c_k)

    g_r = gl_r + hb_c_r

    return g_r


def compute_corr_pot(beta, vl_r, g_r):
    return -(g_r + vl_r) / beta


def compute_gamma_long_range(ifft, hb, sm, vl_k, vl_r, zeta):
    """

    Args:
        ifft: FFT instance
        hb: renormalized correlation function
        sm: molecular structure factor
        vl_k: long range coulomb potential in k-space
        vl_r: long range coulomb potential in k-space
        zeta: long range coulomb potential in k-space

    Returns:

    """

    hb_sm = np.einsum("abn,bcn->acn", hb, sm)

    gl_k = -np.einsum("abn,bn->an", hb_sm, vl_k) - vl_k

    delta_gl_k = gl_k + zeta * vl_k

    gl_r = np.apply_along_axis(ifft.to_rspace, 1, delta_gl_k)

    gl_r = gl_r - zeta * vl_r

    return gl_r, gl_k


def compute_free_energy(beta, rho_0, ifft, vl_r, g_r, h_r, c_k) -> float:
    nsites = g_r.shape[0]
    if nsites != 2:
        raise ValueError

    c_r = np.apply_along_axis(ifft.to_rspace, 1, c_k)
    fe0 = -rho_0 * ifft.integrate_rspace(c_r) / nsites / beta

    # alternative way of calculating fe0
    # fe0 = -rho_0 * np.sum(c_k[:, 0]) / nsites / beta

    phi_r = -(g_r + vl_r) / beta

    fe1 = 0.5 * rho_0 * ifft.integrate_rspace(phi_r * h_r)
    fe_tot = fe0 - fe1
    return fe_tot, (fe_tot - fe1, fe1)


if __name__ == "__main__":
    import time

    # load solvent model
    solvent_name = "s2"
    filename = solvent_model_locate(solvent_name)
    solvent = Solvent.from_file(filename)
    # solvent1 = copy.deepcopy(solvent)
    #
    # solvent1.extend(500)

    pass

    solute = dict(name="Cl", charge=-1.0, sigma=4.83, eps=0.05349244)
    m = Molecule1(aname=["Cl"], charge=[-1.0], sigma=[4.83], eps=[0.05349244])
    # params = dict(diis_iterations=2, tol=1.0e-7, max_iter=500)
    # sim = rsdft_1d(solute, solvent, params=params)

    rmax = solvent.ifft.rgrid[-1]
    params = dict(diis_iterations=2, tol=1.0e-7, max_iter=500)

    start_time = time.process_time()
    f = io.StringIO()

    sim = rsdft_1d(solute, solvent, params=params,
                   print_level={'header', 'parameters'})

    # t = time.process_time() - start_time
    # print(rmax, sim.fe_tot,t)
    #
    # analyze_rdf_peaks_sim(sim)

    # for rmax in np.linspace(100, 1900, 10):
    #     params["rmax"] = rmax
    #     start_time = time.process_time()
    #     sim = rsdft_1d(solute, solvent, params=params, quiet=True)
    #     t = time.process_time() - start_time
    #     print(rmax, sim.fe_tot,t)

    # analyze_rdf_peaks_sim(sim)

    # solvent.extend(1000)
    # sim = rsdft_1d(solute, solvent, params=params)
    #
    # analyze_rdf_peaks_sim(sim)
    # write_rdf_sim(sim)