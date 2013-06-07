import numpy as np
from constants import *

def sigma_psi(psi_P_bar, kappa_a):
    """Compute the dimensionless particle charge density (sigma bar) using the
    dimensionless particle electrical potential (psi_P_bar) and the product of
    the Debye parameter kappa and the particle radius a (kappa_a)."""

    t1 = 2 * np.sinh(psi_P_bar / 2.) - psi_P_bar
    t2 = 4 * np.tanh(psi_P_bar / 4.) - psi_P_bar

    return psi_P_bar + psi_P_bar/kappa_a - t1**2 * kappa_a / (t2 - 
        t1 * kappa_a)

def yukawa_coeff_Magan(a, epsilon, psi_bar, kappa, T):
    """Calculate the Yukawa coefficient with units of (kB * T)
    Arguments:
    a           radius (distance units)
    epsilon     permittivity of solvent (m**-3 kg**-1 s**4 A**2)
    psi_bar     dimensionless electrical potential
    kappa       Debye parameter (1/distance units)
    T           Temperature (K)
    """
    gamma = np.tanh(psi_bar / 4)
    Omega = (psi_bar - 4 * gamma) / (2 * gamma**3)
    return (4 * np.pi * epsilon * a) * (kB * T / e)**2 * ((psi_bar + 4 * gamma * Omega * kappa * a) \
        / (1.0 + Omega * kappa * a))**2

def yukawa_Magan_2006(r, Bpp, kappa, a):
    """Calculate the Yukawa potential with units of (kB * T)
    Arguments:
    r       center-center distance between two particles (units of distance)
    Bpp     Yukawa coefficient (units of kB * T)
    kappa   Debye parameter
    a       particle radius
    Returns: Yukawa potential (units of energy)
    """
    return Bpp * a / r * np.exp(-kappa * (r - 2 * a))

# Dimensionless calculations
def gamma(psi):
    return np.tanh(psi/4.)

def omega(psi, gamma):
    return (psi-4.*gamma)/2./gamma**3

def Bpp(k_T, epsilon, a, psi_P, kappa):
    """Return the dimensionless Yukawa coefficient for particle-particle electrostatic
    interaction potential."""
    gamma = np.tanh(psi_P / 4)
    omega = (psi_P - 4 * gamma) / (2 * gamma**3)

    return (4 * np.pi * k_T * epsilon * a / e0**2)*( (psi_P + 4 * gamma * omega * kappa * a) \
     / (1.0 + omega * kappa * a))**2

def Bps(k_T, epsilon, a, psi_S, psi_P, kappa):
    """Return the dimensionless Yukawa coefficient for the particle-surface electrostatic
    interaction potential."""
    gamma = np.tanh(psi_P / 4)
    omega = (psi_P - 4 * gamma) / (2 * gamma**3)

    return (4 * np.pi * k_T * epsilon * a / e0**2) * (psi_P + 4 * gamma * omega * kappa * a) \
     / (1.0 + omega * kappa * a) * 4 * np.tanh(psi_S / 4.)

def Upp_el(Bpp, kappa, a, r):
    """Return the dimensionless electrical potential between two particles
    Bpp: Yukawa coefficient
    kappa: Debye parameter of medium
    a: particle radius
    r: dimensionless distance between particle centers"""

    return Bpp / r * np.exp(-kappa * a * ( r - 2.0))

def Ups_el(Bps, kappa, a, h):
    """Return the dimensionless electrical potential between a particle and a surface.
    Bps: Yukawa coefficient
    kappa: Debye parameter of medium    
    a: particle radius
    h: dimensionless distance from surface of particle to surface"""
    return Bps * np.exp(-kappa * a * h)
