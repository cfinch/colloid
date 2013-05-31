#!/usr/bin/env python

from numpy import *

# ----- Available Surface Functions and Blocking Functions -----
def phi_theta_Langmuir(theta, theta_max=1.0):
    return 1.0 - theta / theta_max

def phi_theta_Schaaf(theta):
    """Approximates the blocking function phi(theta) as shown in Eqn.31 of
    Schaaf and Talbot, 1989."""
    alpha3 = 40.0 / sqrt(3.0) / pi - 176.0 / (3 * pi**2)
    return 1.0 - 4*theta + 6*sqrt(3.0) / pi * theta**2 + alpha3 * theta**3

def phi_theta_Rubi(gamma, theta):
    """Computes the available surface function phi(gamma, theta) according
    to the approximation shown in Wojtaszczyk, Avalos, and Rubi, Europhys. 
    Lett., 40 (3), pp. 299-304, eqn. 9"""
    x = theta / 0.547
    return (1.0 - x * (1.0 - gamma**2))**3 / (1.0 - 0.812 * x * (1-gamma**2) \
            + 0.237 * x**2 * (1.0 - gamma**2)**2 + 0.084 * x**3 \
            * (1.0 - gamma**2)**3 + 0.233 * x**4 * (1.0 - gamma**2)**4)

def phi_theta_fit2(theta):
    """Approximates the blocking function phi(theta) as shown in Eqn.41 of
    Schaaf and Talbot, 1989."""
    theta_bar = theta / 0.547
    return (1.0 + 0.812 * theta_bar + 0.4258 * theta_bar**2 + 0.0716 \
            * theta_bar**3) * (1.0 - theta_bar)**3

def phi_theta_fit3(theta):
    """Approximates the blocking function phi(theta) as shown in Eqn.42 of
    Schaaf and Talbot, 1989."""
    theta_bar = theta / 0.547
    return (1.0 - theta_bar)**3 / (1.0 - 0.812 * theta_bar + 0.2336 \
            * theta_bar**2 + 0.0845 * theta_bar**3)

def phi_theta_Brownian(theta):
    """Fit to my own Brownian Dynamics simulations.  Values came from fitting
    volume fraction 0.05, D=3.77e-12, box size 120a"""
    theta_bar = theta / 0.547
    return (1.0 - theta_bar)**3 / (1.0 - 0.7906 * theta_bar + -0.0017 \
            * theta_bar**2 - 0.1943 * theta_bar**3)

def ASF_spherocylinders(theta, alpha):
    """From Viot, Tarjus, Ricci, and Talbot 1992"""
    theta_max = 0.554    # alpha = 4  Need to make this a function of alpha.
    theta_bar = theta / theta_max
    d1 = -1.240
    d2 = 0.736
    return (1.0 - theta_bar)**4 / (1.0 + d1 * theta_bar + d2 * theta_bar**2)

# Van Tassel's two-stage model
def Phi_alpha(theta_a, theta_b, ratio):
    """Blocking function for adsorption of particles in the model by Brusatori
    and Van Tassel.
    Arguments:
        theta_a   dimensionless number density of molecules in state a 
        theta_b   dimensionless number density of molecules in state b
        ratio   (area of particle in state b)/(area of particle in state a)
    """
    R = sqrt(ratio)
    rho_a = theta_a
    rho_b = theta_b/ratio
    theta = rho_a + R**2 * rho_b

    return (1.0 - theta) * (exp(-2 * (rho_a + R * rho_b) 
        / (1.0-theta) - (rho_a + rho_b + (R - 1)**2 * rho_a
        * rho_b) / (1.0 - theta)**2))

def Psi_alpha_beta(theta_a, theta_b, ratio):
    """Blocking function for adsorption of particles in the model by Brusatori
    and Van Tassel.
    Arguments:
        theta_a   dimensionless number density of molecules in state a
        theta_b   dimensionless number density of molecules in state b
        ratio   (area of particle in state b)/(area of particle in state a)
    """
    R = sqrt(ratio)
    rho_a = theta_a
    rho_b = theta_b/ratio
    theta = rho_a + R**2 * rho_b

    return exp(-2 * (R - 1.0) * (rho_a + R * rho_b) / (1.0
        -theta) - (R**2 - 1.0) * (rho_a + rho_b + 
            (R - 1.0)**2 * rho_a * rho_b) / (1.0 - theta)**2)

def Phi_alpha_theta(theta_a, theta_b, R):
    """Blocking function for adsorption of particles in the model by Brusatori
    and Van Tassel, in terms of fractional surface coverage.
    Arguments:
        theta_a   fraction of surface covered by molecules in state a 
        theta_b   fraction of surface covered by molecules in state b
        R   (area of particle in state a)/(area of particle in state b)
    """
    theta = theta_a + theta_b
    return (1.0 - theta) * exp((-R**2 * theta_a * theta_b + \
        4 * R * theta_a * theta_b + 2 * R * theta_b**2 + \
        2 * theta_a**2 + theta_a * theta_b - R**2 * theta_b \
        - 2 * R * theta_b - 3 * theta_a) / (theta - 1.0)**2)

def Psi_alpha_beta_theta(theta_a, theta_b, r_a, r_b):
    """Blocking function for adsorption of particles in the model by Brusatori
    and Van Tassel, in terms of fractional surface coverage.
    Arguments:
        theta_a   fraction of surface covered by molecules in state a 
        theta_b   fraction of surface covered by molecules in state b
        r_a       radius of particle a
        r_b       radius of particle b
    """
    theta = theta_a + theta_b
    return exp(2 * (theta_a / r_a + theta_b / r_b) * (r_b - r_a) / (theta - 1.0) + \
            ((r_a - r_b)**2 * theta_a * theta_b / (r_a**2 * r_b**2) + \
            theta_a/r_a**2 + theta_b/r_b**2) * (r_a**2 - r_b**2)/(theta-1.0)**2)

# Langmuir-like two-stage model
def Phi_alpha_Langmuir(theta_a, theta_b, theta_max=1.0):
    """Langmuir blocking function for two-stage adsorption.
    Arguments:
        rho_a   dimensionless number density of molecules in state a
        rho_b   dimensionless number density of molecules in state b
        ratio   (area of particle in state b)/(area of particle in state a)
    """
    theta = theta_a + theta_b
    return 1.0 - theta / theta_max

### Single-stage adsorption ###
def generalized_RSA_kinetics(ka, kd, times, ASF, const_C=None, C_t=None,
        blocking_fn_args=None):
    """Arguments: 
            ka      association constant (1/sec * 1/(units of C))
            kd      dissociation constant (1/sec)
            times   time value, or array of time points (sec)
            ASF     available surface function
            const_C  near-surface bulk analyte concentration (amount/volume)
            C_t     SciPy univariate spline of near-surface concentration over
                    time (amount/volume)
            blocking_fn_args   Optional blocking_fn_args to ASF; default is None
       Returns:
            theta   array of fractional surface coverage values at each time point
    """
    from scipy.integrate import odeint

    def dtheta_dt(theta, t, ka, kd, ASF, blocking_fn_args, const_C, C_t):
        """Schaaf and Talbot 1989, eqn. 2"""
        if C_t is None:
            CA = const_C
        else:
            CA = C_t(t)[0]
        
        if blocking_fn_args is not None:
            return ka * CA * ASF(theta, blocking_fn_args) - kd * theta
        else:
            return ka * CA * ASF(theta) - kd * theta

    theta_initial = 0.0
    theta = odeint(dtheta_dt, theta_initial, times, args=(ka, kd, ASF, 
        blocking_fn_args, const_C, C_t))

    return theta.reshape(1,len(theta))[0]

def RSA_kinetics(ka, kd, times, const_C=None, C_t=None):
    """Solves eqn. 2 from Schaaf and Talbot, 1989
    Arguments:
        ka      association constant (1/sec * 1/(units of C))
        kd      dissociation constant (1/sec)
        times   time value, or array of time points (sec)
        const_C  near-surface bulk analyte concentration (amount/volume)
        C_t     SciPy univariate spline of near-surface concentration over
                time (amount/volume)
    Returns:
        theta   fractional surface coverage
    """
    return generalized_RSA_kinetics(ka, kd, times, phi_theta_fit3, 
            const_C=const_C, C_t=C_t)

def Langmuir_kinetics(ka, kd, times, const_C=None, C_t=None, theta_max=1.0):
    """Compute fractional surface coverage (theta) using a numerical solution
    to the Langmuir evolution equation.
    Arguments:
        ka      association constant (1/sec * 1/(units of C))
        kd      dissociation constant (1/sec)
        times   time value, or array of time points (sec)
        const_C  near-surface bulk analyte concentration (amount/volume)
        C_t     SciPy univariate spline of near-surface concentration over
                    time (amount/volume)
    Returns:
        theta   fractional surface coverage
    """
    return generalized_RSA_kinetics(ka, kd, times, phi_theta_Langmuir, 
            const_C=const_C, C_t=C_t, blocking_fn_args=theta_max)

def Langmuir_analytical(ka, kd, t, CA):
    """Compute fractional surface coverage (theta) using the analytical
    solution to the Langmuir evolution equation.
    Arguments:
        ka      association constant (1/sec * 1/(units of C))
        kd      dissociation constant (1/sec)
        t       time point or numpy array of times (sec)
        CA      analyte concentration near the surface (amount/volume)
    Returns:
        theta   fractional surface coverage
    """
    from numpy import exp

    return ka * CA / (ka * CA + kd) * (1.0 - exp(-t * (ka * CA + kd)))

### Two-layer Adsorption ###
def Langmuir_kinetics_twolayer_analytical(ka1, ka2, CAB_max, t, CA):
    """
    """
    CB = CAB_max * exp(-ka1 * CA * t)
    CAB = ka1 * CAB_max / (ka2 - ka1) * (exp(-ka1 * CA * t) - 
            exp(-ka2 * CA * t))
    CAAB = ka1 * CAB_max / (ka2 - ka1) * exp(-ka2 * CA * t) - \
            ka2 * CAB_max / (ka2 - ka1) * exp(-ka1 * CA * t) + CAB_max

    return CB, CAB, CAAB

def Langmuir_kinetics_twolayer(ka1, ka2, kd1, kd2, times, const_C=None, C_t=None):
    """Compute fractional surface coverage for a two-layer Langmuir adsorption
        model.
    Arguments:
        ka1     Association constant for first layer (1/s * 1/(units of CA))
        ka2     Association constant for second layer
        kd1     Association constant for first layer (1/s * 1/(units of CA))
        kd2     Association constant for second layer
        times   NumPy array of time values
        const_C     Concentration near surface
        C_t         SciPy UnivariateSpline of near-surface concentration values

    Returns:
        theta1  Fraction of surface covered by complex CAB
        theta2  Fraction of surface covered by complex CAAB
    """
    from scipy.integrate import odeint

    def dtheta_dt(thetas, t, ka1, ka2, kd1, kd2, ASF, const_C, C_t):
        """Compute d/dt[theta1, theta2] at time t"""
        if C_t is None:
            CA = const_C
        else:
            CA = C_t(t)[0]
        
        dtheta1 = ka1 * CA * (1.0 - thetas[0] - thetas[1]) - kd1 * thetas[0]\
                - ka2 * CA * thetas[0]
        dtheta2 = ka2 * CA * thetas[0] - kd2 * thetas[1]

        return dtheta1, dtheta2

    theta_initial = [0.0, 0.0]
    thetas= odeint(dtheta_dt, theta_initial, times, args=(ka1, ka2, kd1, kd2, 
        phi_theta_Langmuir, const_C, C_t))

    theta1 = thetas[:, 0]
    theta2 = thetas[:, 1]
    return theta1, theta2

def Langmuir_kinetics_twolayer_C(ka1, ka2, kd1, kd2, CB_0, times, const_C=None, 
        C_t=None):
    """Compute surface concentration for a two-layer Langmuir adsorption model.
    Arguments:
        ka1     Association constant for first layer (1/s * 1/(units of CA))
        ka2     Association constant for second layer
        kd1     Association constant for first layer (1/s * 1/(units of CA))
        kd2     Association constant for second layer
        CB_0    Initial concentration of surface sites
        times   NumPy array of time values
        const_C     Concentration near surface
        C_t         SciPy UnivariateSpline of near-surface concentration values

    Returns:
        theta1  Fraction of surface covered by complex CAB
        theta2  Fraction of surface covered by complex CAAB
    """
    from scipy.integrate import odeint

    def dYdt(Y, t, ka1, ka2, kd1, kd2, const_C, C_t):

        if C_t is None:
            CA = const_C
        else:
            CA = C_t(t)[0]

        CB = Y[0]
        CAB = Y[1]
        CAAB = Y[2]

        dCB = -ka1 * CA * CB + kd1 * CAB + kd2 * CAAB
        dCAB = ka1 * CA * CB - ka2 * CA * CAB - kd1 * CAB
        dCAAB = ka2 * CA * CAB - kd2 * CAAB
        return dCB, dCAB, dCAAB

    y = odeint(dYdt, y0=[CB_0, 0.0, 0.0], t=times, args=(ka1, ka2, kd1, kd2, const_C,
        C_t), hmax=10.0)

    CB = y[:, 0]
    CAB = y[:, 1]
    CAAB = y[:, 2]

    return CB, CAB, CAAB

### Adsorption with post-adsorption transition ###
def Van_Tassel_kinetics_reference(ka, ks, kd, r_alpha, Sigma, times, c):
    """Two-stage adsorption model derived by Brusatori and Van Tassel,
    Journal of Colloid and Interface Science 219, 333-338 (1999).
    Implemented exactly as described, except that gamma is used instead of
    rho to denote surface number density.
    
    Args:
        ka          association rate constant  (length/time)
        ks          transition rate constant  (1/time)
        kd          dissociation rate constant  (1/time)
        r_alpha     radius of particle in State 1  (length)
        Sigma       r_beta/r_alpha
        times       times (time)
        c           bulk concentration (molecules/length^3)
    Returns:
        gamma1        NumPy array of surface concentration of particles in
                        State 1 (molecules/length^2)
        gamma2        NumPy array of surface concentration of particles in 
                        State 2 (molecules/length^2)
    """
    from scipy.integrate import odeint, ode
    from numpy import exp, pi, array
    from Surface_Reaction_Tools.theoretical import Phi_alpha, Psi_alpha_beta

    def dgamma_dt(gamma_values, t, ka, ks, kd, r_alpha, Sigma, c):
        gamma_alpha = gamma_values[0];   gamma_beta = gamma_values[1];

        gamma_alpha_bar = gamma_alpha * pi * r_alpha**2
        gamma_beta_bar = gamma_beta * pi * r_alpha**2

        dgamma_beta_dt = ks * gamma_alpha * Psi_alpha_beta(gamma_alpha_bar, gamma_beta_bar*Sigma**2, Sigma**2)

        dgamma_alpha_dt = ka * c * Phi_alpha(gamma_alpha_bar, gamma_beta_bar*Sigma**2, Sigma**2) - \
                dgamma_beta_dt - kd * gamma_alpha

        return [dgamma_alpha_dt, dgamma_beta_dt]

    gamma_initial = array([0.0, 0.0])
    gamma_values = odeint(dgamma_dt, gamma_initial, times, args=(ka, ks, kd, r_alpha,
        Sigma, c))

    return gamma_values[:,0], gamma_values[:,1]

def Van_Tassel_kinetics_theta(ka, ks, kd, r_alpha, Sigma, times, c):
    """Two-stage adsorption model derived by Brusatori and Van Tassel,
    Journal of Colloid and Interface Science 219, 333-338 (1999).
    Implemented entirely in terms of fractional surface coverage.
    
    Args:
        ka          association rate constant  (length/time)
        ks          transition rate constant  (1/time)
        kd          dissociation rate constant  (1/time)
        r_alpha     radius of particle in State 1  (length)
        Sigma       r_beta/r_alpha
        times       times (time)
        c           bulk concentration (molecules/length^3)
    Returns:
        theta1        NumPy array of surface concentration of particles in
                        State 1 (molecules/length^2)
        theta2        NumPy array of surface concentration of particles in 
                        State 2 (molecules/length^2)
    """
    from scipy.integrate import odeint, ode
    from numpy import exp, pi, array
    from Surface_Reaction_Tools.theoretical import Phi_alpha_theta, Psi_alpha_beta_theta

    def dtheta_dt(theta_values, t, ka, ks, kd, r_alpha, Sigma, c):
        theta_alpha = theta_values[0];   theta_beta = theta_values[1];

        dtheta_alpha_dt = ka * c * Phi_alpha_theta(theta_alpha, theta_beta,
                1.0/Sigma) - ks * theta_alpha * Psi_alpha_beta_theta(
                        theta_alpha, theta_beta, r_alpha, r_alpha * Sigma) \
                                - kd * theta_alpha

        dtheta_beta_dt = ks * Sigma**2 * theta_alpha * Psi_alpha_beta_theta(
                theta_alpha, theta_beta, r_alpha, r_alpha * Sigma)

        return [dtheta_alpha_dt, dtheta_beta_dt]

    theta_initial = array([0.0, 0.0])
    theta_values = odeint(dtheta_dt, theta_initial, times, args=(ka, ks, kd, r_alpha,
        Sigma, c))

    return theta_values[:,0], theta_values[:,1]

def Van_Tassel_kinetics(ka, ks, kd, Sigma, times, const_C=None, 
        C_t=None):
    """Two-stage adsorption model derived by Brusatori and Van Tassel,
    Journal of Colloid and Interface Science 219, 333-338 (1999)
    SOLVED ENTIRELY IN TERMS OF FRACTIONAL SURFACE COVERAGE(S)

    Args:
        ka          association rate constant  (1/s)/(units of c)
        ks          transition rate constant  (1/s)
        kd          dissociation rate constant  (1/s)
        Sigma       ratio of r_beta/r_alpha
        times       NumPy array of times (seconds)
        const_C     bulk concentration (amount/length^3)
        C_t         SciPy univariate spline of near-surface concentration over
                        time (amount/length^3)
    Returns:
        theta1      NumPy array of fraction of surface covered by molecules in
                        state alpha
        theta2      NumPy array of fraction of surface covered by molecules in
                        state beta
    """
    return generalized_two_stage_kinetics(ka, ks, kd, Sigma**2, Phi_alpha, 
            Psi_alpha_beta, times, const_C=const_C, C_t=C_t, blocking_fn_args=Sigma**2) 

def Garcia_kinetics(k, s, a1, b, r, A_total, c, f, times):
    """Direct implementation of eqns. 1-3 from Michael et. al.
    Langmuir 2003, 19, 8033-8040
    Arguments:
        k   association rate constant   (cm^-1 s^-1)
        s   transition rate constant    (cm^-2 s^-1)
        a1  area occupied by adsorbed molecule in state 1 (cm^2)
        b   area in state 2 / area in state 1
        r   dissociation rate constant (1/s)
        A_total     total area of sensor    (cm^2)
        c   solution concentration (ng/cm^3)
        f   molecules/ng    (NA/MW/1e9)
        times   NumPy array of times (seconds)
    Returns:
        Y1  surface density of adsorbed protein in state 1 (ng/cm^2)
        Y2  surface density of adsorbed protein in state 2 (ng/cm^2)
    """
    from scipy.integrate import odeint, ode
    from numpy import exp, pi, array

    def dYdt(Y_values, t, k, s, a1, b, r, A_total, c, f):
        Y1 = Y_values[0]
        Y2 = Y_values[1]

        A_av = A_total * (1.0 - f * a1 * Y1 - f * b * a1 * Y2)
        dY1dt = (k * c - s * Y1) * A_av - r * Y1
        dY2dt = s * Y1 * A_av
        return (dY1dt, dY2dt)

    Y_initial = array([0.0, 0.0])
    Y_values = odeint(dYdt, Y_initial, times, args=(k, s, a1, b, r, A_total, c, f))

    return Y_values[:,0], Y_values[:,1]

def Langmuir_transition_kinetics(ka, ks, kd, ratio, times, const_C=None, 
    C_t=None, blocking_fn_args=1.0):
    """Two-stage adsorption model derived by Michael et al
    Langmuir 2003, 19, 8033-8040

    Args:
        ka          association rate constant  (1/s)/(units of c)
        ks          transition rate constant  (1/s)
        kd          dissociation rate constant  (1/s)
        ratio       ratio (area of stage beta/area of stage alpha)
        times       NumPy array of times (seconds)
        const_C     bulk concentration (amount/length^3)
        C_t         SciPy univariate spline of near-surface concentration over
                        time (amount/length^3)
    Returns:
        theta1      NumPy array of fraction of surface covered by molecules in
                        state alpha
        theta2      NumPy array of fraction of surface covered by molecules in
                        state beta
    """
    return generalized_two_stage_kinetics(ka, ks, kd, ratio, Phi_alpha_Langmuir, 
            Phi_alpha_Langmuir, times, const_C=const_C, C_t=C_t, 
            blocking_fn_args=blocking_fn_args) 

def generalized_two_stage_kinetics(ka, ks, kd, ratio, Phi, Psi, times, 
        const_C=None, C_t=None, blocking_fn_args=None):
    """Routine to solve two-stage adsorption models
    Args:
        ka          association rate constant  (1/s)/(units of c)
        ks          transition rate constant  (1/s)
        kd          dissociation rate constant  (1/s)
        ratio       (area of state alpha) / (area of state beta)
        Phi         Blocking function for initial adsorption
        Psi         Blocking function for transition
        times       NumPy array of times (seconds)
        const_C     bulk concentration (amount/length^3)
        C_t         SciPy univariate spline of near-surface concentration over
                        time (amount/length^3)
    Returns:
        theta1      NumPy array of fraction of surface covered by molecules in
                        state alpha
        theta2      NumPy array of fraction of surface covered by molecules in
                        state beta
    """
    from scipy.integrate import odeint, ode
    from numpy import exp, pi, array

    def dtheta_dt(theta_values, t, ka, ks, kd, ratio, Phi, Psi, const_C, C_t):
        if C_t is None:
            c = const_C
        else:
            c = C_t(t)[0]

        theta_alpha = theta_values[0];   theta_beta = theta_values[1];

        dtheta_alpha_dt = ka * c * Phi(theta_alpha, theta_beta, 
                blocking_fn_args) - ks * theta_alpha * Psi(theta_alpha, 
                        theta_beta, blocking_fn_args) - kd * theta_alpha

        dtheta_beta_dt = ks * ratio * theta_alpha * Psi(
                theta_alpha, theta_beta, blocking_fn_args)

        return [dtheta_alpha_dt, dtheta_beta_dt]

    theta_initial = array([0.0, 0.0])
    theta_values = odeint(dtheta_dt, theta_initial, times, args=(ka, ks, kd,
        ratio, Phi, Psi, const_C, C_t))

    return theta_values[:,0], theta_values[:,1]

### Diffusion-Limited Adsorption ###
def diffusion_limited_theta(phi, t, D=1.0, r=1.0):
    """Compute fractional surface coverage (theta) for the times contained in 
    t.  The boundary is assumed to be perfectly adsorbing, so the surface 
    concentration is limited only by diffusion.  Because there is no blocking
    due to adsorbed particles, the fractional surface coverage increases 
    without bound, instead of reaching saturation at 1.

    phi   volume fraction
    t     times (single value or NumPy array)
    D     diffusion coefficient (optional, default=1.0)
    r     particle radius (optional, default=1.0)
    """
    from numpy import sqrt, pi
    return 3.0 * r * phi / (2 * pi) * sqrt(pi * D * t)

def diffusion_limited_flux(z, t, c0, D=1.0):
    """Compute flux vs.time J(z,t) for diffusion-limited
    adsorption on a perfect absorbing boundary.
    
    Arguments for dimensioned flux in mol/m^2/sec :
        z       z location (m)
        t       numpy array of times (s)
        c0      bulk concentration (mol/m^3)
        D       diffusion coefficient (optional, default=1.0)

    For non-dimensioned flux:
        z       z location (non-dim)
        t       numpy array of times (non-dim)
        c0      bulk number density (non-dim)
    """
    from numpy import sqrt, pi, exp
    return -c0 * D / sqrt(pi * D * t) * exp(-z**2 / (4 * D * t))
