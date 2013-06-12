#!/usr/bin/env python
import unittest
import numpy
from matplotlib import pyplot as plt

from colloid.adsorption.theoretical import *

class TestTheoreticalModule(unittest.TestCase):
    """Unit tests for module theoretical"""
    NA = 6.022e23

    def setUp(self):
        self.molecular_weight = 160e3
        self.f = self.NA / self.molecular_weight * 1e-9

        dt = 0.1
        Tmax = 60.0

        self.ka = 6e5
        self.kd = 1e-6
        self.CA = 1e-6
        self.times = numpy.arange(0.0, Tmax, dt)

        self.ka1 = 1e6        # cm^3/ng/s
        self.ka2 = 1e5        # cm^3/ng/s
        self.kd1 = 0.0      # 1/s
        self.kd2 = 0.0      # 1/s
        self.area = 2.713e-12   # cm^2

    def test_Langmuir_kinetics(self):
        """Verify that analytical solution and numerical solution to the
        Lanmguir equation return the same result
        """

        analytical_solution = Langmuir_analytical(self.ka, self.kd, self.CA, 
                self.times)
        numerical_solution = Langmuir_kinetics(self.ka, self.kd, self.times, 
                const_C=self.CA, theta_max=1.0)

        self.assertTrue(numpy.allclose(analytical_solution, 
            numerical_solution))

    def test_RSA_kinetics(self):
        """Verify that function RSA_kinetics returns the same result as
        generalized_RSA_kinetics.
        """
        solution1 = RSA_kinetics(self.ka, self.kd, self.times, const_C=self.CA)
        solution2 = generalized_RSA_kinetics(self.ka, self.kd, self.times, 
                phi_theta_fit3, const_C=self.CA)

        self.assertTrue(numpy.allclose(solution1, solution2))

    def test_Langmuir_kinetics_twolayer(self):
        """Verify that the numerical two-layer model returns the same result
        as the analytical solution to the same model when dissociation
        is neglected."""
        theta1, theta2 = Langmuir_kinetics_twolayer(self.ka1, self.ka2, self.kd1, 
                    self.kd2, self.times, const_C=self.CA)
        solution1 = theta1 + 2 * theta2

        CAB_max = 1.0 / self.area / self.f
        CB, CAB, CAAB = Langmuir_kinetics_twolayer_analytical(self.ka1, self.ka2,
                CAB_max, self.times, self.CA)
        solution2 = (CAB + 2 * CAAB) * self.area * self.f
        
#        plt.figure()
#        plt.plot(self.times, solution1)
#        plt.plot(self.times, solution2)
#        plt.title("Two-layer Langmuir model")
#        plt.show()

        self.assertTrue(numpy.allclose(solution1, solution2))

    def test_Van_Tassel_reference(self):
        """ """
        print("""Verify that the plots produced by this routine are equivalent
        to figures 2 and 3 from Brusatori and Van Tassel, Journal of Colloid
        and Interface Science 219, 333-338 (1999)""")

        r_alpha = 1.0
        Sigma = 1.2
        Ks = 0.1
        Kd = 0.1
        ka = 1.0
        c = 1.0
        ks = Ks * ka * c * numpy.pi * r_alpha**2
        kd = Kd * ka * c * numpy.pi * r_alpha**2
        nd_times = numpy.arange(0.0, 16.0, 0.1)
        times = nd_times/(ka*c*numpy.pi*r_alpha**2)

        gamma_alpha, gamma_beta = Van_Tassel_kinetics_reference(ka, ks, kd, r_alpha, Sigma,
                times, c)

        # Figure 2
        plt.figure()
        plt.subplot(2, 1, 1)
        plt.plot(nd_times, gamma_alpha * numpy.pi * r_alpha**2, label='State 1')
        plt.plot(nd_times, gamma_beta * numpy.pi * r_alpha**2, label='State 2')
        plt.plot(nd_times, (gamma_alpha + gamma_beta) * numpy.pi * r_alpha**2, label='Total')
        plt.axis([0.0, 16.0, 0.0, 0.6])
        plt.legend(loc='upper right')
        plt.title('Figure 2')

        Sigma = 1.5
        gamma_alpha, gamma_beta = Van_Tassel_kinetics_reference(ka, ks, kd, r_alpha, Sigma,
                times, c)

        plt.subplot(2, 1, 2)
        plt.plot(nd_times, gamma_alpha * numpy.pi * r_alpha**2, label='State 1')
        plt.plot(nd_times, gamma_beta * numpy.pi * r_alpha**2, label='State 2')
        plt.plot(nd_times, (gamma_alpha + gamma_beta) * numpy.pi * r_alpha**2, label='Total')
        plt.axis([0.0, 16.0, 0.0, 0.6])
        plt.xlabel('t')

        # Figure 3
        Ks = 0.01
        ks = Ks * ka * c * numpy.pi * r_alpha**2
        Sigma = 1.2

        plt.figure()
        plt.subplot(2, 1, 1)
        plt.plot(nd_times, gamma_alpha * numpy.pi * r_alpha**2, label='State 1')
        plt.plot(nd_times, gamma_beta * numpy.pi * r_alpha**2, label='State 2')
        plt.plot(nd_times, (gamma_alpha + gamma_beta) * numpy.pi * r_alpha**2, label='Total')
        plt.axis([0.0, 8.0, 0.0, 0.6])
        plt.legend(loc='upper right')
        plt.title('Figure 3')

        Ks = 0.1
        Kd = 0.01
        ks = Ks * ka * c * numpy.pi * r_alpha**2
        kd = Kd * ka * c * numpy.pi * r_alpha**2

        gamma_alpha, gamma_beta = Van_Tassel_kinetics_reference(ka, ks, kd, r_alpha, Sigma,
                times, c)

        plt.subplot(2, 1, 2)
        plt.plot(nd_times, gamma_alpha * numpy.pi * r_alpha**2, label='State 1')
        plt.plot(nd_times, gamma_beta * numpy.pi * r_alpha**2, label='State 2')
        plt.plot(nd_times, (gamma_alpha + gamma_beta) * numpy.pi * r_alpha**2, label='Total')
        plt.axis([0.0, 8.0, 0.0, 0.6])
        plt.xlabel('t')

        plt.show()
        self.assertTrue(True)

    def test_Van_Tassel_kinetics(self):
        """ """
        from colloid.adsorption.theoretical import Van_Tassel_kinetics, \
        Van_Tassel_kinetics_reference, Van_Tassel_kinetics_theta

        ### Parameters ###
        c = 100 * 1000            # ng/cm^3

        times = numpy.arange(0.0, 5000.0, 1.0)

        ka = 1e-5       # cm/s
        ks = 0.05        # 1/s
        kd = 1e-3      # 1/s
        r_alpha = 2e-7      # cm
        r_beta = 5.5e-7      # cm
        Sigma = r_beta / r_alpha

        # Reference implementation
        gamma_alpha_ref, gamma_beta_ref = Van_Tassel_kinetics_reference(ka, ks,
                kd, r_alpha, Sigma, times, c * self.f)

#        plt.figure()
#        plt.plot(times, gamma_alpha_ref/self.NA, label='State 1')
#        plt.plot(times, gamma_beta_ref/self.NA, label='State 2')
#        plt.plot(times, (gamma_alpha_ref + gamma_beta_ref)/self.NA, label='Total')
#        plt.legend(loc='best')
#        plt.title('Reference')
#        axis_limits = plt.axis()

        # New implementation purely in terms of theta
        theta_alpha, theta_beta = Van_Tassel_kinetics_theta(
                ka * numpy.pi * r_alpha**2, ks, kd, r_alpha, Sigma, times, c * self.f)
        gamma_alpha_theta = theta_alpha/(numpy.pi*r_alpha**2)
        gamma_beta_theta = theta_beta/(numpy.pi*r_beta**2)

#        plt.figure()
#        plt.plot(times, gamma_alpha_theta/self.NA, label='State 1')
#        plt.plot(times, gamma_beta_theta/self.NA, label='State 2')
#        plt.plot(times, (gamma_alpha_theta + gamma_beta_theta)/self.NA, label='Total')
#        plt.legend(loc='best')
#        plt.axis(axis_limits)
#        plt.title('All theta')

        # New implementation using generalized two-stage kinetics
        theta_alpha, theta_beta = Van_Tassel_kinetics(
                ka * numpy.pi * r_alpha**2, ks, kd, Sigma, times, c * self.f)
        gamma_alpha = theta_alpha/(numpy.pi*r_alpha**2)
        gamma_beta = theta_beta/(numpy.pi*r_beta**2)

#        plt.figure()
#        plt.plot(times, gamma_alpha/self.NA, label='State 1')
#        plt.plot(times, gamma_beta/self.NA, label='State 2')
#        plt.plot(times, (gamma_alpha + gamma_beta)/self.NA, label='Total')
#        plt.legend(loc='best')
#        plt.axis(axis_limits)
#        plt.title('Generalized two-stage')

        self.assertTrue(numpy.allclose(gamma_alpha_ref + gamma_beta_ref,
            gamma_alpha_theta + gamma_beta_theta) & numpy.allclose(
                gamma_alpha_ref + gamma_beta_ref, gamma_alpha + gamma_beta))
#        plt.show()

    def test_Langmuir_transition_kinetics(self):
        """Verify that the Langmuir adsorption model with post-adsorption
        transition gives the same result as the direct implementation of
        the same model from Langmuir 2003, 19, 8033-8040"""

        from colloid.adsorption.theoretical import Langmuir_transition_kinetics

        # Direct implementation of model published by Garcia
        M = 39e3   # g/mole
        A_total = 0.0075   # cm^2
        f = self.NA/M/1e9
        c = 10 * 1000   # ng/cm^2
        times = arange(0.0, 120.0, 1.0) # sec
        k = 0.19    # /cm /sec
        s = 7.0     # /cm^2 /sec
        a1 = 0.33e-14   # cm^2/molecule
        b = 150.0   # ratio
        r = 0.60    # 1/s

        Y1, Y2 = Garcia_kinetics(k, s, a1, b, r, A_total, c, f, times)
        ka = k * f * a1 * A_total
        ks = s * A_total
        kd = r
        ratio = b
        theta1, theta2 = Langmuir_transition_kinetics(ka, ks, kd, ratio, times,
                const_C=c)
        self.assertTrue(numpy.allclose(theta1 / a1 / f + theta2 / (a1 * b) / f,
            Y1 + Y2))

if __name__ == '__main__':
    unittest.main()
