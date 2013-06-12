#!/usr/bin/env python
from colloid.adsorption.theoretical import RSA_kinetics, Langmuir_kinetics,\
        Langmuir_kinetics_twolayer, Van_Tassel_kinetics, Van_Tassel_kinetics
import numpy as np

dashes = [(12,6), (12,6,4,6), (24,6), (8,8), (2,4)]
colors = ['k', 'b', 'g', 'r', 'm', 'c']
NA = 6.022e23

### Data handling ###
def load_average_data(fileName, reduce=False, max_points=50):
    """Load experimental data from an HDF5 file."""
    import tables
    from Data_Tools.process import rebin

    def get_conc(group):
        """Used to sort list of concentration groups by concentration"""
        return group._v_attrs.concentration

    h5file = tables.openFile(fileName, mode="r")
    
    concentrations = []
    exp_labels = []
    exp_surf_conc = []
    exp_std = []
    exp_times = []

    groupList = h5file.listNodes(h5file.root)
    groupList.sort(key=get_conc, reverse=False)    # List from HDF5 file is in no particular order

    for group in groupList:
        print group
        exp_labels.append(group._v_title)

        concentrations.append(group._v_attrs.concentration)   # ug/ml
        t = group.average_time.read()
        avg = group.average.read()
        std = group.standard_deviation.read()

        if reduce and len(t) > max_points:
            t = rebin(t, (max_points,))
            avg = rebin(avg, (max_points,))
            std = rebin(std, (max_points,))

        exp_times.append(t)
        exp_surf_conc.append(avg)
        exp_std.append(std)

        runList = h5file.listNodes(group)

    h5file.close()
    return np.array(concentrations), exp_times, exp_surf_conc, exp_std, exp_labels

def load_experimental_data(fileName):
    """Load experimental data from an HDF5 file."""
    import tables
    
    def get_conc(group):
        """Used to sort list of concentration groups by concentration"""
        return group._v_attrs.concentration

    def get_run_number(run):
        return run._v_attrs.run

    h5file = tables.openFile(fileName, mode="r")
    
    concentrations = []
    exp_labels = []
    exp_surf_conc = []
    exp_times = []
    exp_injection_times = []

    groupList = h5file.listNodes(h5file.root)
    groupList.sort(key=get_conc, reverse=False)    # List from HDF5 file is in no particular order

    # Iterate over concentration groups
    for i in range(len(groupList)):
        group = groupList[i]
        exp_labels.append([])

        concentrations.append(group._v_attrs.concentration)   # ug/ml

        # Get run groups
        runList = h5file.listNodes(group, classname="Group")
        runList.sort(key=get_run_number)    # sort by run number

        times = []
        surf_conc = []
        injection_times = []

        for run in runList:
            print "    " + str(run)
            surf_conc.append(run.surface_concentration.read())
            times.append(run.time.read())
            injection_times.append(run._v_attrs.proteinInjectionTime)
            exp_labels[i].append(group._v_title + " " + run._v_title)

        exp_surf_conc.append(surf_conc)
        exp_times.append(times)
        exp_injection_times.append(injection_times)

    h5file.close()

    # Eliminate data before the injection of protein, which occurs at t=0
    for i in range(len(concentrations)):
        for run in range(len(exp_surf_conc[i])):
            injection_index = abs(exp_times[i][run]).argmin()   # Find t=0
            exp_surf_conc[i][run] = exp_surf_conc[i][run][injection_index:]
            exp_times[i][run] = exp_times[i][run][injection_index:]
    
    return np.array(concentrations), exp_times, exp_surf_conc, exp_labels

def saturation_surf_conc(exp_surf_conc, saturation_criteria):
    """Quantitatively computes the saturation surface concentration.

    Arguments:
    exp_surf_conc       list of arrays of average surface concentration over time
    saturation_criteria     fraction of maximum value at which saturation will be
                            assumed to have occurred (0 to 1.0)
    """
    saturation_surf_conc = []
    saturation_start_times = []
    
    for i in range(len(exp_surf_conc)):
        # Saturation is assumed to begin when surf. conc. reaches a certain fraction of its maximum value
        saturation_indices, = np.where(exp_surf_conc[i] >= exp_surf_conc[i].max()*saturation_criteria)
        saturation_start_times.append(saturation_indices.min())
        saturation_surf_conc.append(exp_surf_conc[i][saturation_start_times[i]:].mean())

    return saturation_surf_conc

def average_kinetic_data(surf_conc, exp_times, delta_t):
    """Average two or more data sets.

    Arguments:
    surf_conc       list of arrays.  Each array is one run.
    exp_times       list of arrays of time points, each corresponding to
                    one of the arrays in surf_conc
    delta_t         time interval for constructing time array for the average

    Returns:
    interp_times    array of times
    avg             array of average values
    stdev           array of standard deviations
    """
    from scipy.interpolate import UnivariateSpline
    from numpy import arange, empty, mean, std
    
    data = []
    max_time = 0.0
    max_times = []

    for i in range(len(surf_conc)):
        if exp_times[i].max() > max_time:
            max_time = exp_times[i].max()
        data.append(UnivariateSpline(exp_times[i], surf_conc[i], k=1, s=0))
        max_times.append(exp_times[i].max())

    # Summarize
    interp_times = arange(0.0, max_time+delta_t/10, delta_t)
    avg = empty(len(interp_times), dtype=float)
    stdev = empty(len(interp_times), dtype=float)

    summarized_data = []
    for i in range(len(interp_times)):
        t = interp_times[i]
        
        temp = []
        for run in range(len(surf_conc)):
            if t <= max_times[run]:
                temp.append(data[run](t))

        avg[i] = mean(temp)

        if len(temp) > 1:
            stdev[i] = std(temp)
        else:
            stdev[i] = 0.0

    return interp_times, avg, stdev

### Master fitting class ###
class Fitter():
    """A class for fitting adsorption models to adsorption data."""

    def load_CFD_concentrations(self):
        from Data_Tools.process import rebin
        import os
        from scipy.interpolate import UnivariateSpline
        from colloid.adsorption.data_utilities import read_ACE_near_surface_concentration
        from matplotlib import pyplot as plt

        self.cfd_concentrations = []
        for i in range(len(self.concentrations)):
            (cfd_conc_times, cfd_conc) = read_ACE_near_surface_concentration(self.cfd_path +
                    os.sep + self.cfd_fileNames[i])
            cfd_conc = cfd_conc*self.molecular_weight*1e6        # mol/L to ng/cm^3

            if len(cfd_conc_times) > 10000:
                c = rebin(cfd_conc, (1000,))
                t = rebin(cfd_conc_times, (1000,))
            else:
                c = cfd_conc.copy()
                t = cfd_conc_times.copy()

            # Calculate linear approach to bulk concentration
            slope = (c[-1] - c[-100]) / (t[-1] - t[-100])
            t_sat = t[-1] + (self.concentrations[i] - c[-1]) / slope
            t = np.concatenate([t, np.array([t_sat])])
            c = np.concatenate([c, np.array([self.concentrations[i]])])

            # Append final point at bulk concentration
            t = np.concatenate([t, np.array([10 * t_sat])])
            c = np.concatenate([c, np.array([self.concentrations[i]])])

            spline = UnivariateSpline(t, c, k=1, s=0)
            self.cfd_concentrations.append(spline)

    def plot_CFD_concentrations(self):
        from matplotlib import pyplot as plt

        plt.figure()
        for i in range(len(self.cfd_concentrations)):
            plt.plot(self.times[i], self.cfd_concentrations[i](self.times[i]))
            plt.axhline(y=self.concentrations[i], color='k', ls='--')
        plt.xlabel('time/sec')
        plt.ylabel('ng/cm^3')
        plt.title('CFD near-surface concentrations')
        plt.show()

    def load_experimental_data(self):
        (self.concentrations, self.times, self.mean_surf_conc, self.exp_std,
               self.exp_labels) = load_average_data(self.experimentalFileName)

    def compute_error(self, parameters, concentrations, exp_surf_conc, exp_times):
        error = self.compute_error_vector(parameters, concentrations,
                exp_surf_conc, exp_times)
        return sum(pow(error,2))/len(error)

    def fit_adsorption_kinetics(self):
        from scipy.optimize import leastsq

        errors = self.compute_error_vector(self.initial_parameters, 
                self.concentrations, self.mean_surf_conc, self.times)
        self.initial_sse = sum(errors**2)/len(errors)

        (self.fitted_parameters, cov_x, infodict, mesg, ier) = leastsq(
                self.compute_error_vector, self.initial_parameters, 
                args=self.model_parameters, full_output=True)

        errors = self.compute_error_vector(self.fitted_parameters,  
                self.concentrations, self.mean_surf_conc, self.times)
        self.final_sse = sum(errors**2)/len(errors)

    def fit_adsorption_kinetics_slsqp(self):
        from scipy.optimize import fmin_slsqp

        errors = self.compute_error_vector(self.initial_parameters, 
                self.concentrations, self.mean_surf_conc, self.times)
        self.initial_sse = sum(errors**2)/len(errors)

        (self.fitted_parameters, cov_x, infodict, mesg, ier) = fmin_slsqp(
                self.compute_error, self.initial_parameters, 
                args=self.model_parameters, full_output=True)

        errors = self.compute_error_vector(self.fitted_parameters,  
                self.concentrations, self.mean_surf_conc, self.times)
        self.final_sse = sum(errors**2)/len(errors)

    def fit_adsorption_kinetics_fmin_l_bgfs_b(self):
        from scipy.optimize import fmin_l_bfgs_b

        error = self.compute_error(self.initial_parameters, 
                self.concentrations, self.mean_surf_conc, self.times)
        self.initial_sse = error

        (self.fitted_parameters, self.final_sse, infodict) = fmin_l_bfgs_b(
                self.compute_error, self.initial_parameters, 
                args=self.model_parameters, approx_grad=True, 
                bounds=[(0.0, None), (0.0, None), (0.0, None)], 
                epsilon=min(self.initial_parameters)/50,iprint=2)

    def fit_adsorption_kinetics_fmin_cobyla(self):
        from scipy.optimize import fmin_cobyla

        def constr0(params):
            return params[0]
        def constr1(params):
            return params[1]
        def constr2(params):
            return params[2]

        error = self.compute_error(self.initial_parameters, 
                self.concentrations, self.mean_surf_conc, self.times)
        self.initial_sse = error

        self.fitted_parameters = fmin_cobyla(
                self.compute_error, self.initial_parameters, 
                [], args=self.model_parameters,
                consargs=(), rhobeg=min(self.initial_parameters), rhoend=1e-9,
                iprint=2)

        error = self.compute_error(self.fitted_parameters,  
                self.concentrations, self.mean_surf_conc, self.times)
        self.final_sse = error

    def fit_adsorption_kinetics_brute(self):
        from scipy.optimize import brute

        error = self.compute_error(self.initial_parameters, 
                self.concentrations, self.mean_surf_conc, self.times)
        self.initial_sse = error

        (self.fitted_parameters, self.final_sse, self.grid, self.f_values) = brute(
                self.compute_error, self.search_limits, 
                args=self.model_parameters, full_output=True)

    def plot_offset(self, offset_step=50, label_positions=[]):
        from matplotlib import pyplot as plt
        offset = 0.0

        plt.figure(facecolor='white')
        ax = plt.axes(frameon=False)
        ax.arrow(0.0, 0.0, 10000, 0, color='k', linewidth=2.0)

        for i in range(self.start_conc, self.end_conc):
            label = '$' + str(self.concentrations[i]/1000) + '\, \mu g/ml$'

            results = self.total_adsorbed_protein(i)
            ax.plot(self.times[i], results + offset, color='k', dashes=dashes[1], 
                    antialiased=True, linewidth=2, label=label)

            if len(label_positions) == 0:
                ax.text(self.times[i][-1] + 50, results[-1] + offset - 5, label, fontsize='large')
            else:
                ax.text(label_positions[i][0], label_positions[i][1], label, fontsize='large')

            label = '$' + str(self.concentrations[i]/1000) + '\, \mu g/ml \; experiment$'

            ax.plot(self.times[i], self.mean_surf_conc[i] + offset, color='k',
                    ls='-', label=label)
            ax.plot(self.times[i], self.mean_surf_conc[i] + self.exp_std[i] +
                    offset, color='k', 
                    linestyle=':')
            ax.plot(self.times[i], self.mean_surf_conc[i] - self.exp_std[i] + 
                    offset, color='k',
                    linestyle=':')
            offset += offset_step

        # scale bar in place of y axis
        ax.text(-800, 40, r"$100 \, ng/cm^2$", fontsize='large')
        ax.arrow(-50, 0, 0, 100, linewidth=3)
        plt.xlim(xmin=-850)

        # For PEG/OEG
#        # scale bar in place of y axis
#        ax.text(-600, 4, r"$10 \, ng/cm^2$", fontsize='large')
#        ax.arrow(-50, 0, 0, 10, linewidth=3)
#        plt.xlim(xmin=-650)

        # Tick formatting
        plt.ylim(ymin=0)
        locs = ax.xaxis.get_ticklocs()
        new_locs = []
        for loc in locs:
            if loc >= 0:
                new_locs.append(loc)
        ax.xaxis.set_ticks(new_locs)

        ax.set_yticks([])   # Turn of all y ticks
        for tick in ax.xaxis.get_major_ticks():     # Turn off upper x ticks
            tick.tick2On = False

        plt.xlabel(r"time $(s)$")

    def plot(self):
        from matplotlib import pyplot as plt

        plt.figure(facecolor='white')
        ax = plt.axes()

        for i in range(self.start_conc, self.end_conc):
            label = '$' + str(self.concentrations[i]/1000) + '\, \mu g/ml$'

            results = self.total_adsorbed_protein(i)
            ax.plot(self.times[i], results, color=colors[i], ls='-', 
                    antialiased=True, linewidth=2, label=label)

            label = '$' + str(self.concentrations[i]/1000) + '\, \mu g/ml \, exp$'

            ax.plot(self.times[i], self.mean_surf_conc[i], color=colors[i],
                    dashes=dashes[i], antialiased=True)
            ax.plot(self.times[i], self.mean_surf_conc[i] + self.exp_std[i], color=colors[i], 
                    linestyle=':')
            ax.plot(self.times[i], self.mean_surf_conc[i] - self.exp_std[i], color=colors[i],
                    linestyle=':')

        plt.ylim(ymin=0)
        plt.xlabel(r"time $(s)$")
        plt.ylabel(r"adsorbed density $(ng/cm^2)$")
        plt.legend(loc='upper right', handlelength=4, labelspacing=0)

class Fit_RSA(Fitter):
    """Fit the RSA model to experimental data."""
    def compute_error_vector(self, parameters, concentrations, exp_surf_conc, exp_times):
        self.ka = parameters[0]
        self.kd = parameters[1]
        self.r = parameters[2]

        error = np.array([])
        for i in range(self.start_conc, self.end_conc):
            error = np.concatenate([error, self.run_adsorption_model(i) - exp_surf_conc[i]])
        return error


    def run_adsorption_model(self, i):
        if self.constant_concentration:
            return RSA_kinetics(self.ka, self.kd, self.times[i], 
                    const_C=self.concentrations[i]) / (self.f * np.pi * self.r**2)
        else:
            return RSA_kinetics(self.ka, self.kd, self.times[i], 
                    C_t=self.cfd_concentrations[i]) / (self.f * np.pi * self.r**2)

    def total_adsorbed_protein(self, i):
        rho = self.run_adsorption_model(i)
        return rho

    def print_parameters(self, parameters):
        print("ka = {0:.3e} cm^3/ng/s".format(parameters[0]))
        print("kd = {0:.3e} 1/s".format(parameters[1]))
        print("r = {0:.3e} m A = {1:.0f} nm^2".format(parameters[2] / 100,
            np.pi * (parameters[2]/ 100)**2 * 1e18))

#    def plot(self):
#        from matplotlib import pyplot as plt
#
#        Fitter.plot(self)
#        plt.axis([0.0, max(self.times[0].max(), self.times[1].max()), 0.0, 160.0])
#        plt.xlabel(r"time $(s)$")
#        plt.ylabel(r"adsorbed density $(ng/cm^2)$")
#        plt.title('RSA model: ' + self.plot_title)

class Fit_Langmuir(Fitter):
    """Fit the Langmuir model to experimental data."""
    def compute_error_vector(self, parameters, concentrations, exp_surf_conc, exp_times):
        self.ka = parameters[0]
        self.kd = parameters[1]
        self.area = parameters[2]
#        self.theta_max = parameters[3]     # Requires constrained opt.

        error = np.array([])
        for i in range(self.start_conc, self.end_conc):
            error = np.concatenate([error, self.run_adsorption_model(i) - exp_surf_conc[i]])
        return error

    def run_adsorption_model(self, i):
        if self.constant_concentration:
            theta = Langmuir_kinetics(self.ka, self.kd, self.times[i], 
                    const_C=self.concentrations[i], theta_max=self.theta_max)
        else:
            theta = Langmuir_kinetics(self.ka, self.kd, self.times[i], 
                    C_t=self.cfd_concentrations[i], theta_max = self.theta_max)
        return theta / self.area / self.f

    def total_adsorbed_protein(self, i):
        rho = self.run_adsorption_model(i)
        return rho

    def print_parameters(self, parameters):
        print("ka = {0:.3e} cm^3/ng/s".format(parameters[0]))
        print("kd = {0:.3e} 1/s".format(parameters[1]))
        print("A = {0:.0f} nm^2".format(parameters[2] * 1e14))

    def plot(self):
        from matplotlib import pyplot as plt

        Fitter.plot(self)
        plt.axis([0.0, max(self.times[0].max(), self.times[1].max()), 0.0, 160.0])
        plt.xlabel(r"time $(s)$")
        plt.ylabel(r"adsorbed density $(ng/cm^2)$")
        plt.title('Langmuir model: ' + self.plot_title)
        plt.legend(loc='best')

class Fit_Langmuir_TwoLayer(Fitter):
    """Fit the Langmuir model to experimental data."""
    def compute_error_vector(self, parameters, concentrations, exp_surf_conc, exp_times):
        self.ka1 = parameters[0]
        self.ka2 = parameters[1]
        self.kd1 = parameters[2]
        self.kd2 = parameters[3]
        self.area = parameters[4]
        
        if self.kd1 < 0:
            self.kd1 = 0.0
        if self.kd2 < 0:
            self.kd2 = 0.0

        error = np.array([])
        for i in range(self.start_conc, self.end_conc):
            error = np.concatenate([error, self.total_adsorbed_protein(i) - exp_surf_conc[i]])
        
        if parameters[2] < 0 or parameters[3] < 0:
            error *= 100
        return error

    def run_adsorption_model(self, i):
        if self.constant_concentration:
            theta1, theta2 = Langmuir_kinetics_twolayer(self.ka1, self.ka2, self.kd1,
                    self.kd2, self.times[i], const_C=self.concentrations[i])
        else:
            theta1, theta2 = Langmuir_kinetics_twolayer(self.ka1, self.ka2, self.kd1, 
                    self.kd2, self.times[i], C_t=self.cfd_concentrations[i])
        return theta1 / self.area / self.f, theta2 / self.area / self.f

    def total_adsorbed_protein(self, i):
        rho1, rho2 = self.run_adsorption_model(i)
        return rho1 + 2 * rho2

    def active_adsorbed_protein(self, i):
        rho1, rho2 = self.run_adsorption_model(i)
        return rho1 + rho2

    def print_parameters(self, parameters):
        print("ka1 = {0:.3e} cm^3/ng/s".format(parameters[0]))
        print("ka2= {0:.3e} cm^3/ng/s".format(parameters[1]))

        print("kd1 = {0:.3e} 1/s".format(parameters[2]))
        print("kd2 = {0:.3e} 1/s".format(parameters[3]))

        print("A = {0:.0f} nm^2".format(parameters[4] * 1e14))

    def plot_mechanism(self, i):
        from matplotlib import pyplot as plt
        rho1, rho2 = self.run_adsorption_model(i)

        plt.figure()
        plt.plot(self.times[i], rho1, label='Layer 1')
        plt.plot(self.times[i], 2 * rho2, label='Layer 2')
        plt.plot(self.times[i], rho1 + 2 * rho2, label='Total')

        rho_active = self.active_adsorbed_protein(i)
        plt.plot(self.times[i], rho_active, label='Active')

        plt.axis([0.0, max(self.times[0].max(), self.times[1].max()), 0.0, 160.0])
        plt.xlabel(r"time $(s)$")
        plt.ylabel(r"adsorbed density $(ng/cm^2)$")
        plt.title('Two layer Langmuir model: ' + self.plot_title)

        plt.legend(loc='best')

    def plot(self):
        from matplotlib import pyplot as plt

        Fitter.plot(self)
        plt.axis([0.0, max(self.times[0].max(), self.times[1].max()), 0.0, 160.0])
        plt.xlabel(r"time $(s)$")
        plt.ylabel(r"adsorbed density $(ng/cm^2)$")
        plt.title('Two-layer Langmuir model: ' + self.plot_title)

def Langmuir_equilibrium(t, ka, kd, Ca, Pmax):
    return ka*Ca*Pmax(ka*Ca+kd)

### Two-stage adsorption ###
class Fit_Langmuir_Transition(Fitter):
    """Fit a two-stage adsorption model, based on the Langmuir blocking
    function, to experimental data."""

    def compute_error_vector(self, parameters, concentrations, exp_surf_conc, exp_times):
        self.ka = parameters[0]/self.scaling_factors[0]
        self.ks = parameters[1]/self.scaling_factors[1]
        self.kd = parameters[2]/self.scaling_factors[2]
        self.area_1 = parameters[3]/self.scaling_factors[3]
        self.ratio = parameters[4]/self.scaling_factors[4]

        error = np.array([])
        for i in range(self.start_conc, self.end_conc):
            error = np.concatenate([error, self.total_adsorbed_protein(i) - exp_surf_conc[i]])

        return error

    def run_adsorption_model(self, i):
        from colloid.adsorption.theoretical import Langmuir_transition_kinetics
        if self.constant_concentration:
            theta1, theta2 = Langmuir_transition_kinetics(self.ka, self.ks, self.kd,
                    self.ratio, self.times[i], const_C=self.concentrations[i], 
                    blocking_fn_args=self.theta_max)
        else:
            theta1, theta2 = Langmuir_transition_kinetics(self.ka, self.ks, self.kd, 
                    self.ratio, self.times[i], C_t=self.cfd_concentrations[i],
                    blocking_fn_args=self.theta_max)
        rho1 = theta1 / self.area_1 / self.f
        rho2 = theta2 / (self.area_1 * self.ratio) / self.f
        return rho1, rho2

    def total_adsorbed_protein(self, i):
        rho1, rho2 = self.run_adsorption_model(i)
        return rho1 + rho2

    def active_adsorbed_protein(self, i):
        rho1, rho2 = self.run_adsorption_model(i)
        return rho1

    def print_parameters(self, parameters):
        print("ka = {0:.3e} cm^3/ng/s".format(parameters[0]/self.scaling_factors[0]))
        print("ks= {0:.3e} 1/s".format(parameters[1]/self.scaling_factors[1]))
        print("kd = {0:.3e} 1/s".format(parameters[2]/self.scaling_factors[2]))
        area_1 = parameters[3] / self.scaling_factors[3]
        print("area 1 = {0:.0f} nm^2".format(self.area_1 * 1e14))
        print("ratio = {0:.2f}".format(parameters[4]/self.scaling_factors[4]))
        area_2 = parameters[3] / self.scaling_factors[3] * parameters[4] / self.scaling_factors[4]
        print("area 2 = {0:.0f} nm^2".format(area_2 * 1e14))

    def plot_mechanism(self, i):
        from matplotlib import pyplot as plt

        rho1, rho2 = self.run_adsorption_model(i)

        plt.figure()
        plt.plot(self.times[i], rho1, label='State 1')
        plt.plot(self.times[i], rho2, label='State 2')
        plt.plot(self.times[i], rho1 + rho2, label='Total')

        plt.axis([0.0, max(self.times[0].max(), self.times[1].max()), 0.0, 160.0])
        plt.xlabel(r"time $(s)$")
        plt.ylabel(r"adsorbed density $(ng/cm^2)$")
        plt.title('Langmuir two-stage model: ' + self.plot_title)

        plt.legend(loc='best')
    pass

class Fit_VanTassel(Fitter):
    """Fit the VanTassel multi-stage model to experimental data."""
    def compute_error_vector(self, parameters, concentrations, exp_surf_conc, exp_times):
        self.ka = parameters[0]/self.scaling_factors[0]
        self.ks = parameters[1]/self.scaling_factors[1]
        self.kd = parameters[2]/self.scaling_factors[2]
        self.r_alpha = parameters[3]/self.scaling_factors[3]
        self.Sigma = parameters[4]/self.scaling_factors[4]

        error = np.array([])
        for i in range(self.start_conc, self.end_conc):
            error = np.concatenate([error, self.total_adsorbed_protein(i) - exp_surf_conc[i]])

        return error

    def run_adsorption_model(self, i):
        if self.constant_concentration:
            theta1, theta2 = Van_Tassel_kinetics(self.ka, self.ks, self.kd,
                    self.Sigma, self.times[i],
                    const_C=self.concentrations[i])
        else:
            theta1, theta2 = Van_Tassel_kinetics(self.ka, self.ks, self.kd, 
                    self.Sigma, self.times[i],
                    C_t=self.cfd_concentrations[i])
        rho1 = theta1 / (np.pi * self.r_alpha**2) / self.f
        rho2 = theta2 / (np.pi * (self.r_alpha * self.Sigma)**2) / self.f
        return rho1, rho2

    def total_adsorbed_protein(self, i):
        rho1, rho2 = self.run_adsorption_model(i)
        return rho1 + rho2

    def active_adsorbed_protein(self, i):
        rho1, rho2 = self.run_adsorption_model(i)
        return rho1

    def print_parameters(self, parameters):
        print("ka = {0:.3e} cm^3/ng/s".format(parameters[0]/self.scaling_factors[0]))
        print("ks= {0:.3e} 1/s".format(parameters[1]/self.scaling_factors[1]))
        print("kd = {0:.3e} 1/s".format(parameters[2]/self.scaling_factors[2]))
        r_alpha = parameters[3] / self.scaling_factors[3]
        print("r_alpha = {0:.3e} m  A = {1:.0f} nm^2".format(r_alpha / 100, np.pi * (r_alpha * 1e7)**2))
        print("Sigma = {0:.3e}".format(parameters[4]/self.scaling_factors[4]))
        r_beta = parameters[3] / self.scaling_factors[3] * parameters[4] / self.scaling_factors[4]
        print("r_beta = {0:.3e} m  A = {1:.0f} nm^2".format(r_beta / 100, np.pi * (r_beta * 1e7)**2))

    def plot_mechanism(self, i):
        from matplotlib import pyplot as plt

        rho1, rho2 = self.run_adsorption_model(i)

        plt.figure()
        plt.plot(self.times[i], rho1, label='State 1')
        plt.plot(self.times[i], rho2, label='State 2')
        plt.plot(self.times[i], rho1 + rho2, label='Total')

        plt.axis([0.0, max(self.times[0].max(), self.times[1].max()), 0.0, 160.0])
        plt.xlabel(r"time $(s)$")
        plt.ylabel(r"adsorbed density $(ng/cm^2)$")
        plt.title('VanTassel model: ' + self.plot_title)

        plt.legend(loc='best')
