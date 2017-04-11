import numpy as np
from VesselAnalysisPkg import analysis_functions
import csv
import datetime
from VesselAnalysisPkg import plotting_functions
import matplotlib.pyplot as plt
from VesselAnalysisPkg import data_tools

# Constants:
RHO_FUEL = 3179  # estimated fuel density gram/gallon
KW_PER_HP = 0.746  # unit conversion factor
COLORS = [[0,  0, 0], [0.5, 0.5, 0.5], [0, 0, 1], [1, 0, 0], [0, 1, 0], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.]]
class Boat:
    """
        Accumulates all data for a boat
    """
    def __init__(self, file_path, **kwargs):
        """
        Reads data from file_path as well as any keyword arguments
        :param file_path:   Path to a text file containing basic data about the boat. Must have the form label: data
        :param kwargs:      Optional input arguments...generally dicts of data constructed elsewhere.
        """
        self.sea_trials = None
        self.vdr = "No VDR data loaded"
        with open(file_path) as f:
            for row in f:
                temp = row.split(':')
                try:
                    val = float(temp[1].strip())
                except ValueError:
                    val = temp[1].strip()
                setattr(self, temp[0].replace(' ', '_').lower(), val)
            self.hull_speed = 1.34 * np.sqrt(self.length)
            self.bsfc_coeffs_main = None
            self.max_measured_power = None
            self.max_measured_speed = None
            self.min_measured_speed = None
            # Read in torque data if provided
            if 'torque' in kwargs:
                self.sea_trials = kwargs['torque']
            if 'vdr' in kwargs:
                self.vdr = kwargs['vdr']

    def calc_bsfc(self, trials="all", engines = ['main'], checkfit=True):
        """
        Calcuates the brake specific fuel consumption in gal/hp-hr and grams/kWh
        :param trials: a list of trials to be included in the BSFC calculation
        """
        plotting_functions.display_settings(interactive=True, inline=True)
        self.max_measured_power = 0 
        if trials == "all":
            trials = self.sea_trials.keys()
        # prop_curve indicates that the BSFC should be calculated based on manufacturer data, rather than measured data
        ind1 = 0
        if 'bsfc_coeffs' not in dir(self):
            self.bsfc_coeffs = {}
        for engine in engines:
            fuel = []
            power = []
            if trials == "prop_curve":
                fuel.extend(self.prop_curve['fuel'])
                power.extend(self.prop_curve['prop_power'])
            else:
                for trial in trials:
                    # For single propulsion engine vessels, the sea trials dict has no engine keys.
                    # In that case, fuel and power are read from the "trial" level
                    # If there are multiple propulsion engines, then the sea trial dict has more
                    # than one engine key for each trial that must be read separately.
                    if engines == ['main']:
                        temp = self.sea_trials[trial]
                    else:
                        temp = self.sea_trials[trial]['engines'][engine]
                    for ind in range(0,len(temp['fuel'])):
                        if not np.isnan(temp['fuel'][ind]) and not np.isnan(temp['power'][ind]):
                            fuel.append(temp['fuel'][ind])
                            powerind = temp['power'][ind]
                            try:
                                powerind += temp['phantom'][ind]
                            except KeyError:
                                pass
                                
                            power.append(powerind)
                # The bsfc is calculated based on a quadratic fit to the power vs fuel consumption (gal/hr) relationship
                p = data_tools.fit_plotter(power, fuel, deg=2, xlabel='Shaft power (hp)', ylabel = 'Fuel gal/hr', checkfit=checkfit)
                self.bsfc_coeffs[engine] = p
                self.max_measured_power = max(max(power),self.max_measured_power)
                ind1 += 1

    def calc_speed_power(self, trials='all', fit_type = 'default'):
        """
        Calculates a cubic correlation between the shaft power and boat speed
        :param trials: list of sea trials. For example, many boats have "tanked" and "untanked" data
        :param fit_type: Defines the type of fit to use. Can be 'default' or 'constrained'
        :return: None
        """
        plotting_functions.display_settings(interactive=True, inline=True)
        # format trials as a list
        trials = self.select_trials(trials)
        # labels for plots:
        xlab = 'Speed (kn)'
        ylab = 'Power (hp)'
        for trial in trials:
            temp = self.sea_trials[trial]
            if trial.lower() == 'bollard':
                continue
            print('Currently considering the sea trial "' + trial + '"')
            speed = temp['speed']
            # if there are multiple engines, sum their power. Otherwise, simply use 'power'
            if 'engines' in temp:
                power = np.zeros(len(temp['speed']))
                for engine in temp['engines']:
                    power += temp['engines'][engine]['power']
            else:
                power = temp['power']
            self.sea_trials[trial]['speed_power_coeffs'] = data_tools.fit_plotter(speed, power, deg=3, xlabel=xlab, ylabel=ylab, fit_type=fit_type)
            if 'engines' in temp:
                for engine in temp['engines']:
                    power = temp['engines'][engine]['power']
                    temp['engines'][engine]['speed_power_coeffs'] = data_tools.fit_plotter(speed, power, deg=3, xlabel=xlab, ylabel=ylab, fit_type=fit_type)
            if self.max_measured_speed is None:
                self.max_measured_speed = max(speed)
            else:
                self.max_measured_speed = max(self.max_measured_speed, max(speed))
            if self.min_measured_speed is None:
                self.min_measured_speed = min(speed)
            else:
                self.min_measured_speed = min(self.min_measured_speed,min(speed))

    def bsfc_plot(self, engine, save=False, file_out=None, display=True, interactive=True, inline=False, units='gal/hp-hr', mfgdata = True):
        """
        Plot the brake specific fuel consumption curve for the engine
        :param engine: String identifying the engine data to be plotted
        :param save: Boolean determines whether to automatically save the plot
        :param file_out: File to save the figure to if save is True
        :param display: Boolean determines whether to display the plot
        :param inline: If using Jupyter, this can be used to display the plot inline
        :param units: String defining units used for BSFC
        :param mfgdata: include plot of manufacturer data
        :return: None
        """
        plotting_functions.display_settings(interactive=interactive, inline=inline)
        needs_legend = False
        color_ind = 0
        if engine == 'main':
            plist = [self.bsfc_coeffs_main]
        elif engine == 'aux':
            plist = [self.bsfc_coeffs_aux]
        elif type(engine) is list:
            plist = []
            needs_legend = True
            for engine_name in engine:
                plist.append(self.bsfc_coeffs[engine_name])
        else:
            plist = [self.bsfc_coeffs[engine]]
        
        plotting_functions.plot_prep()
        eng_ind = -1
        for p in plist:
            psample = np.linspace(1, self.max_measured_power, 1000)
            bsfc = np.polynomial.polynomial.polyval(psample, np.array(list(reversed(p)))) / psample
            if units.lower() == 'gram/kw-hr':
                bsfc = bsfc * RHO_FUEL / KW_PER_HP
            if needs_legend:
                eng_ind += 1
                label = engine[eng_ind]
            else:
                label = "Measured data"
            plt.plot(psample, bsfc, color=COLORS[color_ind], linewidth=4, label=label)
            color_ind += 1
        if hasattr(self,'prop_curve') and mfgdata:
            needs_legend = True
            power = np.array(self.prop_curve['prop_power'])
            plt.plot(power, np.array(self.prop_curve['fuel'])/power, color=COLORS[color_ind], label="Manufacturer data")
        plt.xlabel('Power (hp)', fontsize=14, fontweight='bold')
        plt.ylabel('BSFC (' + units + ')', fontsize=14, fontweight='bold')
        plt.ylim([0.0, 0.25])
        #try:
        #    plt.plot(self.prop_curve['prop_power'],np.array(self.prop_curve['fuel']) /
        #             np.array(self.prop_curve['prop_power']), color=COLORS[1], label="Manufacturer data")
        #    needs_legend = True
        #except NameError:
        #    pass
        plotting_functions.plot_fixer()
        if needs_legend:
            plt.legend(loc='best')
        plotting_functions.display_save(display, save, file_out)

    def fuel_speed_power_plot(self, trials='all', write_results = True, bsfc_coeffs=None):
        """
        Plot the fuel to speed and power to speed relationships on a plot with two y axes. The correlation for each
        trial is plotted separately.
        :param trials: List of trials to be considered
        :return: None
        """
        f_p = 3  # fuel price $/gal
        plotting_functions.plot_prep()
        fig, ax1 = plt.subplots()
        ax2 = ax1.twinx()
        # format trials as a list
        trials = self.select_trials(trials)
        counter = -1
        if write_results:
            text_file = open("Output.txt", "w")

        for trial in trials:
            counter += 1
            if trial.lower() == 'bollard':
                continue
            temp = self.sea_trials[trial]
            speed = np.linspace(min(temp['speed']), max(temp['speed']), 1000)
            p = list(reversed(temp['speed_power_coeffs']))
            power = np.polynomial.polynomial.polyval(speed, p)
            if 'engines' not in temp:
                neng = 1 # number of engines used in trial
                if bsfc_coeffs is None:
                    fp = list(reversed(self.bsfc_coeffs['main']))
                else:
                    fp = list(reversed(bsfc_coeffs))
                fuel = np.polynomial.polynomial.polyval(power, fp) / speed
            else:
                fuel = np.zeros(speed.shape)
                neng = 0 # counter for number of engines in trial
                # In a trail shaft scenario two engines may be listed, but one should be omitted
                engines = list(temp['engines'].keys())
                for ind in reversed(range(0,len(engines))):
                    if np.isnan(temp['engines'][engines[ind]]['fuel']).all():
                        del engines[ind]
                    else:
                        neng += 1
                for engine in engines:
                    if bsfc_coeffs is None:
                        fp = list(reversed(self.bsfc_coeffs[engine]))
                    else:
                        fp = list(reversed(bsfc_coeffs))
                    fuel += np.polynomial.polynomial.polyval(power/neng, fp) / speed
            if write_results:
                for i in range(0, len(fuel)):
                    if speed[i] - int(speed[i]) <= speed[1]-speed[0]:
                        text_file.write("{:3.0f}".format(speed[i])+ ' & ' + "{:3.1f}".format(power[i]) + ' & ' +
                                        "{:3.2f}".format(fuel[i]) + ' & ' + "{:3.2f}".format(fuel[i] * f_p * 100)
                                        + '\\\\ \n')
            ax1.plot(speed, power, '-', color=COLORS[counter], linewidth=4, label=trial)
            ax2.plot(speed, fuel, '--', color=COLORS[counter], linewidth=4, label=trial)
        ax1.set_xlabel('Speed (kn)', fontsize=14, fontweight='bold')
        ax1.set_ylabel('Power (hp)', fontsize=14, fontweight='bold')
        ax1.grid(linestyle='-', b=True)
        ax2.set_ylabel('Fuel (gal/nm)', fontsize=14, fontweight='bold')
        ax2.grid(linestyle='--', b=True)
        if len(trials) > 1:
            handles, labels = ax1.get_legend_handles_labels()
            ax2.legend(handles, labels, bbox_to_anchor=(0.7, 1.02))
        plotting_functions.plot_fixer(ax1)
        plotting_functions.plot_fixer(ax2)

    def prop_curve_plot(self, trials='all', engines=['main']):
        # Plot manufacturer data if it exists
        plotting_functions.plot_prep()
        # Format trials as a list
        trials = self.select_trials(trials)
        if type(engines) is list:
            for ind in range(0,len(engines)):
                plt.figure(ind)
                try:
                    self.mfg_prop_curve_plt()
                except NameError:
                    pass
        counter = 2
        for trial in trials:
            temp = self.sea_trials[trial]
            if 'engines' not in temp:
                try:
                    self.mfg_prop_curve_plt()
                except NameError:
                    pass
                tach = np.array(self.sea_trials[trial]['RPM']) * self.transmission_ratio
                sortargs = tach.argsort()
                power = np.array(self.sea_trials[trial]['power'])
                plt.plot(tach[sortargs], power[sortargs], label=trial, color=COLORS[counter])
                plt.xlabel('RPM', fontsize=14, fontweight='bold')
                plt.ylabel('Power [hp]', fontsize=14, fontweight='bold')
                plt.legend(bbox_to_anchor=(0.4, 0.95))
                plotting_functions.plot_fixer()
            else:
                ind = -1
                for engine in engines:
                    ind += 1
                    plt.figure(ind)
                    tach = np.array(temp['engines'][engine]['RPM']) * self.transmission_ratio
                    sortargs = tach.argsort()
                    power = np.array(temp['engines'][engine]['power'])
                    plt.plot(tach[sortargs], power[sortargs], label=trial, color=COLORS[counter])
                    counter += 1
        for ind in range(0,len(engines)):
            engine = engines[ind]
            plt.figure(ind)
            plt.xlabel('RPM', fontsize=14, fontweight='bold')
            plt.ylabel('Power [hp]', fontsize=14, fontweight='bold')
            plt.legend(loc='best') #bbox_to_anchor=(0.4, 0.95)
            plt.title(engine)
            plotting_functions.plot_fixer()
            
    def mfg_prop_curve_plt(self):
        plt.plot(self.prop_curve['RPM'], self.prop_curve['prop_power'],
                     label='Mfg prop curve', color=COLORS[0])
        plt.plot(self.prop_curve['RPM'], self.prop_curve['rated_shaft_power'],
                     label="Mfg shaft power", color=COLORS[1])

    def select_trials(self, trials):
        """
        Takes a variety of input types for trials and formats them as a list
        :param trials: string or list identifying sea trials.
        :return: a list of trial identifiers
        """
        if type(trials) is str and trials.lower() == 'all':
            trials = list(self.sea_trials.keys())
        if type(trials) is not list:
            trials = [trials]
        return trials
