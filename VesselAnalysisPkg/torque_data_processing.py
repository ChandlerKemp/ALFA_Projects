import numpy as np
from VesselAnalysisPkg import analysis_functions
import csv
import datetime
import matplotlib.pyplot as plt
from .plotting_functions import display_settings
from VesselAnalysisPkg import data_tools


class SeaTrial:
    """
    SeaTrial objects accommodate torque, shaft rpm, tach, power, speed and fuel consumption data from a sea trial
    """
    def __init__(self, file_path, n_eng=1, port=1, plot_eng=1, delimiter='\t'):
        """
        Input params:
            file_in     a raw data file to read
            n_eng       number of engines recorded in this file (must be either 1 or 2)
            port        position of the port engine data in the raw file (1 or 2)
        """
        if n_eng == 1:
            n_headerlines = 11
            self.prop_curve_dict = {'RPM': [], 'power': []}
        elif n_eng == 2:
            n_headerlines = 14
            self.prop_curve_dict = {'engines':{'port':{'RPM':[], 'power': []}, 'starboard':{'RPM': [], 'power': []}}}
        else:
            raise ValueError('n_eng must be 1 or 2')
        if port > n_eng:
            raise ValueError('port must be less than or equal to n_eng')
        self.plot_eng = plot_eng
        self.n_eng = n_eng
        file_var = open(file_path)
        input_data = csv.reader(file_var, delimiter=delimiter)
        n_timesteps = sum(1 for row in input_data) - n_headerlines
        file_var.seek(0)
        input_data = csv.reader(file_var, delimiter=delimiter)
        self.units = "Time: seconds, Torque: ft-lbs, Power: hp, RPM: RPM"
        self.time = np.zeros(n_timesteps)
        self.torque = np.zeros(n_timesteps)
        self.power = np.zeros(n_timesteps)
        self.rpm = np.zeros(n_timesteps)
        self.date_time = []
        if n_eng == 2:
            self.torque2 = np.zeros(n_timesteps)
            self.power2 = np.zeros(n_timesteps)
            self.rpm2 = np.zeros(n_timesteps)
        for row in input_data:
            try:
                if row[0][0:6] == 'Units:':
                    self.units = row[0][7:]
                elif row[0][0:24] == 'TorqueTrak 1 full scale:':
                    self.full_scale_torque = float(row[0][25:-6])
                elif analysis_functions.is_number(row[-1]):
                    data_num = input_data.line_num - n_headerlines - 1
                    self.date_time.append(data_tools.str_to_datetime(row[0]))
                    self.time[data_num] = float(row[1])
                    if port == 1:
                        self.torque[data_num] = float(row[2])
                        self.power[data_num] = float(row[3])
                        self.rpm[data_num] = float(row[4])
                        if n_eng == 2:
                            self.torque2[data_num] = float(row[5])
                            self.power2[data_num] = float(row[6])
                            self.rpm2[data_num] = float(row[7])
                    else:
                        self.torque[data_num] = float(row[5])
                        self.power[data_num] = float(row[6])
                        self.rpm[data_num] = float(row[7])
                        self.torque2[data_num] = float(row[2])
                        self.power2[data_num] = float(row[3])
                        self.rpm2[data_num] = float(row[4])
            except IndexError:
                continue
        file_var.close() 

    def prop_curve(self, inline=True):
        display_settings(inline=inline)
        if self.plot_eng == 1:
            plt_power = self.power
        elif self.plot_eng == 2:
            plt_power = self.power2
        while True:
            plt.plot(self.time, plt_power, label='power')
            plt.title('Power')
            plt.xlabel('Time (s)')
            plt.ylabel('Power (hp)')
            plt.minorticks_on()
            plt.gca().xaxis.grid(b=True, which='major', color='0.65',linestyle='-')
            plt.gca().xaxis.grid(b=True, which='minor', color='r', linestyle='--')
            plt.show()
            xaxis_limits = input('Type the minimum and maximum values for the x axis separated by a comma')
            xaxis_limits = xaxis_limits.split(',')
            a = float(xaxis_limits[0])
            b = float(xaxis_limits[1])
            yaxis_limits = input('Type the minimum and maximum values for the y axis separated by a comma')
            yaxis_limits = yaxis_limits.split(',')
            c = float(yaxis_limits[0])
            d = float(yaxis_limits[1])
            plt.plot(self.time, plt_power, label='power')
            plt.title('Power')
            plt.xlabel('Time (s)')
            plt.ylabel('Power (hp)')
            plt.minorticks_on()
            plt.gca().xaxis.grid(b=True, which='major', color='0.65',linestyle='-')
            plt.gca().xaxis.grid(b=True, which='minor', color='r', linestyle='--')
            plt.xlim([a, b])
            plt.ylim([c, d])
            plt.show()
            plt.figure()
            plt.plot(self.time, self.rpm)
            plt.grid(True)
            plt.title('RPM')
            plt.xlabel('Time (s)')
            plt.ylabel('RPM')
            plt.minorticks_on()
            plt.gca().xaxis.grid(b=True, which='major', color='0.65', linestyle='-')
            plt.gca().xaxis.grid(b=True, which='minor', color='r', linestyle='--')
            plt.xlim([a, b])
            plt.ylim([0, 500])
            plt.show()

            while True:
                print("Type 'quit' to exit")
                start_times = input('Enter the start time for each stage separated by commas (",")')
                if start_times.lower() == 'quit':
                    return
                start_times = start_times.split(',')
                end_times = input('Enter the end time for each stage separated by commas (",")')
                if end_times.lower() == 'quit':
                    return
                end_times = end_times.split(',')
                if len(start_times) == len(end_times):
                    break
                else:
                    print('The number of start times must be equal to the number of end times')
            plt.close('all')
            plt.plot(self.time, plt_power, label='power')
            self.start_times = []
            self.end_times = []
            counter = 0
            delta_t = 0
            while delta_t == 0:
                # look for repeated times
                delta_t = self.time[counter+1] - self.time[counter]
                counter += 1
                if counter > 30:
                    raise ValueError("Input time series is invalid.")
            for ind in range(0, len(start_times)):
                start_times[ind] = float(start_times[ind])
                temp_a = start_times[ind]
                timeind = data_tools.time_match(self.time, temp_a, delta_t)
                self.start_times.append(self.date_time[timeind])
                plt.plot(temp_a, plt_power[timeind], color='r', marker='o')
                end_times[ind] = float(end_times[ind])
                temp_b = end_times[ind]
                timeind = data_tools.time_match(self.time, temp_b, delta_t)
                self.end_times.append(self.date_time[timeind])
                plt.plot(temp_b, float(plt_power[timeind]), color='y', marker='o')
            plt.ylim([c,d])
            plt.xlim([a,b])
            plt.show()
            cont = input("Are the start and end times acceptable? (yes or no)")
            if cont.lower() == 'yes':
                break
            elif cont.lower() == 'quit':
                return
        for ind in range(0, len(start_times)):
            start_ind = data_tools.time_match(self.time, start_times[ind], delta_t)
            end_ind = data_tools.time_match(self.time, end_times[ind], delta_t)
            self.prop_curve_dict['start_times'] = self.start_times
            self.prop_curve_dict['end_times'] = self.end_times
            if self.n_eng == 1:
                self.prop_curve_dict['RPM'].append(np.average(self.rpm[start_ind:end_ind]))
                self.prop_curve_dict['power'].append(np.average(self.power[start_ind:end_ind]))
            elif self.n_eng == 2:
                self.prop_curve_dict['engines']['port']['RPM'].append(np.average(self.rpm[start_ind:end_ind]))
                self.prop_curve_dict['engines']['port']['power'].append(np.average(self.power[start_ind:end_ind]))
                self.prop_curve_dict['engines']['starboard']['RPM'].append(np.average(self.rpm2[start_ind:end_ind]))
                self.prop_curve_dict['engines']['starboard']['power'].append(np.average(self.power2[start_ind:end_ind]))
        self.start_times = start_times
        self.end_times = end_times
