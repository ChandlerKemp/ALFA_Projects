import numpy as np
import csv
from VesselAnalysisPkg import data_tools
import datetime

class VDR:
        """
        VDR objects accommodate typical values saved in a vdr
        """
        def __init__(self, file_path, disp=None, n_headerlines=1, **kwargs):
            """
            Initialize a VDR object
            :param file_path: Path to a vdr database exported to a csv
            :param disp: Displacment per revolution (cm^3)
            :param n_headerlines: number of rows above the first row of data.
            :param kwargs: enter an identifying label followed by the index of each column to be read into the object
            :return:
            """

            # The following set of data labels enforces consistent labeling for all VDR objects
            data_labs = ['rpm', 'speed', 'fuel', 'port_fuel', 'star_fuel', 'port_rpm', 'star_rpm', 'port_rpm']
            data_labs.extend(['aux_rpm', 'aux_fuel', 'hyd', 'hyd_1', 'hyd_2', 'date_time', 'date', 'time'])
            data_labs.extend(['supply_temp', 'return_temp', 'port_supply_temp', 'star_supply_temp', 'aux_supply_temp'])
            data_labs.extend(['port_return_temp', 'star_return_temp', 'aux_return_temp'])
            for data_label in kwargs:
                if data_label not in data_labs:
                    raise ValueError('Data label ' + data_label + ' is not a permissible label. Permissible labels are: ' + ', '.join(data_labs))
            with open(file_path) as file_var:
                self.hyd_disp = disp
                # Iterate over each data type specified in kwargs
                # dialect = csv.Sniffer().sniff(file_var.read(2024), delimiters="\t,")
                for data_label in kwargs:
                    file_var.seek(0)
                    input_data = csv.reader(file_var, delimiter = ',')
                    for i in range(0, n_headerlines):
                        next(input_data)
                    temp = []
                    for row in input_data:
                        if row[kwargs[data_label]] == '':
                            temp.append(np.nan)
                            continue
                        if data_label == 'date_time':
                            temp.append(data_tools.str_to_datetime(row[kwargs[data_label]]))
                        else:
                            temp.append(float(row[kwargs[data_label]]))
                    setattr(self, data_label.lower(), np.array(temp))
            if disp is not None:
                self.hyd_power = self.hyd * self.rpm * disp / 60 * 6894.76 * 1e-6 / 746
                
class ACLOGGER:
        def __init__(self, file_path, n_headerlines=1, **kwargs):
            """
            Initialize an AClogger object
            :param file_path: Path to a vdr database exported to a csv
            :param n_headerlines: number of rows above the first row of data.
            :param kwargs: enter an identifying label followed by the index of each column to be read into the object
            :return:
            """
            with open(file_path) as file_var:
                # Iterate over each data type specified in kwargs
                # dialect = csv.Sniffer().sniff(file_var.read(2024), delimiters="\t,")
                for data_label in kwargs:
                    file_var.seek(0)
                    input_data = csv.reader(file_var, delimiter = ',')
                    for i in range(0, n_headerlines):
                        next(input_data)
                    temp = []
                    for row in input_data:
                        if row[kwargs[data_label]] == '':
                            temp.append(np.nan)
                            continue
                        if data_label == 'date_time':
                            temp.append(data_tools.str_to_datetime(row[kwargs[data_label]]))
                        else:
                            temp.append(float(row[kwargs[data_label]]))
                    setattr(self, data_label.lower(), np.array(temp))

class CSV_General:
    """ A general class for storing data from a csv file in a python object"""
    def __init__(self, file_path, n_headerlines=1, **kwargs):
        """
        Initialize a general csv object
        :param file_path: Path to a csv file
        :param n_headerlines: number of rows above the first row of data.
        :param kwargs: enter an identifying label followed by the index of each column to be read into the object
        :return:
        """
        with open(file_path) as file_var:
            # Iterate over each data type specified in kwargs
            # dialect = csv.Sniffer().sniff(file_var.read(2024), delimiters="\t,")
            for data_label in kwargs:
                dl = data_label.lower()
                file_var.seek(0)
                input_data = csv.reader(file_var, delimiter = ',')
                for i in range(0, n_headerlines):
                    next(input_data)
                temp = []
                for row in input_data:
                    if row[kwargs[data_label]] == '':
                        temp.append(np.nan)
                        continue
                    if dl == 'date_time' or dl == 'date' or dl == 'time':
                        temp.append(data_tools.str_to_datetime(row[kwargs[data_label]]))
                    else:
                        temp.append(float(row[kwargs[data_label]]))
                setattr(self, data_label, np.array(temp))

    def date_time_comb(self, dateattr='date', timeattr='time'):
        """
            Combine date and time attributes to form a datetime attributes
            :param dateattr: name of the date attribute (string)
            :param timeattr: name of the time attribute (string)
        """
        date = getattr(self, dateattr)
        time = getattr(self, timeattr)
        setattr(self, 'date_time', [])
        for ind in range(0,len(date)):
             self.date_time.append(datetime.datetime.combine(date[ind], time[ind]))