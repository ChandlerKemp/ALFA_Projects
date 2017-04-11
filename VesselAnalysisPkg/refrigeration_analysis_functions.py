""" This document includes all of the functions called by the refrigeration
	analysis notebook for Eigil B.
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import csv

# Constants:
WATER_HEAT_CAP = 8.35 # btu/gal-F
JOULE_PER_BTU = 1055 
SEC_PER_HR = 3600
WATT_PER_KW = 1000
KW_PER_HP = 0.746

class Vessel_RSW:
	""" An object for storing all vessel specific information that is 
		independent of the RSW system """
	def __init__(self):
		self.sea_temp = None # deg F
		self.source_temp = None # deg F
		self.target_temp = None
		self.tank_capacity = None # gallons
		self.n_tank_pulldown = None # pulldowns per year
		self.lbs_fish_cooled = None #lbs
		self.hrs_holding_temp = None # hrs
		self.ref_hours = None # hrs
		self.condenser_pump_cap = None # kW
		self.compressor_cool_pump = None # kW
		# rate at which the holds warm when at temperature (deg F/hr)
		self.heat_infiltration_33 = None
		self.efficiency_penalty = None

class Compressor:
	""" An object for storing all specifications that characterize a Compressor 
		system"""
	def __init__(self, sst, capacity, power, aux_load, hold_temps,\
				n_comp=1, switch=False, unloaders=False):
		""" Initialization function for the class
			:param sst: a list of saturated suction temperatures (deg F)
			:param capacity: a list of compressor capacities at the given S.S.
				temps (btu/hr)
			:param power: a list of compressor input powers at the given S.S. 
				temps (kW)
			:param aux_load: constant auxiliary loads to the compressor (kW)
			:param n_comp: number of compressors
			:param switch: Boolean to determine ability to run just one 
				compressor
			:param unloaders: If the compressor has an unloader, unloaders
				should be set equal to a dict with two entries: 
				'power': ratio of power with unloader to power without
				'capacity': ratio of capacity with unloader to capacity without
				Both 'power' and 'capacity' may be floats or lists or arrays
		"""
		self.sst_data = np.array(sst)
		self.capacity_data = np.array(capacity)
		self.power_data = np.array(power)
		self.power = interp1d(sst, power, kind='linear')
		self.capacity = interp1d(sst, capacity, kind='linear')
		self.sst = interp1d(hold_temps, sst, kind='linear')
		self.aux_load = aux_load
		factor = JOULE_PER_BTU / SEC_PER_HR / WATT_PER_KW
		self.cop = self.capacity_data/ self.power_data * factor
		self.n = n_comp
		self.switch = switch
		# store unloaders['power'] and ['capacity'] as arrays if they are
		# entered as scalars
		if not unloaders:
			self.unloaders = {'power':np.array([1]), 'capacity':np.array([1])}
		else:
			try:
				len(unloaders['power'])
				self.unloaders = unloaders
			except TypeError:
				# settings are always stored from greatest to least
				unloaders['power'] = -np.sort(-np.array(unloaders['power']))
				unloaders['capacity'] = -np.sort(-np.array(\
					unloaders['capacity']))
				self.unloaders = unloaders
		# index of unloader setting
		self.unload_ind = 0

	def print_cop(self):
		strout = 'The COP is '
		for cop in self.cop:
			strout += '{:1.2f}'.format(cop) + ', '
		print_string(strout[:-2])
		strout = "For Saturated Suction Temps (deg F) of "
		for sst in self.sst_data:
			strout += '{:1.1f}'.format(sst) + ', '
		print_string(strout[:-2])

	def print_capacity(self):
		strout = "The capacity (kBtu/hr) is "
		for cap in self.capacity_data:
			strout += '{:3.0f}'.format(cap/1000 * self.n) + ', '
		print_string(strout[:-2])

class Evaporator:
	""" An object for storing all specifications that characterize an evaporator
	"""
	def __init__(self, temp_delta):
		""" Initialization function
			:param evap_temp_delta: the maximum allowable water temperature rise
				throught the evaporator.
		"""
		self.temp_delta = temp_delta

class CircPump:
	""" An object for storing all specification that characterize a circulation
		pump"""
	def __init__(self, capacity, power, vfd=False, n_pump=1, switch=False, \
		fr=1, motor_eff=1):
		""" Initizialization function
			:param capacity: Flow through the pump when powered at 60 Hz 
				(gal/hr)
			:param power: rated power for the pump (kW)
			:param vfd: boolean to enable variable frequency operation
			:param n_pump: integer number of pumps
			:param switch: boolean to enable operation of just one pump
			:param fr: ratio of pump speed while filling tanks to rated
				pump speed
			:param motor_eff: ratio of power input to motor over power consumed
				by pump
		"""
		self.capacity = capacity # gal/hr
		# the power stored is the power input to the motor, not the rated power
		# of the pump
		self.power = power / motor_eff # kW
		self.vfd = vfd
		self.n_pumps = n_pump
		self.switch = switch
		self.fill_ratio = fr
		self.motor_eff = motor_eff


class Generator:
	""" An object for storing all specification that characterize a generator"""
	def __init__(self, eng_power, fuel, motor_eff):
		""" Initizialization function based on engine mfg. specifications. 
			motor_eff defines the fraction of power produced by the engine that
			is converted to electricity.
			:param eng_power: power produced by the engine
			:param fuel: fuel consumed by the engine at the powers given by
				eng_power
			:param motor_eff: the fraction of engine power converted to
				electricity by the generator.
		"""
		self.eng_power = eng_power
		self.fuel = fuel
		self.motor_eff = motor_eff
		power_kw = eng_power * motor_eff
		self.bsfc_coeffs = np.polyfit(power_kw, fuel, 1)

	def fuel_rate(self, power, marginal=False):
		""" a function for calculating the rate of fuel consumption for a load.
			:param power: load on the generator (kW)
			:param marginal: If true, it's assumed that the generator will run
				regardless of whether or not the refrigeration system is on.
			:return fuel: the rate of fuel consumption (gal/hr)
		"""
		fuel_a = self.bsfc_coeffs[0]
		fuel_b = self.bsfc_coeffs[1]
		if marginal:
			fuel_b = 0
		fuel = (fuel_a*power + fuel_b)
		return fuel

	def plot_fit(self):
		""" function for plotting the BSFC fit to the input data"""
		powertest = np.linspace(0, max(self.eng_power) * self.motor_eff, 100)
		fueltest = self.fuel_rate(powertest)
		plt.plot(self.eng_power * self.motor_eff, self.fuel, 'o',\
			label='Mfg data')
		plt.plot(powertest, fueltest, label='fit')
		plt.ylabel('Fuel (gal/hr)')
		plt.xlabel('Engine Power (kW)')
		plt.legend(loc='best')


class Season:
	""" An object for storing all data from a season simulation """
	def __init__(self, ves, comp, circ, evap, gen, name=None,\
			override=False, marginal=False):
		self.fuel = {}
		self.fuel['filling tanks'] = fuel_filling_tanks(ves, circ, gen,\
			marginal=marginal)
		vol = ves.tank_capacity
		time, fuel = fuel_pulldown(vol, ves, comp, circ, evap, gen, \
						marginal=marginal)
		self.fuel_per_load = fuel
		self.time_per_load = time
		self.fuel['pulldowns']= fuel * ves.n_tank_pulldown
		if ves.hrs_holding_temp is None and ves.ref_hours is not None:
			hrs = max(0, ves.ref_hours - time * ves.n_tank_pulldown)
		else:
			hrs = ves.hrs_holding_temp
		dutycycle, fuel = fuel_maintain(vol, hrs, ves, comp, circ, evap, gen,\
			override=override, marginal=marginal)
		self.dutycycle = dutycycle
		self.fuel['steady temp'] = fuel
		self.fuel_total = 0
		for fuel in self.fuel.values():
			self.fuel_total += fuel
		self.name = name

def unloaders(volume, heat_in_rate, ves, comp, circ, evap, maint=True):
	"""
	This function calculates the index of the unload condition.

	:param volume: volume of the holds being cooled (gal)
	:param heat_in_rate: Rate at which heat enters the hold from ambient 
		(btu/hr/gal). An array with the same length as comp.sst
	:param ves: vessel object
	:param comp: a compressor object
	:param circ: A circulation pump object
	:param evap: An evaporator object
	:param maint: boolean to determine whether the hold is being maintained at
		at temp or pumped down
	:return : total heating load (btu/hr)
	"""
	if not maint:
		return 0
	else:
		for ind in reversed(range(0,len(comp.unloaders['power']))):
			comp.unload_ind = ind
			try:
				heat_load_fcn(volume, heat_in_rate, ves, comp, circ, evap, \
					ves.target_temp, maint=maint)
				return ind
			except ValueError:
				pass
		# if the function gets to here, there's a problem
		raise ValueError('Insufficient capacity to maintain desired hold0 \
			temperature.')

def heat_load_fcn(volume, heat_in_rate, ves, comp, circ, evap, hold_temp, \
					maint=False):
	"""
	This function calculates the total heating load on the 
	system at each sst given in comp. This function allows for VFDs on the circ
	pump, and for one circ pump to be switched off.

	:param volume: volume of the holds being cooled (gal)
	:param heat_in_rate: Rate at which heat enters the hold from ambient 
		(btu/hr/gal). An array with the same length as comp.sst
	:param ves: vessel object
	:param comp: a compressor object
	:param circ: A circulation pump object
	:param evap: An evaporator object
	:param maint: boolean to determine whether the hold is being maintained at
		at temp or pumped down
	:return : total heating load (btu/hr)
	"""
	cond = type(heat_in_rate) is np.ndarray or type(heat_in_rate) is list
	if cond and maint:
		raise ValueError("heat_in_rate must be a scalar if maint is True")

	amb_load = volume * heat_in_rate # btu/hr
	if maint:
		# error message raised if capacity is not sufficient to maintain hold
		# temp.
		sst = comp.sst(ves.target_temp)
		err_msg = 'Insufficient capacity to maintain desired hold temperature.'
		min_flow = comp.capacity(sst)/ WATER_HEAT_CAP / evap.temp_delta
	    # Allows one compressor to be turned off if one is enough to maintain 
		# hold temp
		if comp.switch is True:
			n_comp = 1
		else:
			n_comp = comp.n
		min_flow = min_flow * n_comp
		circ_pow_ratio = circ_power_ratio(min_flow, circ)
		circ_load = circ.power * circ_pow_ratio * circ.n_pumps * \
					circ.motor_eff
		circ_load = circ_load / JOULE_PER_BTU * WATT_PER_KW * SEC_PER_HR
		if amb_load <= comp.capacity(sst) - circ_load:
			return amb_load + circ_load
		elif comp.switch is True:
			min_flow = comp.capacity(sst) / WATER_HEAT_CAP / evap.temp_delta \
				* comp.n
			circ_pow_ratio = circ_power_ratio(min_flow, circ)
			circ_load = circ_pow_ratio * circ.power * circ.n_pumps \
				* circ.motor_eff
			circ_load = circ_load / JOULE_PER_BTU * WATT_PER_KW * SEC_PER_HR
			if amb_load <= comp.capacity(sst) * comp.n - circ_load:
				return amb_load + circ_load
			else:
				raise ValueError(err_msg)
		else:
			raise ValueError(err_msg)
	else:
		# the compressor runs at full load unless the hold is at temp.
		comp.unload_ind = 0
		min_flow = min_flow_fcn(comp, evap, hold_temp)
		circ_load = circ.power * circ.n_pumps * circ_power_ratio(min_flow, circ)
		circ_load = circ_load / JOULE_PER_BTU * WATT_PER_KW * SEC_PER_HR
		return circ_load + amb_load


def fuel_filling_tanks(ves, circ, gen, marginal=False):
	""" A function for calculating the fuel used to fill tanks per season
	:param circ: a circ pump object
	:param gen: a generator object
	:param ves: a vessel object
	:return fuel: fuel used to fill tanks per season
	"""
	load_kw = circ.n_pumps * circ.power
	volume = ves.tank_capacity * ves.n_tank_pulldown # gal
	time = volume / circ.capacity * 2 # hrs factor of 2 allows for emptying hold
	if circ.vfd:
		load_kw *= circ.fill_ratio**3
		time /= circ.fill_ratio
	fuel = gen.fuel_rate(load_kw, marginal=marginal) * time
	return fuel

def fuel_pulldown(water_vol, ves, comp, circ, evap, gen, marginal=False):
	""" A function for calculating fuel consumption of an RSW 
		system during pull down
	:param water_vol: volume of water to be cooled (gallons)
	:param ves: a vessel object
	:param comp: a compressor object
	:param circ: a circulation pump object
	:param evap: an evaporator object
	:param marginal: If true, only the marginal fuel consumption is considered
		(idle fuel consumption is not)
	:return time: time to pull down hold (hrs)
	:return fuel: fuel burned during pull down (gallon)
	"""
	ntemps = 100
	temps = np.linspace(ves.source_temp, ves.target_temp, ntemps)
	avtemps = (temps[:-1]+temps[1:])/2
	sst = comp.sst(avtemps)
	btus_removed = water_vol * WATER_HEAT_CAP * (temps[:-1]-temps[1:])
	source_tratio = (ves.sea_temp - ves.source_temp)/(ves.sea_temp - \
		ves.target_temp)
	heat_in_rate = np.linspace(source_tratio, ves.heat_infiltration_33, \
			ntemps-1) * WATER_HEAT_CAP #btu/hr/gal
	heat_load = heat_load_fcn(water_vol, heat_in_rate, ves, comp, circ, evap, \
				avtemps)
	times = btus_removed / (comp.capacity(sst) * comp.n - heat_load)
	min_flow = min_flow_fcn(comp, evap, avtemps)
	cpr = circ_power_ratio(min_flow, circ)
	other_loads = comp.aux_load * comp.n + circ.power * circ.n_pumps * cpr
	total_load = comp.power(sst) * comp.n+other_loads
	fuel = sum(gen.fuel_rate(total_load, marginal=marginal)*times)
	time = sum(times)
	return time, fuel

def fuel_maintain(vol, hrs, ves, comp, circ, evap, gen, override=False,\
			marginal=False):
	""" Function for calculating fuel consumed to maintain hold temp 
		:param vol: volume of hold maintained
		:param hrs: hours to maintain hold
		:param ves: vessel object
		:param comp: compressor object
		:param circ: circulation pump object
		:param evap: evaporator object
		:param gen: generator object
		:return fuel: fuel consumed to maintain hold temperature
	"""
	temp = ves.target_temp
	sst = comp.sst(temp)
	heat_in = vol * ves.heat_infiltration_33 #btu/hr
	hr = ves.heat_infiltration_33
	# Set the unload condition
	comp.unload_ind = unloaders(vol, hr, ves, comp, circ, evap)
	# power requirement [kW] and capacity [btu/hr] of the compressor
	comp_pow = comp.power(sst) * comp.unloaders['power'][comp.unload_ind]
	cap = comp.capacity(sst) * comp.unloaders['power'][comp.unload_ind]
	hl_circ= heat_load_fcn(vol, hr, ves, comp, circ, evap, temp, maint=True)
	# Note: circ pumps are assumed to be always running, whether or not the
	# compressor is on.
	if comp.switch and cap-hl_circ > 0:
		# The dutycycle is equal to the ratio of the rate that heat infiltrates
		# the system to the rate that heat is removed. Heat is removed at a 
		# rate of comp.capacity[-1] minus the heat load of the circ pumps.
		# the heat load of the circ pumps is hl_circ-heat_in
		min_flow = cap / WATER_HEAT_CAP / evap.temp_delta
		cpr = circ_power_ratio(min_flow, circ)
		circpow = circ.power * circ.n_pumps * cpr
		power = comp_pow + circpow + comp.aux_load
		duty_cycle = (heat_in + hl_circ)/ cap
	else:
		min_flow = cap / WATER_HEAT_CAP / evap.temp_delta \
			* comp.n
		cpr = circ_power_ratio(min_flow, circ)
		circpow = circ.power * circ.n_pumps * cpr
		power = comp_pow * comp.n + circpow + comp.aux_load * comp.n
		duty_cycle = (heat_in + hl_circ)/ cap
	if override is not False:
		duty_cycle = override
	comp.unload_ind = 0
	fuel = gen.fuel_rate(power, marginal=marginal) * hrs * duty_cycle
	return duty_cycle, fuel


def circ_power_ratio(min_flow, circ):
	""" Function for calculating the ratio of power used to rated power of a 
		circ pump.
		:param min_flow: minimum allowable flow through the evaporator.
			can be an array
		:param circ: a circ pump object
	"""
	if circ.vfd is True:
		circ_pow_ratio = (min_flow / (circ.capacity * circ.n_pumps))**3
	elif circ.n_pumps == 2 and circ.switch is True:
		circ_pow_ratio = []
		try:
			for f in min_flow:
				if f < circ.capacity:
					circ_pow_ratio.append(0.5)
				else:
					circ_pow_ratio.append(1)
			circ_pow_ratio = np.array(circ_pow_ratio)
		except TypeError:
			if min_flow < circ.capacity:
				circ_pow_ratio = 0.5
			else:
				circ_pow_ratio = 1
	else:
		circ_pow_ratio = 1
	return circ_pow_ratio

def min_flow_fcn(comp, evap, hold_temp):
	""" A function for calculating the minimum required flow through an
		evaporator
		:param comp: a compressor object
		:param evap: an evaporator object
		:return min_flow: an array of minimum flows (gal/hr) 
	"""
	sst = comp.sst(hold_temp)
	min_flow = comp.capacity(sst) * comp.n / WATER_HEAT_CAP / evap.temp_delta
	return min_flow		

def print_string(string):
	""" Prints a string with wrapping to avoid run-off lines when exporting to
		pdf.

		:param string: a string to be printed
	"""
	linelength = 79 # characters
	str1 = ''
	while len(string) > linelength:
		str1 += string[0:linelength+1] + '\n'
		string = string[linelength+1:]
	str1 += string

	print(str1)

def results_table(seasons):
	""" prints a results table """

	print("Total fuel consumed in a year (gallons):")
	for season in seasons:
		print(season.name, '{:4.0f}'.format(season.fuel_total))


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