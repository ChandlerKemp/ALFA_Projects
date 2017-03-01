""" This document includes all of the functions called by the refrigeration
	analysis notebook for Eigil B.
"""

import numpy as np
import matplotlib.pyplot as plt

# Constants:
WATER_HEAT_CAP = 8.35 # btu/gal-F
JOULE_PER_BTU = 1055 
SEC_PER_HR = 3600
WATT_PER_KW = 1000
KW_PER_HP = 0.746

class Vessel_RSW:
	""" An object for storing all vessel specific information that is 
		indpendent of the RSW system """
	def __init__(self):
		self.sea_temp = None
		self.target_temp = None
		self.aft_tank_capacity = None # gallons
		self.forward_tank_capacity = None # gallons
		self.n_aft_tank_pulldown = None # pulldowns per year
		self.n_all_tank_pulldown = None # pulldowns per year
		self.lbs_fish_cooled = None #lbs
		self.hrs_holding_aft_temp = None # hrs
		self.hrs_holding_all_temp = None # Herring season only--no data available
		self.condenser_pump_cap = None # kW
		self.compressor_cool_pump = None # kW
		# rate at which the holds warm when at temperature (deg F/hr)
		self.heat_infiltration_33 = None
		self.efficiency_penalty = None

class Compressor:
	""" An object for storing all specifications that characterize a Compressor 
		system"""
	def __init__(self, sst, capacity, power, aux_load, n_comp=1, switch=False):
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
		"""
		self.sst = np.array(sst)
		self.capacity = np.array(capacity)
		self.power = np.array(power)
		self.aux_load = aux_load
		factor = JOULE_PER_BTU / SEC_PER_HR / WATT_PER_KW
		self.cop = self.capacity / self.power * factor
		self.n = n_comp
		self.switch = switch

	def print_cop(self):
		strout = 'The COP is '
		for cop in self.cop:
			strout += '{:1.2f}'.format(cop) + ', '
		print_string(strout[:-2])
		strout = "For Saturated Suction Temps (deg F) of "
		for sst in self.sst:
			strout += '{:1.1f}'.format(sst) + ', '
		print_string(strout[:-2])

	def print_capacity(self):
		strout = "The capacity (kBtu/hr) is "
		for cap in self.capacity:
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
	def __init__(self, capacity, power, vfd=False, n_pump=1, switch=True, fr=1):
		""" Initizialization function
			:param capacity: Flow through the pump when powered at 60 Hz 
				(gal/hr)
			:param power: rated power for the pump (kW)
			:param vfd: boolean to enable variable frequency operation
			:param n_pump: integer number of pumps
			:param switch: boolean to enable operation of just one pump
			:param fr: ratio of pump speed while filling tanks to rated
				pump speed
		"""
		self.capacity = capacity # gal/hr
		self.power = power # kW
		self.vfd = vfd
		self.n_pumps = n_pump
		self.switch = switch
		self.fill_ratio = fr


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

	def fuel_rate(self, power):
		""" a function for calculating the rate of fuel consumption for a load.
			:param power: load on the generator (kW)
			:return fuel: the rate of fuel consumption (gal/hr)
		"""
		fuel_a = self.bsfc_coeffs[0]
		fuel_b = self.bsfc_coeffs[1]
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
	def __init__(self, ves, comp, circ, evap, gen, name=None, override=False):
		self.fuel = {}
		self.fuel['filling tanks'] = fuel_filling_tanks(ves, circ, gen)
		vol = ves.aft_tank_capacity
		time, fuel = fuel_pulldown(vol, ves, comp, circ, evap, gen)
		self.salmon_fuel_per_load = fuel
		self.salmon_time_per_load = time
		self.fuel['salmon pulldowns']= fuel * ves.n_aft_tank_pulldown
		hrs = ves.hrs_holding_aft_temp
		dutycycle, fuel = fuel_maintain(vol, hrs, ves, comp, circ, evap, gen,\
			override=override)
		self.salmon_dutycycle = dutycycle
		self.fuel['salmon maintenance'] = fuel
		vol = ves.aft_tank_capacity + ves.forward_tank_capacity
		time, fuel = fuel_pulldown(vol, ves, comp, circ, evap, gen)
		self.herring_fuel_per_load = fuel
		self.herring_time_per_load = time
		self.fuel['herring pulldowns']= fuel * ves.n_all_tank_pulldown
		hrs = ves.hrs_holding_all_temp
		dutycycle, fuel = fuel_maintain(vol, hrs, ves, comp, circ, \
			evap, gen, override=override)
		self.herring_dutycycle = dutycycle
		self.fuel['herring maintenance'] = fuel
		self.fuel_total = 0
		for fuel in self.fuel.values():
			self.fuel_total += fuel
		self.name = name


def heat_load_fcn(volume, heat_in_rate, ves, comp, circ, evap, maint=False):
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
		min_flow = comp.capacity[-1] / WATER_HEAT_CAP / evap.temp_delta
	 # Allows one compressor to be turned off if one is enough to maintain hold
		# temp
		if comp.switch is True:
			n_comp = 1
		else:
			n_comp = comp.n
		min_flow = min_flow * n_comp
		circ_pow_ratio = circ_power_ratio(min_flow, circ)
		circ_load = circ.power * circ_pow_ratio * circ.n_pumps
		circ_load = circ_load / JOULE_PER_BTU * WATT_PER_KW * SEC_PER_HR
		if amb_load <= comp.capacity[-1] - circ_load:
			return amb_load + circ_load
		else:
			min_flow = comp.capacity[-1] / WATER_HEAT_CAP / evap.temp_delta \
				* comp.n
			circ_pow_ratio = circ_power_ratio(min_flow, circ)
			circ_load = circ_pow_ratio * circ.power * circ.n_pumps
			circ_load = circ_load / JOULE_PER_BTU * WATT_PER_KW * SEC_PER_HR
			return amb_load + circ_load
	else:
		min_flow = min_flow_fcn(comp, evap)
		circ_load = circ.power * circ.n_pumps * circ_power_ratio(min_flow, circ)
		circ_load = circ_load / JOULE_PER_BTU * WATT_PER_KW * SEC_PER_HR
		return circ_load + amb_load


def fuel_filling_tanks(ves, circ, gen):
	""" A function for calculating the fuel used to fill tanks per season
	:param circ: a circ pump object
	:param gen: a generator object
	:param ves: a vessel object
	:return fuel: fuel used to fill tanks per season
	"""
	load_kw = circ.n_pumps * circ.power
	volume1 = ves.aft_tank_capacity * ves.n_aft_tank_pulldown # gal
	volume2 = ves.aft_tank_capacity + ves.forward_tank_capacity
	volume = volume1 + volume2 * ves.n_all_tank_pulldown
	time = volume / circ.capacity * 2 # hrs factor of 2 allows for emptying hold
	if circ.vfd:
		load_kw *= circ.fill_ratio**3
		time /= circ.fill_ratio
	fuel = gen.fuel_rate(load_kw) * time
	return fuel

def fuel_pulldown(water_vol, ves, comp, circ, evap, gen):
	""" A function for calculating fuel consumption of an RSW 
		system during pull down
	:param water_vol: volume of water to be cooled (gallons)
	:param ves: a vessel object
	:param comp: a compressor object
	:param circ: a circulation pump object
	:param evap: an evaporator object
	:return time: time to pull down hold (hrs)
	:return fuel: fuel burned during pull down (gallon)
	"""
	btus_removed = water_vol * WATER_HEAT_CAP * (ves.sea_temp - ves.target_temp)
	if len(comp.sst) > 1:
		delta_sst = comp.sst[0] - comp.sst[1]
		btu_per_sst = btus_removed/(delta_sst * len(comp.sst))
	else:
		btu_per_sst = btus_removed
		delta_sst = 1
	heat_in_rate = np.linspace(0, ves.heat_infiltration_33, \
			len(comp.sst)) * WATER_HEAT_CAP #btu/hr/gal
	heat_load = heat_load_fcn(water_vol, heat_in_rate, ves, comp, circ, evap)
	times = btu_per_sst * delta_sst / (comp.capacity * comp.n - heat_load)
	min_flow = min_flow_fcn(comp, evap)
	cpr = circ_power_ratio(min_flow, circ)
	other_loads = comp.aux_load * comp.n + circ.power * circ.n_pumps * cpr
	fuel = sum(gen.fuel_rate(comp.power * comp.n+other_loads)*times)
	time = sum(times)
	return time, fuel

def fuel_maintain(vol, hrs, ves, comp, circ, evap, gen, override=False):
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
	heat_in = vol * ves.heat_infiltration_33 #btu/hr
	hr = ves.heat_infiltration_33
	hl_circ= heat_load_fcn(vol, hr, ves, comp, circ, evap, maint=True)
	if comp.switch and comp.capacity[-1]-hl_circ > 0:
		# The dutycycle is equal to the ratio of the rate that heat infiltrates
		# the system to the rate that heat is removed. Heat is removed at a 
		# rate of comp.capacity[-1] minus the heat load of the circ pumps.
		# the heat load of the circ pumps is hl_circ-heat_in
		min_flow = comp.capacity[-1] / WATER_HEAT_CAP / evap.temp_delta
		cpr = circ_power_ratio(min_flow, circ)
		circpow = circ.power * circ.n_pumps * cpr
		power = comp.power[-1] + circpow + comp.aux_load
		duty_cycle = heat_in / (comp.capacity[-1] - (hl_circ-heat_in))
	else:
		min_flow = comp.capacity[-1]/ WATER_HEAT_CAP / evap.temp_delta * comp.n
		cpr = circ_power_ratio(min_flow, circ)
		circpow = circ.power * circ.n_pumps * cpr
		power = comp.power[-1] * comp.n + circpow + comp.aux_load * comp.n
		duty_cycle = heat_in / (comp.capacity[-1] * comp.n - (hl_circ-heat_in))
	if override is not False:
		duty_cycle = override
	return duty_cycle, gen.fuel_rate(power) * hrs * duty_cycle


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

def min_flow_fcn(comp, evap):
	""" A function for calculating the minimum required flow through an
		evaporator
		:param comp: a compressor object
		:param evap: an evaporator object
		:return min_flow: an array of minimum flows (gal/hr) 
	"""
	min_flow = comp.capacity * comp.n / WATER_HEAT_CAP / evap.temp_delta
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