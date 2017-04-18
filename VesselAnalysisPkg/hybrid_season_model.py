import numpy as np
from scipy.interpolate import interp1d
from numpy.polynomial.polynomial import polyval
import copy

FUEL_DENSITY = 3179 # gram/gallon
KW_PER_HP = 0.746
PASCAL_PER_PSI = 6984


def hybrid_fuel_savings(load_fractions, load_sizes, main_bsfc, aux_bsfc, aux_cap, total_hrs):
    """
        Calculates total fuel savings with a hybrid propulsion system
        :param load_fractions: array of fractional time in each load size
        :param load_sizes: array of loads associated with each fraction
        :param main_bsfc: bsfc parameters for the main engine
        :param aux_bsfc: bsfc parameters for the auxiliary power source
        :param aux_cap: capacity of the auxiliary power source.
        :param total_hrs: total number of hours to be modeled
        :return main_fuel: total fuel burned with only the main engine
        :return hybrid_fuel: total fuel burned with the hybrid system
    """
    main_fuel = sum(polyval(load_sizes, main_bsfc)*load_fractions) * total_hrs
    cond = (load_sizes > aux_cap)
    hybrid_fuel = sum(polyval(load_sizes[cond], main_bsfc) * load_fractions[cond])
    cond = np.invert(cond)
    hybrid_fuel += sum(polyval(load_sizes[cond], aux_bsfc) * load_fractions[cond])
    hybrid_fuel *= total_hrs
    return main_fuel, hybrid_fuel


def hybrid_fuel_savings(load_fractions, load_sizes, main_bsfc, aux_bsfc, aux_cap, total_hrs):
    """
        Calculates total fuel savings with a hybrid propulsion system
        :param load_fractions: array of fractional time in each load size
        :param load_sizes: array of loads associated with each fraction
        :param main_bsfc: bsfc parameters for the main engine
        :param aux_bsfc: bsfc parameters for the auxiliary power source
        :param aux_cap: capacity of the auxiliary power source.
        :param total_hrs: total number of hours to be modeled
        :return main_fuel: total fuel burned with only the main engine
        :return hybrid_fuel: total fuel burned with the hybrid system
    """
    main_fuel = sum(polyval(load_sizes, main_bsfc)*load_fractions) * total_hrs
    cond = (load_sizes > aux_cap)
    hybrid_fuel = sum(polyval(load_sizes[cond], main_bsfc) * load_fractions[cond])
    cond = np.invert(cond)
    hybrid_fuel += sum(polyval(load_sizes[cond], aux_bsfc) * load_fractions[cond])
    hybrid_fuel *= total_hrs
    return main_fuel, hybrid_fuel


def power_of_speed(sp, sp_pow_coeffs, speed_data, power_data):
    """
        Calculates the propulsion power based on vessel speed
        For speeds between min(speed_data) and max(speed_data),
        a cubic spline is used to interpolate between points
        For speeds less than min(speed_data), power is assumed to
        be proportional to the speed cubed.
        For speeds greater than max(speed_data), speed is assumed
        to follow the curve determined by sp_pow_coeffs
        :param sp: speed at which to calculate required power (kn)
        :param sp_pow_coeffs: speed to power coefficients to be \
            in upper extrapolation
        :param speed_data: list or array of recorded speeds
        :param power_data: list or array of recorded powers
        :return power: power required to achieve desired speeds (hp)
    """
    if type(sp) is int or type(sp) is float:
        sp = [sp]
    if type(sp) is list:
        sp = np.array(sp)
    power = np.zeros(sp.shape)
    p = sp_pow_coeffs
    x = np.array(speed_data)
    y = np.array(power_data)
    power_fcn_internal = interp1d(x, y, kind='cubic')
    a = power_fcn_internal(x[0])/ x[0]**3
    cutoff = x[0]
    upper = x[-1]
    power[sp <= cutoff] = a * sp[sp <= cutoff]**3
    power[ (upper > sp) & (sp > cutoff)] = power_fcn_internal(sp[(upper > sp) & (sp > cutoff)])
    power[sp >= upper] = polyval(sp[sp >= upper],p)
    return power


def power_of_fuel(fuel, bsfc_coeffs, pmax=150, track_err = True):
    """
        Returns the power based on measured fuel consumption
        Uses an iterative method to calculate implied power
        based on a polynomial equation for fuel consumption
        (gal/hr) as a function of power.
        :param fuel: rate of fuel consumption (gal/hr)
        :param bsfc_coeffs: bsfc coefficients (gal, hp, hr)
        :return power: power implied by fuel consumption (hp)
    """
    p = np.array(list(reversed(bsfc_coeffs)))
    power = []
    low_fuel_counter = 0
    for f in fuel:
        if f > p[0]:
            phigh = pmax
            plow = 0.01
            err = 100
            counter = 0
            while err > 0.1 and counter < 100:
                counter+=1
                if counter==100 and track_err:
                    print("counter limit reached \n",\
                      "fuel is: ", f, "\n High power is: ", phigh,\
                      "\n Low power is: ", plow)
                ptemp = (phigh + plow)/2
                fuelguess = polyval(ptemp,p)
                if fuelguess > f:
                    phigh = copy.copy(ptemp)
                else:
                    plow = copy.copy(ptemp)
                err = np.abs(f-fuelguess)
            power.append(ptemp)
        else:
            low_fuel_counter += 1
            power.append(0)
    if low_fuel_counter > 0:
        print("Warning: fuel flow less than idle flow was recorded at  " + str(low_fuel_counter) + " timesteps.")
    return np.array(power)


def power_of_tach(tach, tach_power_coeffs, pmax=150, track_err = True, speed=None):
    """
        Calculates the propulsion power based on engine rpm
    :param tach: recording of propulsion engine speed (rpm)
    :param tach_power_coeffs: polynomial coefficients relating
        engine speed to propulsion power (rpm, hp)
    :param pmax: maximum possible power (hp)
    :return power: propulsion power (hp)
    """
    try:
        a = len(tach)
    except TypeError:
        tach = np.array([tach])
        speed = [speed]
    p = np.array(list(reversed(tach_power_coeffs)))
    power = []
    for i in range(0,len(tach)):
        if speed is None or speed[i] > 1:
            t = tach[i]
            phigh = pmax
            plow = 0
            err = 100
            counter = 0
            while err > 0.1 and counter < 100:
                counter+=1
                if counter==100 and track_err:
                    print("counter limit reached \n",\
                          "tach is: ", t, "\n High power is: ", phigh,\
                          "\n Low power is: ", plow)
                ptemp = (phigh + plow)/2
                tachguess = polyval(ptemp,p)
                if tachguess > t:
                    phigh = copy.copy(ptemp)
                else:
                    plow = copy.copy(ptemp)
                err = np.abs(t-tachguess)
            power.append(ptemp)
        else:
            power.append(0)
    power = np.array(power)
    return power # hp


def hyd_pow(press, hz, d0, pmax=3000*PASCAL_PER_PSI,\
            pmin=3000*PASCAL_PER_PSI):
    """
        Calculates the power consumed by a hydraulic pump
        Can be pressure compensating
    :param press: hydraulic pressure (pascal)
    :param hz: rate of pump rotation (hz)
    :param d0: displacement per revolution at 0 hz (m^3)
    :param pmax: Maximum pump pressure=pressure at which \
        displacement is zero
    :param pmin: pressure at which the swash plate begins \
        to change angle in a pressure compensating pump.
        To model a fixed displacement pump, let pmin=pmax
    :return power: power consumed by pump (kW)
    """
    d = np.zeros(press.shape)
    d[press <= pmin] = d0
    d[press > pmin] = d0 * (pmax - press[press > pmin])/(pmax - pmin)
    power = press * d * hz / 1000 # kW
    return power # kW


def loads(vessel, elecfreezer, start, end, track_err=True, cond_in=None):
    """
    A function for calculating loads at each time step between 'start' and 'end.'
    :param vessel: a vessel object with vdr data (designed for the object saved as myriad.p)
    :param elecfreezer: a general csv object storying ac logger recordings for a freezer system (amps)
    :param start: earliest data to be considered from myriad
    :param end: latest data to be considered from myriad
    :param track_err: If True, print warnings
    :return load_res: dict of load categories and associated load arrays (kW)
    :return cond: array indicating timesteps between start and end when the main engine was running
    """
    # DC -------------------------------------
    dcload = 1 * KW_PER_HP
    # Deck -----------------------------------
    myr_dt = vessel.vdr['date_time']
    timecond = (end >= myr_dt) & (myr_dt >= start)
    if cond_in is None:
        cond_in = np.ones(timecond.shape, dtype=bool)
    cond = cond_in & timecond
    press = vessel.vdr['hyd_deck'][cond] * PASCAL_PER_PSI
    hz = vessel.vdr['tach'][cond] / 60
    d0 = getattr(vessel, 'deck_hydraulics_displacement_[cm3/rev]')
    d0 *= 1e-6 # m^3
    deckload = hyd_pow(press, hz, d0, pmin=250*PASCAL_PER_PSI)
    # Hydraulic Freezer ----------------------------------
    press = vessel.vdr['hyd_freezer'][cond] * PASCAL_PER_PSI
    d0 = getattr(vessel, 'freezer_hydraulics_displacement_[cm3/rev]')
    d0 *= 1e-6 # m^3
    hydfreezeload = hyd_pow(press, hz, d0)
    # Electric freezer -----------------------------------
    # note: this is a bit tricky.
    # First, the distribution of AC loads for Sea Miner is computed
    pf = 0.84
    volt = 211
    power = elecfreezer.total_amps * pf * np.sqrt(3) * volt/1000
    # The loads from both systems are indexed from least to greatest
    # using argsort
    power = power[power>0]
    hydpow = hydfreezeload[hydfreezeload > 0]
    hydsort = np.argsort(hydpow)
    esort = np.argsort(power)
    elecfreezeload = np.zeros(hydfreezeload.shape)
    # Finally, for each time step an electric freezer load is chosen
    # that has an equivalent position in the electric freezer load
    # cumulative distribution to the position of the hydrualic load
    # in the hydrualic load cumulative distribution.
    # This preserves the time correlation recorded for the hydraulic vessel
    # and the load distribution recorded for electric vessel
    # The difference in the mean freezer load for FV Sea Miner
    # using this method versus using the measurements directly is
    # 0.00017 kW
    powind = 0
    for ind in range(0, len(hydfreezeload)):
        if hydfreezeload[ind] != 0:
            fractional_ind = (hydsort[powind]+1)/len(hydsort)
            elecfreezeload[ind] = power[esort[int(round(fractional_ind*len(esort)-1))]]
            powind += 1
    # Propulsion -----------------------------------------
    sp = vessel.vdr['speed'][cond]
    tach = vessel.vdr['tach'][cond]
    tach_power_coeffs = vessel.tach_power_coeffs
    propload = power_of_tach(tach, tach_power_coeffs, track_err=track_err, speed=sp) * KW_PER_HP

    load_res = {'dc': dcload, 'deck': deckload, 'hfreeze': hydfreezeload, 'efreeze': elecfreezeload, 'prop': propload}
    return load_res, cond


def load_timeseries(load_res):
    """
    A function to return total load arrays for hydraulic and electric freezer, and non-refrigeration vessels
    :param load_res: dict of load categories and associated load arrays (kW)
    :return hydtotal: timeseries of total loads on a hydraulic freezer vessel
    :return electotal: timeseries of total loads on an electric freezer vessel
    """
    hydcats = ['dc', 'deck', 'hfreeze', 'prop']
    hydtotal = np.zeros(load_res['deck'].shape)
    for cat in hydcats:
        hydtotal += load_res[cat]

    ecats = ['dc', 'deck', 'efreeze', 'prop']
    electotal = np.zeros(load_res['deck'].shape)
    for cat in ecats:
        electotal += load_res[cat]

    icecats = ['dc', 'deck', 'prop']
    icetotal = np.zeros(load_res['deck'].shape)
    for cat in icecats:
        icetotal += load_res[cat]
    return hydtotal, electotal, icetotal


def night_day(vessel, cond):
    """
    A function to define night and day conditions and hours
    :param vessel: A vessel object with necessary vdr fields
    :param cond: a boolean array, True for times to be included in the analysis
    :return daycond: a boolean array True for times defined as 'day'
    :return nightcond: a boolean array True for times defines as 'night'
    :return day_hrs: total number of hours defined as 'day' in the set
    :return night_hrs: total number of hours defined as 'night' in the set
    """
    daycond = (vessel.vdr['tach'][cond] > 200)
    day_hrs = 2000
    nightcond = (vessel.vdr['tach'][cond] < 200)
    night_ratio = sum(nightcond) / sum(daycond)
    night_hrs = day_hrs * night_ratio
    return daycond, nightcond, day_hrs, night_hrs


def vessel_docked(vessel, cond, timestep = 1/30):
    """
    A function to define when a vessel was docked
    :param vessel: A vessel object with necessary vdr fields
    :param cond: a boolean array, True for times to be included in the analysis
    :param timestep: length of time step in recording (hrs)
    :return docked: a boolean array True for times when the vessel was docked
    """
    t = vessel.vdr['tach'][cond]
    docked = np.zeros(cond.shape, dtype=bool)
    counter = 0
    for ind in range(0, len(t)):
        if t[ind] < 200:
            counter += 1
        else:
            counter = 0
        if counter * timestep > 24:
            docked[max(0, ind-counter):ind] = True
    return docked

def fuel_of_load(load, bsfc_coeffs, hrs):
    """
    A function to calculate the fuel consumption of an engine given its bsfc coeffs,
    hrs of operation, and the applied load.
    :param load: Distribution of total load on the engine [hp--should match bsfc_coeffs units]
    :param bsfc_coeffs: coefficients of BSFC curve for the engine (gal/hr = polyval(load, bsfc_coeffs))
    :param hrs: total hours put on the engine
    :return fuel: total gallons consumed in the given time
    """
    counts, bins = np.histogram(load, bins=100)
    counts = counts / sum(counts)
    bincenters = (bins[:-1] + bins[1:])/2
    fuel = sum(polyval(bincenters, bsfc_coeffs) * counts) * hrs
    return fuel


def season_fuel(vessel, battery, pmain_in, paux_in, aux_cutoff, start, end, load_res, cond,
                track_err=True, print_results=False):
    """
    A function to calculate the fuel consumption of hybrid and basic
    propulsion systems with hybrid, electric, and ice-based cooling systems.
    :param load_res: dict of load categories and loads (kW)
    :param cond: boolean array True for datetimes in vessel to be used in the analysis
    :param vessel: a vessel object with vdr data (designed for the object saved as myriad.p)
    :param battery: a dict defining a battery system
    :param pmain_in: Main engine bsfc coeffs (reversed by fcn)
    :param paux_in: Auxiliary engine bsfc coeffs (reversed by fcn)
    :param aux_cutoff: Maximum power supplied by auxiliary
    :param start: Earliest datetime to be included from Myriad data.
    :param end: Latest datetime to be included from Myriad data.
    :param track_err: Print warnings if True
    :pram print_results: Print a summary table at end of function if True
    :return fuel: Dict of fuel consumed in each scenario (gal)
    :return auxhrs: Array of hours put on the auxiliary in each scenario
    :return mainhrs: Array of hours put on the main in each scenario
        note: order of returns is electric basic, electric hybrid,
            hydrualic basic, hydraulic hybrid
            
    Definition of scenarios:
    1. 47’ Troll vessel without refrigeration
        a. Base case—single engine vessel with main engine powering all propulsion, deck hydraulic and AC/DC electric 
            needs.
        b. Battery powered hybrid propulsion— Main engine powers a 5 KW DC generator to charge batteries in addition to 
            all propulsion, deck hydraulic and AC/DC electric needs when operating.  In Hybrid mode, 10 KW/hr Lithium 
            ion batteries power 15KW electric propulsion motor, 4kw hydraulic power pack, and provide AC/DC hotel needs.
            Batteries provide power for approximately X hours of hybrid operation then require XX hours of recharging.

    2. 47’ troll vessel with refrigeration
        a. Single engine with electric powered blast freeze system— Main engine powers a 15 KW DC generator to power 
            refrigeration system in addition to all propulsion, deck hydraulic and AC/DC electric needs.  
            Main engine operates 24 hrs/day.
        b. Single engine with hydraulic powered blast freeze system— Main engine powers a 6 cubic inch pressure 
            compensating hydraulic pump to power refrigeration system in addition to all propulsion, deck hydraulic and 
            AC/DC electric needs.  Main engine operates 24 hrs/day.
        c. Single engine with battery powered hybrid propulsion and electric powered blast freeze system— Main engine 
            powers a 15 KW DC generator to power refrigeration system in addition to all propulsion, deck hydraulic and 
            AC/DC electric needs. In Hybrid mode, 10 KW/hr Lithium ion batteries power 15KW electric propulsion motor, 
            4kw hydraulic power pack, and provide AC/DC hotel needs.  Batteries provide power for approximately X hours 
            of hybrid operation then require XX hours of recharging.

    3.  47’ troll vessel with refrigeration and 30 KW gen set
        a. 30 Kw gen set powers refrigeration system; main engine powers all propulsion, deck hydraulic and AC/DC 
            electric needs.  Gen set operates 24 hrs/day.  Main engine also operating to power all propulsion, deck 
            hydraulic and AC/DC electric needs when in fishing mode
        b. Gen set or main engine power hydraulic refrigeration system; main engine powers all propulsion, deck 
            hydraulic and AC/DC electric needs when in fishing mode, gen set provides hotel AC and refrigeration 
            hydraulic power at night.  Only 1 engine operating at a time.
        c.  Hybrid drive—Gen set powers 15 KW electric propulsion motor when in fishing mode in addition to 
            refrigeration, deck hydraulic, and AC/DC loads.  Main engine used for peak loads and transit to and from 
            fishing grounds.  15 KW propulsion motor acts as AC generator when main engine engaged providing electric 
            power for refrigeration.  This configuration only requires a single engine to operate at a time and provides 
            redundant refrigeration and propulsion power.
        d. Hybrid drive with batteries-- Gen set powers propulsion motor when in fishing mode in addition
            to electric refrigeration system, deck hydraulic, AC/DC loads, and a 16KW/hr of battery charging capacity. 
            Main engine  used for peak loads and transit to and from fishing grounds.  15 KW propulsion motor acts as AC 
            generator when main engine engaged providing electric power for refrigeration.  This configuration only  
            requires a single engine to operate at a time and provides redundant refrigeration and propulsion power.   
            Batteries provide power for approximately X hours of hybrid operation then require XX hours of recharging.
    """
    dc_gen_eff = 0.9
    timestep = 1/30 #hrs

    fuel, auxhrs, mainhrs = {}, {}, {}
    # Create load timeseries ------------------------------------
    hydtotal, electotal, icetotal = load_timeseries(load_res)

    # Put bsfc coeffs in the correct order ---------------------
    pmain = np.array(list(reversed(pmain_in)))
    paux = np.array(list(reversed(paux_in)))

    # Define night and day conditions ---------------------------
    daycond, nightcond, day_hrs, night_hrs = night_day(vessel, cond)
    # Define times when vessel was at the dock -------------------------------------------------------------------------
    docked = vessel_docked(vessel, cond)

    # Troll vessel w/o refrigeration, single engine --------------------------------------------------------------------
    label = '1a'
    fuel[label] = fuel_of_load(icetotal[daycond], pmain, day_hrs)
    mainhrs[label] = day_hrs
    auxhrs[label] = 0

    label = '1b'
    eng_load, batt_charge = battery_load_calculator(icetotal[daycond], battery['batt_cap'], docked=docked)
    maincond = eng_load > 0
    mainhrs[label] = sum(maincond) / len(maincond) * day_hrs
    fuel[label] = fuel_of_load(eng_load[maincond], pmain, mainhrs[label])
    auxhrs[label] = 0
    mainhrs[label] = sum(maincond) * timestep

    # Troll vessel with refrigeration, single engine -------------------------------------------------------------------
    label = '2a'
    total_load = electotal + 0 #load_res['efreeze'] * (1 - dc_gen_eff)
    print(max(total_load), np.average(total_load), min(total_load))
    fuel[label] = fuel_of_load(total_load, pmain, day_hrs+night_hrs)
    mainhrs[label] = day_hrs + night_hrs
    auxhrs[label] = 0

    label = '2b'
    fuel[label] = fuel_of_load(hydtotal, pmain, day_hrs + night_hrs)
    mainhrs[label] = day_hrs + night_hrs
    auxhrs[label] = 0

    label = '2c'
    eng_load, _ = battery_load_calculator(electotal, battery['batt_cap'], docked=docked)
    maincond = eng_load > 0
    print(max(eng_load[maincond]), np.average(eng_load[maincond]), min(eng_load[maincond]))
    mainhrs[label] = sum(maincond)/len(maincond) * (day_hrs+night_hrs)
    fuel[label] = fuel_of_load(eng_load[maincond], pmain, mainhrs[label])
    auxhrs[label] = 0

    # Troll vessel with refrigeration and auxiliary engine -------------------------------------------------------------
    label = '3a'
    auxfuel1 = fuel_of_load(load_res['efreeze'][daycond], paux, day_hrs)
    auxfuel2 = fuel_of_load(electotal[nightcond], paux, night_hrs)
    mainload = load_res['dc'] + load_res['deck'] + load_res['prop']
    mainfuel = fuel_of_load(mainload[daycond], pmain, day_hrs)
    fuel[label] = auxfuel1 + auxfuel2 + mainfuel
    mainhrs[label] = day_hrs
    auxhrs[label] = day_hrs + night_hrs

    label = '3b'
    auxfuel = fuel_of_load(hydtotal[nightcond], paux, night_hrs)
    mainfuel = fuel_of_load(hydtotal[daycond], pmain, day_hrs)
    fuel[label] = auxfuel + mainfuel
    mainhrs[label] = day_hrs
    auxhrs[label] = night_hrs

    label = '3c'
    auxcond = total_load < aux_cutoff
    maincond = total_load >= aux_cutoff
    auxhrs[label] = (night_hrs + day_hrs) * sum(auxcond)/len(auxcond)
    auxfuel = fuel_of_load(electotal[auxcond], paux, auxhrs[label])
    mainhrs[label] = night_hrs + day_hrs - auxhrs[label]
    mainload = electotal[maincond] + load_res['efreeze'][maincond] * (1 - dc_gen_eff)
    mainfuel = fuel_of_load(mainload, pmain, mainhrs[label])
    fuel[label] = mainfuel + auxfuel

    label = '3d'
    eng_load, _ = battery_load_calculator(electotal, battery['batt_cap'], docked=docked)
    auxcond = (eng_load <= aux_cutoff) & (eng_load > 0)
    maincond = (eng_load >= aux_cutoff) & (eng_load > 0)
    auxhrs[label] = sum(auxcond) / len(auxcond) * (day_hrs + night_hrs)
    mainhrs[label] = sum(maincond) / len(maincond) * (day_hrs + night_hrs)
    auxfuel = fuel_of_load(eng_load[auxcond], paux, auxhrs[label])
    mainfuel = fuel_of_load(eng_load[maincond], pmain, mainhrs[label])
    fuel[label] = mainfuel + auxfuel
    return fuel, auxhrs, mainhrs, batt_charge, eng_load


def battery_load_calculator(load, batt_cap, batt_cutoff=None, docked=None):
    """
    A function for calculating the load placed on an engine through time based on a vessel's total load time series and
    battery properties.
    :param load: a time series of total load (kW)
    :param batt_cap: battery capacity (kWh)
    :param batt_cutoff: max discharge rate of the batteries--calculated if left as None (kW)
    :param docked: boolean array, True if vessel is at the dock.
    :return eng_load: time series of the load placed on the engine
    """
    # Define battery characteristics --------------------------------
    # The max conditions enforces a max charging rate of 50 kW
    charge_time = max(1, batt_cap / 50)  # hr
    discharge_time = 1  # hr
    charge_rate = batt_cap / charge_time  # kW
    if batt_cutoff is None:
        batt_cutoff = batt_cap / discharge_time  # kW
    batt_eff = 0.9

    # Initialize variables ---------------------------------------------------------------------------------------------
    batt_charge = np.zeros(load.shape[0]+1)
    eng_load = np.zeros(load.shape)
    charging = False
    timestep = 1 / 30  # hrs
    for ind in range(0, len(load)):
        # Calculating battery state and fuel consumption ------------
        if docked is not None and docked[ind]:
            eng_load[ind] = 0
            if batt_charge[ind] < batt_cap:
                charging = True
                batt_charge[ind+1] = batt_charge[ind] + charge_rate * timestep
                if batt_charge[ind + 1] > batt_cap:
                    charging = False
        elif load[ind] < batt_cutoff:
            if batt_charge[ind] > 0 and not charging:
                batt_charge[ind+1] = (batt_charge[ind] - load[ind] * timestep)
                eng_load[ind] = 0
            else:
                charging = True
                batt_charge[ind+1] = batt_charge[ind] + charge_rate * timestep
                eng_load[ind] = charge_rate / batt_eff + load[ind]
                if batt_charge[ind+1] >= batt_cap:
                    charging = False
        else:
            if batt_charge[ind] < batt_cap:
                eng_load[ind] = load[ind] + charge_rate / batt_eff
                batt_charge[ind+1] = batt_charge[ind] + charge_rate * timestep
                if batt_charge[ind+1] > batt_cap:
                    charging = False
            else:
                eng_load[ind] = load[ind]
    return eng_load, batt_charge

def sf_battery(vessel, efreeze, pmain_in, batt_cap, start, end, load_res, cond,
               batt_cutoff=None, track_err=True, print_results=False):
    """
    A function for estimating the seasonal fuel consumption of a battery-hybrid system.
    :param vessel: a vessel object with vdr data (designed for the object saved as myriad.p)
    :param elecfreezer: a general csv object storying ac logger recordings for a freezer system (amps)
    :param load_res: dict of load categories and loads (kW)
    :param cond: boolean array True for datetimes in vessel to be used in the analysis
    :param pmain_in: Main engine bsfc coeffs (reversed by fcn)
    :param start: Earliest datetime to be included from Myriad data.
    :param end: Latest datetime to be included from Myriad data.
    :param track_err: Print warnings if True
    :param print_results: Print a summary table at end of function if True
    :param batt_cap: battery capacity (kWh)
    :param batt_cutoff: max discharge rate of the batteries--calculated if left as None (kW)
    :return fuel: Array of fuel consumed in each scenario (gal)
    :return auxhrs: Array of hours put on the auxiliary in each scenario
    :return mainhrs: Array of hours put on the main in each scenario
    """

    # Create load timeseries ------------------------------------
    hydtotal, electotal = load_timeseries(load_res)

    # Put bsfc coeffs in the correct order ---------------------
    pmain = np.array(list(reversed(pmain_in)))

    # Define night and day conditions ---------------------------
    daycond, nightcond, day_hrs, night_hrs = night_day(vessel, cond)

    # Define battery characteristics --------------------------------
    # The max conditions enforces a max charging rate of 50 kW
    charge_time = max(1, batt_cap / 50)  # hr
    discharge_time = 1  # hr
    charge_rate = batt_cap / charge_time  # kW
    if batt_cutoff is None:
        batt_cutoff = batt_cap / discharge_time  # kW
    batt_eff = 0.9
    mainhrs = 0
    batt_charge = batt_cap
    charging = False

    # Initiate variables -----------------------------------------
    tach_in = vessel.vdr['tach'][cond]
    fuel = []
    timestep = 1 / 30  # hrs
    dockcharge = 0
    run_hrs = sum(daycond) * timestep
    for ind in range(0, sum(cond)):
        # Calculating battery state and fuel consumption ------------
        if electotal[ind] < batt_cutoff:
            if batt_charge > 0 and charging == False:
                fuel.append(0)
                batt_charge -= electotal[ind] * timestep
            else:
                charging = True
                fuel.append(polyval(electotal[ind] + charge_rate / batt_eff, pmain) * timestep)
                batt_charge += charge_rate * timestep
                if batt_charge > batt_cap:
                    charging = False
        else:
            if batt_charge < batt_cap:
                fuel.append(polyval(electotal[ind] + charge_rate / batt_eff, pmain) * timestep)
                batt_charge += charge_rate * timestep
                if batt_charge > batt_cap:
                    charging = False
            else:
                fuel.append(polyval(electotal[ind], pmain) * timestep)
        if fuel[-1] > 0:
            mainhrs += timestep
        test_ind = max(ind - int(24 / timestep), 0)
        if sum(tach_in[test_ind:ind]) < 200:
            if batt_charge < batt_cap * 0.5:
                dockcharge += 1
            batt_charge = batt_cap
            charging = False
            # if ind % 20 == 0:
            #    print(batt_charge, charging, fuel[-1])
    # Multiplying by day_hrs / run_hrs gives the fuel consumed throughout a season
    # with the number of propulsion hours specified by day_hrs
    main_fuel = sum(polyval(electotal, pmain) * timestep) * day_hrs / run_hrs
    fuel = [sum(fuel) * day_hrs / run_hrs]
    batthrs = [night_hrs + day_hrs]
    mainhrs = [mainhrs * day_hrs / run_hrs]
    return fuel, batthrs, mainhrs

