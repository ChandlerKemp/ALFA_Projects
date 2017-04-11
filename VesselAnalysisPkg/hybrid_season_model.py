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
    A function to return total load arrays for hydraulic and electric freezer vessels
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
    return hydtotal, electotal


def night_day(vessel, cond):
    """
    A function to define night and day conditions and hours
    :param vessel: A vessel object with necessary vdr fields
    :param cond: a boolean array True for times to be included in the analysis
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


def season_fuel(vessel, efreeze, pmain_in, paux_in, aux_cutoff, start, end, load_res, cond,
                track_err=True, print_results=False):
    """
    A function to calculate the fuel consumption of hybrid and basic
    propulsion systems with hybrid and electric freezers.
    :param load_res: dict of load categories and loads (kW)
    :param cond: boolean array True for datetimes in vessel to be used in the analysis
    :param vessel: a vessel object with vdr data (designed for the object saved as myriad.p)
    :param elecfreezer: a general csv object storying ac logger recordings for a freezer system (amps)
    :param pmain_in: Main engine bsfc coeffs (reversed by fcn)
    :param paux_in: Auxiliary engine bsfc coeffs (reversed by fcn)
    :param aux_cutoff: Maximum power supplied by auxiliary
    :param start: Earliest datetime to be included from Myriad data.
    :param end: Latest datetime to be included from Myriad data.
    :param track_err: Print warnings if True
    :pram print_results: Print a summary table at end of function if True
    :return fuel: Array of fuel consumed in each scenario (gal)
    :return auxhrs: Array of hours put on the auxiliary in each scenario
    :return mainhrs: Array of hours put on the main in each scenario
        note: order of returns is electric basic, electric hybrid,
            hydrualic basic, hydraulic hybrid
    """
    # Create load timeseries ------------------------------------
    hydtotal, electotal = load_timeseries(load_res)

    # Put bsfc coeffs in the correct order ---------------------
    pmain = np.array(list(reversed(pmain_in)))
    paux = np.array(list(reversed(paux_in)))

    # Define night and day conditions ---------------------------
    daycond, nightcond, day_hrs, night_hrs = night_day(vessel, cond)

    # Aux load with electric freezer during no-propulsion times -----------------------
    counts, bins = np.histogram(electotal[nightcond])
    counts = counts / sum(counts)
    bincenters = (bins[:-1] + bins[1:]) / 2
    elec_night = sum(polyval(bincenters, paux) * counts) * night_hrs

    # Aux load with hydraulic freezer during no-propulsion times -----------------------
    counts, bins = np.histogram(hydtotal[nightcond])
    counts = counts / sum(counts)
    bincenters = (bins[:-1] + bins[1:]) / 2
    hyd_night = sum(polyval(bincenters, paux) * counts) * night_hrs

    # Fuel consumed by the main engine in an electric freezer boat -------------------------
    # No hybrid propulsion, freezer always powered by auxiliary
    counts, bins = np.histogram(electotal[daycond] - load_res['efreeze'][daycond], bins=100)
    counts = counts / sum(counts)
    bincenters = (bins[:-1] + bins[1:]) / 2
    main = sum(polyval(bincenters, pmain) * counts) * day_hrs

    # Fuel consumed by the auxiliary to power an electric freezer while the main is running
    # No hybrid propulsion, freezer always powered by auxiliary
    counts, bins = np.histogram(load_res['efreeze'][daycond], bins=100)
    counts = counts / sum(counts)
    bincenters = (bins[:-1] + bins[1:]) / 2
    aux_default = sum(polyval(bincenters, paux) * counts) * day_hrs
    elecbasic = main + aux_default + elec_night

    # Fuel consumed by a hybrid propulsion system with an electric freezer
    # This is fuel consumed by the main:
    cond = (electotal > aux_cutoff) & daycond
    counts, bins = np.histogram(electotal[cond] - load_res['efreeze'][cond], bins=100)
    counts = counts / sum(counts)
    bincenters = (bins[:-1] + bins[1:]) / 2
    mainhrs_eh = day_hrs * sum(cond) / sum(daycond)
    main = sum(polyval(bincenters, pmain) * counts) * mainhrs_eh
    # This is fuel consumed by the auxiliary while the main is running:
    counts, bins = np.histogram(load_res['efreeze'][cond], bins=100)
    counts = counts / sum(counts)
    bincenters = (bins[:-1] + bins[1:]) / 2
    aux_hrs1 = day_hrs * sum(cond) / sum(daycond)
    aux1 = sum(polyval(bincenters, paux) * counts) * aux_hrs1

    # This is fuel consumed by the auxiliary while the main is off:
    cond = (electotal <= aux_cutoff) & daycond
    counts, bins = np.histogram(electotal[cond], bins=100)
    counts = counts / sum(counts)
    bincenters = (bins[:-1] + bins[1:]) / 2
    aux_hrs2 = day_hrs * sum(cond) / sum(daycond)
    aux2 = sum(polyval(bincenters, paux) * counts) * aux_hrs2

    elechybrid = main + aux1 + aux2 + elec_night

    # Hydraulic freezer load
    cond = (hydtotal < aux_cutoff) & daycond
    counts, bins = np.histogram(hydtotal[daycond], bins=100)
    counts = counts / sum(counts)
    bincenters = (bins[:-1] + bins[1:]) / 2
    main2, hybrid = hybrid_fuel_savings(counts, bincenters, pmain, paux, aux_cutoff, day_hrs)
    hydbasic = main2 + hyd_night
    hydhybrid = hybrid + hyd_night

    auxhrs_eb = night_hrs + day_hrs
    auxhrs_eh = night_hrs + aux_hrs1 + aux_hrs2
    auxhrs_hb = night_hrs
    auxhrs_hh = night_hrs + day_hrs * sum(cond) / sum(daycond)
    mainhrs_eb = day_hrs
    # mainhrs_eh is already defined
    mainhrs_hb = day_hrs
    mainhrs_hh = day_hrs * (1 - sum(cond) / sum(daycond))

    fuel = np.array([elecbasic, elechybrid, hydbasic, hydhybrid])
    auxhrs = np.array([auxhrs_eb, auxhrs_eh, auxhrs_hb, auxhrs_hh])
    mainhrs = np.array([mainhrs_eb, mainhrs_eh, mainhrs_hb, mainhrs_hh])
    if print_results:
        print("Myriad main bsfc, 20 kW auxiliary: ")
        print('Electric freezer basic & ' + '{:0.0f}'.format(elecbasic) + ' \\\\')
        print('Electric freezer hybrid & ' + '{:0.0f}'.format(elechybrid) + ' \\\\')
        print('Hydraulic freezer basic & ' + '{:0.0f}'.format(hydbasic) + ' \\\\')
        print('Hydraulic freezer hybrid & ' + '{:0.0f}'.format(hydhybrid) + ' \\\\')
    return fuel, auxhrs, mainhrs


def sf_main_ref(vessel, efreeze, pmain_in, paux_in, aux_cutoff, start, end, load_res, cond,
                track_err=True, print_results=False):
    """
    A function to calculate the fuel consumption of a system in which the main
    can provide eletricity.
    :param vessel: a vessel object with vdr data (designed for the object saved as myriad.p)
    :param elecfreezer: a general csv object storying ac logger recordings for a freezer system (amps)
    :param load_res: dict of load categories and loads (kW)
    :param cond: boolean array True for datetimes in vessel to be used in the analysis
    :param pmain_in: Main engine bsfc coeffs (reversed by fcn)
    :param paux_in: Auxiliary engine bsfc coeffs (reversed by fcn)
    :param aux_cutoff: Maximum power supplied by auxiliary
    :param start: Earliest datetime to be included from Myriad data.
    :param end: Latest datetime to be included from Myriad data.
    :param track_err: Print warnings if True
    :param print_results: Print a summary table at end of function if True
    :return fuel: Array of fuel consumed in each scenario (gal)
    :return auxhrs: Array of hours put on the auxiliary in each scenario
    :return mainhrs: Array of hours put on the main in each scenario
        Note: Order of returns is hybrid, no auxiliary
    """

    # Create load timeseries ------------------------------------
    hydtotal, electotal = load_timeseries(load_res)

    # Put bsfc coeffs in the correct order ---------------------
    pmain = np.array(list(reversed(pmain_in)))
    paux = np.array(list(reversed(paux_in)))

    # Define night and day conditions ---------------------------
    daycond, nightcond, day_hrs, night_hrs = night_day(vessel, cond)

    # Fuel with electric freezer during no-propulsion times
    counts, bins = np.histogram(electotal[nightcond])
    counts = counts / sum(counts)
    bincenters = (bins[:-1] + bins[1:]) / 2
    main_night = sum(polyval(bincenters, pmain) * counts) * night_hrs
    elec_night = sum(polyval(bincenters, paux) * counts) * night_hrs

    # Fuel consumed by the main engine in an electric freezer boat
    # with freezer run on the main
    cond_temp = (electotal > aux_cutoff) & daycond
    main_hrs = sum(cond_temp) / sum(daycond) * day_hrs
    counts, bins = np.histogram(electotal[cond_temp], bins=100)
    counts = counts / sum(counts)
    bincenters = (bins[:-1] + bins[1:]) / 2
    main_day1 = sum(polyval(bincenters, pmain) * counts) * main_hrs
    main = sum(polyval(bincenters, pmain) * counts) * main_hrs

    # Fuel consumed by the auxiliary during propulsion times
    cond_temp = (electotal <= aux_cutoff) & daycond
    aux_hrs = sum(cond_temp) / sum(daycond) * day_hrs
    counts, bins = np.histogram(electotal[cond_temp], bins=100)
    counts = counts / sum(counts)
    bincenters = (bins[:-1] + bins[1:]) / 2
    main_day2 = sum(polyval(bincenters, pmain) * counts) * aux_hrs
    aux_day = sum(polyval(bincenters, paux) * counts) * aux_hrs

    total_fuel = elec_night + main + aux_day
    main_fuel = main_night + main_day1 + main_day2
    total_hrs = day_hrs + night_hrs
    aux_hy_hrs = aux_hrs + night_hrs
    main_hy_hrs = total_hrs - aux_hy_hrs
    fuel = np.array([total_fuel, main_fuel])
    auxhrs = np.array([aux_hy_hrs, 0])
    mainhrs = np.array([main_hy_hrs, total_hrs])
    return fuel, auxhrs, mainhrs



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

