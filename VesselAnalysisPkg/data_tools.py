"""
This module contains functions intended to help with data processing.
"""
import datetime
import numpy as np
import matplotlib.pyplot as plt

def month_to_int(str_in):
    """
    Converts a string month name to a number
    param: str_in   the name of a month
    """
    months = ['january', 'february', 'march','april','may','june','july','august','september','october','november','decemeber']
    str_in = str_in.lower()
    nchar = len(str_in)
    try:
        return float(str_in)
    except ValueError:
        for i in range(0,12):
            if months[i][0:nchar] == str_in:
                return i+1
    raise ValueError("Invalid month name")


def str_to_datetime(str_in):
    index = 0
    test = 0
    if '.' in str_in:
        temp = ['']*8
    else:
        temp = ['']*7
    for char in str_in:
        if char != ' ' and char != '/' and char != ':' and char != '-':
            if test == 1:
                index += 1
                test = 0
            try:
                float(temp[index])
                try:
                    float(char)
                except ValueError:
                    if temp[index] != '':
                        index += 1
            except ValueError:
                pass
            temp[index] += char
        else:
            test = 1
    if temp[-1].lower() == 'pm' and int(temp[3]) != 12:
        temp[3] = int(temp[3]) + 12
    elif temp[-1].lower() == 'am' and temp[3] == 12:
        temp[3] = 0
    elif temp[-1].lower() != 'am' and temp[-1].lower() != 'pm' and temp[-1].lower() != '':
        raise ValueError('AM/PM not read correctly from file')
    temp[0] = month_to_int(temp[0])
    second = int(temp[5])
    if '.' in str_in:
        microsec = int(float(temp[6])*1e6)
    else:
        microsec = 0
    return datetime.datetime(int(temp[2]), int(temp[0]), int(temp[1]), int(temp[3]),
                                       int(temp[4]),second,microsec)


def time_match(times, target, delta_t):
    """
    Find the time in times that is closest to the target
    :param times: A list of times
    :param target: the target time
    :param delta_t: the normal separation between times
    :return: index of the correct time
    """
    val = []
    for time in times:
        val.append(abs(time-target))
    index = val.index(min(val))
    if val[index] > delta_t/2:
        raise ValueError("Target time not found")
    return index

def set_kwargs_attrs(obj_in, kwargs, only_existing=True):
    """
    Function for overwriting parameters with key word arguments
    Inputs:
        obj_in          Object with parameters to be updated
        kwargs          Dict containing new parameter values
        only_existing   Determines with new parameters can be created with kwargs
    """
    for key in kwargs.keys():
        # If only_existing is true, only set attributes that already exist
        if only_existing:
            if not hasattr(obj_in, key):
                raise ValueError("Tried to set invalid attribute. Class: ", type(obj_in), 'attempted attribute:', key)
        setattr(obj_in, key, kwargs[key])

def fit_plotter(x, y, deg=2, xlabel='x', ylabel='y', checkfit=True):
    """
    Calculates a polynomial fit of degree 'deg' of x to y
    :param x: list or array of the independent variable 
    :param y: list or array of the dedependent variable
    :param deg: degree of polymonial
    :param xlabel: x-axis label for plots
    :param ylabel: y-axis label for plots
    :return: p: list of coefficients-- p
    """
    p = np.polyfit(x, y, deg)
    if checkfit:
        # Plot the fit to check results
        ssample = np.linspace(min(x)-1, max(x)+1, 100)
        psample = np.polynomial.polynomial.polyval(ssample, list(reversed(p)))
        plt.plot(x, y, 'o')
        plt.plot(ssample, psample)
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.show()
        cont = input("Is the fit acceptable? (yes or no) ")
        if cont.lower() == 'yes':
            return p
        else:
            # Alternative fit with x and y values paired and averaged before fitting
            # Useful for sea trial fits.
            x2 = []
            y2 = []
            for i in range(0,int(len(x)/2)):
                x2.append((x[i]+x[i+int(len(x)/2)])/2)
                y2.append((y[i]+y[i+int(len(y)/2)])/2)
            p = np.polyfit(x2, y2, deg)
            # Plot the fit to check results
            ssample = np.linspace(min(x2)-1, max(x2)+1, 100)
            psample = np.polynomial.polynomial.polyval(ssample, list(reversed(p)))
            plt.plot(x2, y2, 'o')
            plt.plot(ssample, psample)
            plt.xlabel(xlabel)
            plt.ylabel(ylabel)
            plt.show()
            cont = input("Is the fit acceptable? (yes or no) ")
            if cont.lower() == 'yes':
                return p
            else:
                raise ValueError("You must define the fit manually.")
    else:
        return p