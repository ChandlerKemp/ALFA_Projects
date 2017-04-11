"""
This module contains functions intended to help with data processing.
"""
import datetime
import numpy as np
import matplotlib.pyplot as plt
import copy

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
    """
    Convert a string to a datetime Object
    """
    #--------Determine whether the string is a date, a time or both------

    if len(str_in) <= 10 and ':' not in str_in:
        #Assume that time is not included
        strtype = 'date'
        temp = ['']*3
    elif '/' not in str_in and '-' not in str_in:
        # assume that date is not included
        strtype = 'time'
        if '.' in str_in:
            temp = ['']*5
        else:
            temp = ['']*4
    elif '.' in str_in:
        strtype = 'datetime'
        temp = ['']*8
    else:
        strtype = 'datetime'
        temp = ['']*7

    #--------Define the position where the time part of the string begins
    if strtype == 'datetime':
        timestart = 3
    else:
        timestart = 0

    #--------Read the string into a list
    index = 0
    test = 0
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
    #----------Convert the list into a datetime object
    if strtype != 'date':
        if temp[-1].lower() == 'pm' and int(temp[timestart]) != 12:
            temp[timestart] = int(temp[timestart]) + 12
        elif temp[-1].lower() == 'am' and temp[timestart] == 12:
            temp[timestart] = 0
        elif temp[-1].lower() != 'am' and temp[-1].lower() != 'pm' and temp[-1].lower() != '':
            raise ValueError('AM/PM not read correctly from file')
        second = int(temp[timestart+2])
        if '.' in str_in:
            microsec = int(float(temp[timestart+3])*1e6)
        else:
            microsec = 0
    if strtype != 'time':
        temp[0] = month_to_int(temp[0])
    if strtype == 'date':
        return datetime.date(int(temp[2]), int(temp[0]), int(temp[1]))
    elif strtype == 'time':
        return datetime.time(int(temp[0]), int(temp[1]), second, microsec)
    else:
        return datetime.datetime(int(temp[2]), int(temp[0]), int(temp[1]), int(temp[3]), int(temp[4]), second, microsec)



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

def fit_plotter(x, y, deg=2, xlabel='x', ylabel='y', checkfit=True, fit_type = 'default'):
    """
    Calculates a polynomial fit of degree 'deg' of x to y
    :param x: list or array of the independent variable 
    :param y: list or array of the dedependent variable
    :param deg: degree of polymonial
    :param xlabel: x-axis label for plots
    :param ylabel: y-axis label for plots
    :param checkfit: prints a plot of the fit if True
    :param fit_type: defines the fitting method to use. Can be 'default' or 'constrained'
    :return: p: list of coefficients-- p
    """
    if fit_type.lower() == 'default':
        p = np.polyfit(x, y, deg)
    elif fit_type.lower() == 'constrained':
        '''The constrained fit is a basic optimization algorithm. The constraints are dy/dx >= 0 for x > 0, y(x=0)=0
            The constraint is enforced by b >= -sqrt(12 a c), where y=ax + bx**2 + cx**3 (c and a > 0)
            The method is two find a least squares fit of the form ax + bx**2 + cx**3. If the condition is not met by the
            fit, b is set equal to -sqrt(12 a c). The algorithm then walks along that curve by varying b and a until a 
            local minimum is found.
        '''
        if type(x) is list:
            x = np.array(x)
            y = np.array(y)
        W = np.zeros([deg,len(x)])
        for ind in range(0, deg):
            W[ind,:] = x**(ind+1) 
        # Solve for least squares solution
        p = np.linalg.solve(np.dot(W,np.transpose(W)),np.dot(W,y))
        nsteps = 20
        amax = 100
        amin_master = -10
        cmin_master = 0
        cmax = 10
        bmax_master = 50
        bmin = None
        a = np.linspace(amin_master,amax, nsteps)
        c = np.linspace(cmin_master,cmax, nsteps)
        asave = a[0]
        csave = c[0]
        bsave = 0
        def errfcn(a, b, c):
            y0 = np.average(y-(a*x + b*x**2+c*x**3))
            return sum((y-(y0 + a*x + b*x**2+c*x**3))**2)
        errmin = errfcn(a[0], 0, c[0])
        errmax = 2*errmin
        bmin_master = -1e6
        counter = 0
        delta_b = 10
        while (errmax - errmin)/errmin > 0.01:
            errmin = errmax # this is in case the the coarser iteration happened to find a better minimum
            errmax = 0
            for aval in a:
                for cval in c:
                    if bmin is None:
                        bmin = (-3*aval*cval - (min(x)*aval)**2)/(2*min(x))
                    bmin = max(bmin, bmin_master)
                    bmax = max(bmin, bmax_master)
                    for bval in np.linspace(bmin,bmax, nsteps):
                        err = errfcn(aval, bval, cval)
                        if err < errmin:
                            asave = aval
                            bsave = bval
                            csave = cval
                            delta_b = (bmax - bmin)/(nsteps-1)
                            errmin = err
                        errmax = max(err, errmax)
            counter += 1
            print("The summed-squared error after " + str(counter) + " iterations is" + '{:6.0f}'.format(errmin))
            delta_a = a[1] - a[0]
            delta_c = c[1] - c[0]
            a = np.linspace(max(amin_master, asave - delta_a), asave + delta_a, nsteps)
            c = np.linspace(max(cmin_master, csave - delta_c), csave + delta_c, nsteps)
            bmin_master = bsave - delta_b
            bmax_master = bsave + delta_b
        yint = np.average(y-(asave*x + bsave*x**2+csave*x**3))
        p = [csave, bsave, asave, yint]
        p_string = ["%.2f" % member for member in p]
        print("The final coefficiencts are: ", p_string)

    #     a = p[0]
    #     b = p[1]
    #     c = p[2]
    #     print(c)
    #     def errfcn(a, b, c):
    #         y0 = np.average(y-(a*x + b*x**2+c*x**3))
    #         return sum((y-(y0 + a*x + b*x**2+c*x**3))**2)
    #     def constrained_opt(a, b, c):
    #         deltac = c/1e3
    #         deltab = b/1e4
    #         deltaa = a/1e4
    #         #err = errfcn(a, b, c)
    #         #err1 = errfcn(a, b + deltab, c)
    #         err = errfcn(a, b, c)
    #         err1 = errfcn(a+deltaa, b, c)
    #         if err1 > err:
    #             #deltab = -deltab
    #             deltaa = -deltaa
    #         deltaerr = 1
    #         counter = 0
    #         xtest = np.linspace(0,15,100)
    #         while deltaerr > 0:
    #             #b += deltab
    #             #a = b**2/(3*c)
    #             a += deltaa
    #             #if b < -np.sqrt(3*a*c):
    #             if b < -(3*a*c + (min(x)*a)**2)/(2*min(x)):
    #                 b = (-3*a*c - (min(x)*a)**2)/(2*min(x))
    #             err1 = errfcn(a, b, c)
    #             deltaerr = err - err1
    #             counter += 1
    #             if counter > 5000:
    #                 print("Last two normalized errors: ", err/sum(y**2), err1/sum(y**2))
    #                 raise ValueError('Constrained fit did not converge after 5000 iterations.')
    #             err = copy.copy(err1)
    #             yint = np.average(y-(a*x + b*x**2+c*x**3))
    #             # plt.plot(x,y,'o')
    #             # plt.plot(xtest, yint + a*xtest + b*xtest**2 + c*xtest**3)
    #             # plt.show()
    #         err_0 = errfcn(a, b, c)
    #         # if errfcn(a, b, c) > errfcn(a,b,c+deltac):
    #         #     [yint, a, b, c] = constrained_opt(a, b, c+deltac)
    #         # elif errfcn(a, b, c) > errfcn(a, b, c-deltac):
    #         #     [yint, a, b, c] = constrained_opt(a, b, c-deltac)
    #         err_forward = errfcn(a,b,c+deltac)
    #         err_backward = errfcn(a,b,c-deltac)
    #         if not (err_0 < err_forward and err_0 < err_backward):
    #             newc = (err_0 - err_backward) * (err_forward - err_0) / (err_forward - 2 *err_0 + err_backward)
    #             newc *= deltac / 2 * (1/(err_0-err_backward) + 1 / (err_forward - err_backward))
    #             print(newc  )
    #             [yint, a, b, c] = constrained_opt(a, b, newc)
    #         yint = np.average(y-(a*x + b*x**2+c*x**3))
    #         return [yint, a, b, c]
    #     p = constrained_opt(a, b, c)
    #     p = list(reversed(p))
    #     #p.append(0)
    #     p = np.array(p)
    # #elif fit_type.lower() == 'spline':

    if checkfit:
        # Plot the fit to check results
        ssample = np.linspace(min(x)-2, max(x)+1, 100)
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
            print("Since you selected 'No', pairs of points will be averaged before producing the fit.")
            print("The number of data points recorded in each direction should be equal for this fit.")
            x2 = []
            y2 = []
            for i in range(0,int(len(x)/2)):
                x2.append((x[i]+x[i+int(len(x)/2)])/2)
                y2.append((y[i]+y[i+int(len(y)/2)])/2)
            p = fit_plotter(x2, y2, deg=deg, xlabel=xlabel, ylabel=ylabel, checkfit=checkfit, \
                fit_type = fit_type)
            return p
    else:
        return p