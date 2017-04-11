import matplotlib.pyplot as plt
import matplotlib
import numpy as np


def display_settings(interactive=True, inline=True):
    """
    If running in IPYTHON, display_settings uses magic commands to plot interactively
    interactive         determines whether the plot is interactive
    """
    if interactive:
        # display settings depend on the python environment being used
        try:
            __IPYTHON__
            ipy = get_ipython()
            if inline:
                try:
                    ipy.magic("matplotlib inline")
                except:
                    ipy.magic("matplotlib")
            else:
                ipy.magic("matplotlib")
        except:
            plt.ion()


# display_save is called by other plotting functions
def display_save(display, save, file_out, fig=None):
    if fig is None:
        fig = plt.gcf()
    if save:
        if file_out is None:
            file_out = 'Figure.pdf'
        fig.savefig(file_out)

    if display:
        fig.show()
    else:
        fig.close()
    plt.ioff()


def ac_data_plot(data_dict,keys_to_plot,display=True,save=False, file_out=None):
    plt.figure()
    time_list = []
    t0 = data_dict['start_time'][0]
    for t in data_dict['start_time']:
        delta_t = t - t0
        time_list.append(delta_t.days + delta_t.seconds/86400)
    for key in keys_to_plot:
        plt.plot(time_list,data_dict[key],label=key)
    plt.xlabel("Time (days since start of trial)")
    plt.ylabel("Amps")
    plt.legend()
    display_save(display,save,file_out)

def plot_prep():
    font = {'weight' : 'bold',
        'size'   : 14}
    matplotlib.rc('font', **font)

def plot_fixer(ax=None):
    if ax is None:
        ax = plt.gca()
    ax.tick_params(axis='both', which='major', labelsize=14)
    lw = 4
    for ln in ax.lines:
        ln.set_linewidth(lw)
    plt.tight_layout()

def hist_plot(data, bins=None, xlabel=None, ylabel=None):
    """
        A function for producing histograms of data
        :param data: a list or an array of numerical data
        :param bins: an integer number of bins or a list of bin boundaries
        :param xlabel: x-axis label
        :param ylabel: y-axis label
        :return counts: array of counts for each bin
        :return bins: array of bin boundaries
        :return hndl: handle for bar plot
    """
    # define an array of bin boundaries
    if bins is None:
        bins = 30
    try:
        nbins = len(bins)-1
    except TypeError:
        nbins = bins
        bins = np.linspace(min(data),max(data), nbins)

    [counts, bins] = np.histogram(data, bins=bins)
    counts = counts/sum(counts)
    plot_prep()
    width = (bins[-1]-bins[0])/nbins
    hndl = plt.bar(bins[1:]-width/2,counts,width=width)
    if xlabel is None:
        xlabel = 'Bin center'
    if ylabel is None:
        ylabel = 'Fraction of time in bin'
    plt.ylabel(ylabel, fontsize=14, fontweight='bold')
    plt.xlabel(xlabel,fontsize=14, fontweight='bold')
    plot_fixer()
    return counts, bins, hndl