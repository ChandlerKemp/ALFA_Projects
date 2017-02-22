import matplotlib.pyplot as plt
import matplotlib


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
