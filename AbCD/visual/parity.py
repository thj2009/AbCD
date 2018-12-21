import matplotlib.pyplot as plt
import numpy as np

def parity_plot(xdata, ydata, xlabel='', ylabel='',
                parity_line=True, rmse=True, scale='log10',
                expand=0.1, figsize=(3.25, 3.25), dpi=300, figname=''):
    '''
    Parity plot for model fitting
    '''
    if scale == 'log10':
        xdata = np.log10(np.array(xdata))
        ydata = np.log10(np.array(ydata))
    elif scale == 'log':
        xdata = np.log10(np.array(xdata))
        ydata = np.log10(np.array(ydata))
    elif scale == 'norm':
        xdata = np.array(xdata)
        ydata = np.array(ydata)

    min_ = min((min(xdata), min(ydata)))
    max_ = max((max(xdata), max(ydata)))
    d_ = max_ - min_
    r_ = [min_ - expand * d_, max_ + expand * d_]

    if rmse:
        metric = (xdata - ydata).std()
    fig = plt.figure(figsize=figsize, dpi=dpi)
    ax = fig.gca()
    if parity_line:
        ax.plot(r_, r_, 'k:', lw=1.5)
    ax.plot(xdata, ydata, 'ro', ms=2, label='RMSE=%.2e' %(metric))
    ax.set_xlim(r_)
    ax.set_ylim(r_)
    ax.set_xlabel(xlabel, fontsize=9)
    ax.set_ylabel(ylabel, fontsize=9)
    ax.tick_params(axis='all', labelsize=8)
    ax.legend(frameon=False, fontsize=6)
    plt.tight_layout()
    if not figname:
        fig.savefig(figname)
    return fig
