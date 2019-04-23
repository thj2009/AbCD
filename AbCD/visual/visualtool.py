import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import FormatStrFormatter

def parity_plot(xdata, ydata, xlabel='', ylabel='',
                parity_line=True, rmse=True, scale='log10',
                expand=0.1, figsize=(2.5, 2.5), dpi=300, figname=None):
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
    ax.set_xlabel(xlabel, fontsize=8)
    ax.set_ylabel(ylabel, fontsize=8)
    ax.tick_params(axis='all', labelsize=7)
    ax.legend(frameon=False, fontsize=6)
    plt.tight_layout()
    plt.show()
    if figname is not None:
        fig.savefig(figname)
    return fig

def tpd_profile(Tem, Rate, xlabel='', ylabel='',
                figsize=(3.5, 2.5), dpi=300, figname=None,
                r_=None, fig=None, fmt='r-', scale=1):
    if fig is None:
        fig = plt.figure(figsize=figsize, dpi=dpi)
    ax = fig.gca()
    Rate = np.array(Rate) * scale
    ax.plot(Tem, Rate, fmt, linewidth=1)
    if r_ is None:
        r_ = [min(Tem), max(Tem)]
    ax.set_xlim(r_)
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.1e'))
    ax.set_xlabel(xlabel, fontsize=8)
    ax.set_ylabel(ylabel, fontsize=8)
    ax.tick_params(axis='all', labelsize=6)
    plt.tight_layout()
    plt.show()
    if figname is not None:
        fig.savefig(figname)
    return fig


