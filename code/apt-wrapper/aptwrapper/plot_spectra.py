from imp import lock_held
from typing import Tuple, Union
import xarray as xr
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter
import numpy as np


def try_get_iramp(spectrum: Union[xr.DataArray, xr.Dataset]):
    """ 
    If given a Dataset, extract the IR Amp channel
    If given a DataArray, return it, assuming it is the IR channel
    """

    if isinstance(spectrum, xr.DataArray): return spectrum
    else: return spectrum['IR Amplitude (mV)']

def plot_spectrum(spectrum, ax=None, loc_ax=None, xincrease=False, title='', loc_kwargs=None, label=None, **kwargs):
    if ax is None: ax=plt.gca()
    if label is None and 'Label' in spectrum.coords:
        label = spectrum.Label.item()

    try_get_iramp(spectrum).plot(
        ax=ax,
        xincrease=xincrease,
        label=label,
        **kwargs
    )

    ax.set_title(title)

    # Plot location
    if loc_ax is not None:
        loc_kwargs = loc_kwargs or {}
        if 'c' in kwargs and 'c' not in loc_kwargs:
            loc_kwargs['c'] = kwargs['c']
        if 'color' in kwargs and 'color' not in loc_kwargs:
            loc_kwargs['color'] = kwargs['color']
        mark_location(spectrum, ax=loc_ax, **loc_kwargs)


def savgol(spectrum, window=11, degree=3, deriv=0, *args, **kwargs):
    return xr.apply_ufunc(savgol_filter, spectrum, window, degree, deriv, *args, **kwargs, keep_attrs=True)


def mark_location(
    spectrum: Union[xr.DataArray, xr.Dataset],
    ref_map=None,
    drift_speed=None,
    ax=None,
    jitter=None,
    marker="X",
    ec="w",
    lw=1,
    s=81,
    **kwargs
):
    """
    Indicate the location where a spectrum was acquired on a heightmap plotted in
    `ax`. `jitter` is the std of a normally distributed random displacement in um.
    This can be used if multiple spectra were acquired on top of each other.
    """

    if ax is None:
        ax = plt.gca()

    x, y = calculate_xy(spectrum, ref_map, drift_speed)

    if jitter is not None:
        x += np.random.randn() * jitter * 1e-3
        y += np.random.randn() * jitter * 1e-3

    return ax.scatter(x, y, marker=marker, ec=ec, lw=lw, s=s, **kwargs)

def calculate_xy(spectrum, ref_map=None, drift_speed=None) -> Tuple[float, float]:
    """
    Gets the xy location (in um) of a spectrum. Can apply drift correction.
    """

    x = spectrum['Location.X'].item()
    y = spectrum['Location.Y'].item()

    if drift_speed is None:
        return x, y

    else:
        assert ref_map is not None

        t1 = ref_map.TimeStamp.item()
        t2 = spectrum.TimeStamp.item()
        dt = (t2 - t1) * 1e-9 # convert ms into s
        dx = np.array(drift_speed) * dt

        return x - dx[0], y - dx[1]
