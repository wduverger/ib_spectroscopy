from ast import Tuple
from binascii import a2b_hex
from email.policy import default
from typing import Callable, Optional
import xarray as xr
import matplotlib.pyplot as plt
from matplotlib_scalebar.scalebar import ScaleBar 
from matplotlib.offsetbox import AnchoredText
import numpy as np
from anasyspythontools.image import Image

## Basic plotting functions

# TODO: support Video images and calculated maps
def plot_map(
    hmap: xr.DataArray,
    ax=None,
    ticks=False,
    labels=False,
    robust=False,
    colorbar=True,
    cax=None,
    cbar_opts=None,
    cb_label=None,
    cmap=None,
    annotation=None,
    show_sb_val=True,
    sb_color=None,
    sb_loc='lower left',
    flatten=None,
    contrast=(0,1),  # for images: tuple of offset and slope
    force_pcolormesh=False,
    vmin=None,
    vmax=None,
    **kwargs):
    
    # Set default values to arguments that are None
    if not ticks: ticks = []
    elif not isinstance(ticks, list): ticks = None

    if ax is None: ax = plt.gca()
    if flatten is None: flatten = default_flatten
    if flatten is False: flatten = (lambda x: x)
    if cmap is None: cmap = default_cmap(hmap)
    if annotation is None: annotation = default_annotation(hmap)
    if sb_color is None: sb_color = default_sbcolor(hmap)

    if colorbar and (cbar_opts is None): 
        cbar_opts = dict(
            location='left', 
            label=cb_label if cb_label is not None else default_cblabel(hmap)
        )
    
    # Plotting of AFM images
    if isinstance(hmap, Image):

        # Apply flattening and contrast
        z = np.array(hmap.SampleBase64[..., 0:3]/255, dtype=float)  # Delete the alpha channel
        z = flatten(z)
        offset, slope = contrast
        z = np.clip((z-offset)*slope, 0, 1)

        # Plot the image
        mappable = ax.imshow(
            z, 
            extent=get_im_extent(hmap), 
            **kwargs
        )
    
    # Plotting of AFM maps
    else:

        # Set default values for colorbar
        
        # Flatten the data
        hmap = flatten(hmap)

        # Plot using pcolormesh (automatically adds a colorbar)
        if force_pcolormesh:
            mappable = hmap.plot(
                x='X', y='Y', 
                ax=ax,
                cmap=cmap, 
                robust=robust,
                add_colorbar=colorbar,
                cbar_ax=cax,
                cbar_kwargs=cbar_opts,
                vmin=vmin,
                vmax=vmax,
                **kwargs
            )

        # Plot using imshow
        else:

            # Set color boundaries
            vmin=vmin or (np.percentile(hmap.data, 2) if robust else hmap.min())
            vmax=vmax or (np.percentile(hmap.data, 98) if robust else hmap.max())

            # Determine whether to extend the colorbar
            if vmin > hmap.min() and vmax < hmap.max(): extend = 'both'
            elif vmax < hmap.max(): extend = 'max'
            elif vmin > hmap.min(): extend = 'min'
            else: extend = 'neither'

            # Plot the image
            mappable = ax.imshow(
                hmap,
                cmap=cmap or default_cmap(hmap),
                vmin=vmin,
                vmax=vmax,
                extent=get_im_extent(hmap),
                **kwargs
            )

            # Add colorbar
            if colorbar:
                plt.gcf().colorbar(mappable, ax=ax, extend=extend, **cbar_opts)

    # Plot styling
    annotate(ax, annotation, c=sb_color)
    add_scalebar(ax, color=sb_color, show_value=show_sb_val, loc=sb_loc)

    ax.set(aspect=1)
    ax.autoscale()  # automatically set axis x and y limits based on this image and all previous ones
    if not labels:
        ax.set(xlabel='', ylabel='', title='')
    if not ticks: 
        ax.set(xticks=[], yticks=[])

    return mappable

def get_im_extent(hmap):
    try: 
        width = float(hmap.Size.X)
        height = float(hmap.Size.Y)
        X0 = float(hmap.Position.X)
        Y0 = float(hmap.Position.Y)
        return (X0 - width / 2, X0 + width / 2, Y0 - height / 2, Y0 + height / 2)
    except AttributeError:
        return (hmap.X.min(), hmap.X.max(), hmap.Y.min(), hmap.Y.max())

def annotate(ax, text: Optional[str]=None, weight='heavy', c='white', loc='upper left', fontsize=None, ha='left', **kwargs):
    if text is None:
        print('Warning: no text supplied to annotate()')

    at = AnchoredText(
        text, 
        prop=dict(weight=weight, c=c, fontsize=fontsize, ha=ha), 
        loc=loc,
        **kwargs
    )
    at.patch.set(alpha=0)
    ax.add_artist(at)

def add_scalebar(
        ax=None,
        dx=1e-6,
        color='w',
        pad=.7,
        scale_loc='top',
        loc='lower left',
        box_alpha=0,
        font_weight="heavy",
        font_props=None,
        scale_formatter: Optional[Callable[[Tuple], str]] = None,
        show_value=True,
        length_fraction=.26,
        **kwargs
    ):
        if ax is None: ax = plt.gca()
        fontprops = font_props or dict(weight=font_weight)

        # Default scale formatter function. Replace latex mu with unicode version.
        def default_formatter(x, y):
            if y==r'$\mathrm{\mu}$m': y = 'µm'
            return f'{x} {y}'

        if scale_formatter is None: 
            if not show_value:
                scale_formatter = lambda x, y: ''
            else:
                scale_formatter = default_formatter

        sb_artist = ScaleBar(
            dx=dx,
            pad=pad,
            box_alpha=box_alpha,
            color=color,
            font_properties=fontprops,
            location=loc,
            scale_loc=scale_loc,
            scale_formatter=scale_formatter,
            length_fraction=length_fraction,
            **kwargs
        )

        ax.add_artist(sb_artist)

## Extended plotting functions for multiple maps

def plot_maps(fig, ax, doc, mapname, show_sb_val, height_frac=.5):
    mapname_retrace = mapname + ' (1)'
    hmap_trace = doc.HeightMaps[mapname]

    mappable = plot_map(hmap_trace, ax[0, 0], colorbar=False, annotation='', show_sb_val=show_sb_val)
    fig.colorbar(mappable, cax=ax[0, 1], label=default_cblabel(hmap_trace))

    line_i = int(np.floor(height_frac * hmap_trace.shape[0]))

    ax[1, 0].plot(default_flatten(hmap_trace)[line_i], lw=.5)
    ax[1, 0].set_xlim(0, hmap_trace.shape[1])
    if mapname_retrace in doc.HeightMaps:
        ax[1, 0].plot(default_flatten(doc.HeightMaps[mapname_retrace])[line_i], lw=.5)

    ax[1,1].axis('off')
    

def plot_iterations(doc, map_i, height_frac=.5):
    channels = [k for k in doc.HeightMaps if k.endswith(' ' +str( map_i))]
    n = len(channels)

    fig, ax = plt.subplots(2,2*n,figsize=(3.2*n,3.3), constrained_layout=True, gridspec_kw=dict(
        width_ratios=[1,.05]*n, height_ratios=[3,1]
    ))

    for i, c in enumerate(channels):
        j=i*2
        plot_maps(fig, ax[:, j:j+2], doc, c, show_sb_val=(i==0), height_frac=height_frac)

    hmap = doc.HeightMaps[channels[0]]
    ax[0, 0].set_title(f"({map_i}) {hmap['Tags.IRWavenumber'].item():.0f}, {hmap['Tags.IRPolarization'].item():.0f}º", font_properties={'weight': 'bold'})

def plot_all_maps(doc):
    heightmaps = [k.split(' ')[1] for k in doc.HeightMaps 
                  if 'Height' in k and not k.endswith(')') and 'Copy' not in k]
    for i in heightmaps:
        plot_iterations(doc, i)


## Defaults

def subtract_lines(hmap, mask=True, degree=1):
    fit_coeff = hmap.where(mask).polyfit('x', deg=degree).polyfit_coefficients
    fit = xr.polyval(hmap.x, fit_coeff)
    
    return hmap.copy(data=hmap-fit)

def default_flatten(hmap):
    if isinstance(hmap, Image) or isinstance(hmap, np.ndarray): return hmap
    if 'Label' not in hmap.coords: return hmap

    label = hmap.Label.item()
    if   'Height'     in label: return subtract_lines(hmap, degree=1) 
    elif 'Deflection' in label: return subtract_lines(hmap, degree=0) 
    else: return hmap

def default_cmap(hmap):
    if not isinstance(hmap, xr.DataArray) and not isinstance(hmap, xr.Dataset): return None
    if 'Label' not in hmap.coords: return None

    label = hmap.Label.item()
    if   'Height'        in label: return 'afmhot'
    elif 'Deflection'    in label: return 'bwr'
    elif 'IR Amplitude'  in label: return 'viridis'
    elif 'IR Phase'      in label: return 'bwr'
    elif 'PLL Frequency' in label: return 'bwr'
    else: return None

def default_cblabel(hmap):
    if not isinstance(hmap, xr.DataArray) and not isinstance(hmap, xr.Dataset): return ''
    if 'Label' not in hmap.coords: return None

    label = hmap.Label.item()
    if   'Height'        in label: return 'Height (nm)'
    elif 'Deflection'    in label: return 'Deflection (mV)'
    elif 'IR Amplitude'  in label: return 'IR Amplitude (mV)'
    elif 'IR Phase'      in label: return 'IR Phase (deg)'
    elif 'PLL Frequency' in label: return 'PLL Frequency (kHz)'

    return label

def default_annotation(hmap):
    if not isinstance(hmap, xr.DataArray) and not isinstance(hmap, xr.Dataset): return ''
    if 'Label' not in hmap.coords: return ''
    
    if 'IR Amplitude' in hmap.Label.item():
        return str(int(hmap['Tags.IRWavenumber'].item())) + ' cm⁻¹'
    return ''

def default_sbcolor(hmap):
    return 'k' if default_cmap(hmap) == 'bwr' else 'w'

def pixel_size_xy(hmap):
    x = (hmap.X.sel(x=1,y=0) - hmap.X.sel(x=0,y=0)).item()
    y = (hmap.Y.sel(x=0,y=0) - hmap.Y.sel(x=0,y=1)).item()
    return np.array([x,y])

## Utilities

def force_2d_coords(hmap):
    if len(hmap.X.shape)==2:
        return hmap.copy()

    nx, ny = len(hmap.X), len(hmap.Y)

    return hmap.assign_coords(
        X=(('y', 'x'), hmap.X.to_numpy().reshape(1,-1).repeat(ny,axis=0)),
        Y=(('y', 'x'), hmap.Y.to_numpy().reshape(-1,1).repeat(nx,axis=1)),
    )