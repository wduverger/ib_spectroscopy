import subprocess
from functools import cache

import anasyspythontools as apt
import xarray as xr
from anasyspythontools.export import image_to_DataArray, spectrum_to_Dataset
from anasyspythontools.xarray_utils import create_projected_coords

from .map_transformations import *
from .plot_maps import *
from .plot_spectra import *

def read(path, cache=True):
    """
    Load an anasys .axz or .irb file, converting heightmaps and spectra 
    to xarray objects   
    """
    if cache:
        return read_cached(path)
    else:
        return read_nocache(path)

@cache
def read_cached(path):
    return read_nocache(path)

def read_nocache(path):
    """
    Load an anasys .axz or .irb file, converting heightmaps and spectra 
    to xarray objects   
    """

    doc = apt.read(path)

    if path.endswith('.irb'):
        doc_xr = xr.DataArray(
            doc.signal * float(doc.UnitScale) + float(doc.UnitOffset),
            dims='v',
            coords={
                'v': doc.wn,
                'AttenuatorPower': float(doc.AttenuatorPower[0]*100),
                'EnergyRange': float(doc.EnergyRange),
                'FileName': doc.FileName,
                'ID': doc.ID,
                'IRSweepResolution': float(doc.IRSweepResolution),
                'IRSweepSpeed': float(doc.IRSweepSpeed),
        })
        doc_xr.attrs = {'long_name':'Laser power', 'units':doc.Units}
        doc_xr.v.attrs = {'long_name':'Wavenumber', 'units':'cm⁻¹'}
        return doc_xr

    # Parse Heightmaps
    def parse_heightmap(h):
        h = force_2d_coords(create_projected_coords(image_to_DataArray(h)))
        h.attrs = {}
        return h

    doc._HeightMaps = doc.HeightMaps
    doc.HeightMaps = {
        k: parse_heightmap(v) 
        for k, v in doc._HeightMaps.items()
    }

    # Parse spectra
    def parse_spectrum(s):
        try: 
            s = spectrum_to_Dataset(s)
            s = s.rename_dims(wavenumbers='v')
            s = s.rename_vars(wavenumbers='v')
            s.attrs = dict(long_name='IR Amplitude', units='mV')
            s.v.attrs = dict(long_name='Wavenumber', units='cm⁻¹')
            return s
        except Exception as e:
            print(f'Spectrum "{s.Label.item()}" in document "{path}" failed to parse correctly')
            print(e)
            return None
            
    if 'RenderedSpectra' in dir(doc):
        doc._RenderedSpectra = doc.RenderedSpectra
        doc.RenderedSpectra = {
            k: parse_spectrum(v)
            for k, v in doc._RenderedSpectra.items()
        }
    
    return doc

def notify(title, text):
    """ Send a MacOS popup notification """
    
    CMD = '''
    on run argv
    beep 3
    display notification (item 2 of argv) with title (item 1 of argv) with sound
    end run
    '''
    subprocess.call(['osascript', '-e', CMD, title, text])