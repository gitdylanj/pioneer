import uproot
import matplotlib.pyplot as plt
import hist as Hist
import numpy as np


def adjacent_pixelsv2(pixel, n = 2):
    '''
        Returns this pixel and the n pixels on either side in the same plane. Why can't you write it like this instead of the one above?    
    '''
    
    return [x for x in range(pixel - n, pixel + n + 1)]



def build_spline(hist):
    '''
        Builds a cubic spline to interpolate the energy sharing
        distribution histogram
    '''
    from scipy.interpolate import UnivariateSpline, CubicSpline
    
    centers = hist.axes[0].centers

    amps = hist.values()/np.amax(hist.values())

    spline = CubicSpline(centers, amps, extrapolate = False)

    return spline



def compute_energy_sharing(tree, entry, spline, pixels_per_plane = 100, pitch = 200):

    '''
    Returns a dictionary of digitized atar data
    '''

output = {
    'pixel_pdg':[],
    'pixel_edep':[],
    'pixel_time':[],
    'pixel_hits':[],
}


