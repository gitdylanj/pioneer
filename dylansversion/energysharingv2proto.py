import uproot
import matplotlib.pyplot as plt
import hist as Hist
import numpy as np

def get_number_of_slits(spline_width, pitch, pixels_per_plane):

    '''
        Based on slit parameters, returns how many other slits an energy 
        deposit hit can share with on one side assuming the original hit was in
        center. 
    '''

    a = spline_width - ((pixels_per_plane/2) + pitch)
    b = pixels_per_plane + pitch

    return int(np.floor((a/b) + 1))

def adjacent_slits(pixel, n, pitch = 200, width = 100):
    '''
        Returns the pixel positions of adjacent slits out to n on either side of original
        slit. There is edge case control in this code that may need to be optimized. 
    '''
    slits = []
    for i in range(-n, n + 1):
        if abs(i) < n:
            slits.append(pixel + int(i * (pitch + width)))
        elif i < 0:
                slits.append(pixel + int((i + 1) * (pitch + width) - (pitch + width/2)))
        else:
            slits.append(pixel + int((i - 1) * (pitch + width) + (pitch + width/2)))  
                  
    return slits
    


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



def compute_energy_sharing(atar, spline, pixels_per_plane = 100, pitch = 200, splinewidth = 600):

    '''
    Returns a digitized form atar pixel data
    '''


    data = atar.arrays(['pixel_pdg', 'pixel_edep', 'pixel_time', 'pixel_hits'])
    pdg = data['pixel_pdg'][0]
    edep = data['pixel_edep'][0]
    time = data['pixel_time'][0]
    hits = data['pixel_hits'][0]

    for v in ['pixel_edep'].array():  
        
        for i, hit in enumerate(hits):
            #Loop over all hits and calculate energy sharing

            these_energies = []

            places = adjacent_slits(hit, n = 2)
                    
            for pos in places:
                energyi = spline(pos - hit) * edep[i]
                these_energies.append(energyi)
            
            #Normalize the energies to conserve energy
            these_energies = np.array(these_energies) * (edep[i]/np.sum(these_energies))

            # print("total Energy", np.sum(these_energies), edep[i])















    


