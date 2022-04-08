import uproot
import matplotlib.pyplot as plt
import hist as Hist
import numpy as np

class digitizeatar:

    '''
        This class computes energy sharing in the Atar and returns a modified root file.
    '''

    #The following are class variables based on Atar simulation parameters. 

    spline_width = 600
    pixels_per_plane = 200
    pixel_pitch = 200



    #Here is a given path to a spline './BNL_Signal_Response.root'

    def __init__(self, path_to_spline):
        
        #Energysharing spline path must be provided by the user. There is a default path below. 

        self.path_to_spline =  path_to_spline


    #For the paramaters in the atar simulation, the number of slits accessible on either side is two.

    def get_adjacent_strips(self, pixel, pixels_per_plane = pixels_per_plane, n = 2):
        '''
        Returns n adjacent strips within the same plane from a given strip (pixel). 
        N is determined by the spline width and can be adjusted as necessary. 
        Assumed pixel pitch of 200 microns and a spline width of 600 microns.
        '''
        return [(pixel + x) for x in range(-n, n+1) if np.ceil(float((pixel + x)/pixels_per_plane)) == np.ceil(float(pixel/pixels_per_plane))]
        


    def build_spline(self):
        '''
            Builds a cubic spline to interpolate the energy sharing
            distribution histogram. Has a default spline if no other data is given.
        '''

        with uproot.open(self.path_to_spline) as f:

            h = f['pmax_histogram'].to_hist()


            from scipy.interpolate import UnivariateSpline, CubicSpline

            centers = h.axes.centers[0]
            
            amps = h.values()/np.amax(h.values())

            spline = CubicSpline(centers, amps, extrapolate = False)

        return spline


    def compute_energy_sharing(self, atar, spline, entry, pixel_pitch = 200):

        '''
        Returns a digitized form atar pixel data in dictionary form
        '''

        output = {
            'pixel_pdg':[],
            'pixel_edep':[],
            'pixel_time':[],
            'pixel_hits':[]
        }

        

        data = atar.arrays(['pixel_pdg', 'pixel_edep', 'pixel_time', 'pixel_hits'])
        pdg = data['pixel_pdg'][entry]
        edep = data['pixel_edep'][entry]
        time = data['pixel_time'][entry]
        hits = data['pixel_hits'][entry]

        for v in ['pixel_edep'].array():  
            
            for i, hit in enumerate(hits):
                #Loop over all hits and calculate energy sharing

                these_energies = []

                places = self.get_adjacent_strips(hit, n = 2)
                        
                for pos in places:
                    energy_position = (hit - pos) * pixel_pitch
                    energyi = spline(energy_position) * edep[i]
                    these_energies.append(energyi)
                
                #Normalize the energies to conserve energy
                these_energies = np.array(these_energies) * (edep[i]/np.sum(these_energies))

                output['pixel_edep'] += list(these_energies)
                output['pixel_hits'] += self.get_adjacent_strips(hits)
                output['pixel_time'] += list(np.full_like(these_energies, time[i]))
                output['pixel_pdg'] += list(np.full_like(these_energies, pdg[i], dtype=int))



    


