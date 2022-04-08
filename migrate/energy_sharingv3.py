import uproot
import matplotlib.pyplot as plt
import hist as Hist
import numpy as np
import os 
import shutil


class digitizeatar:

    '''
        This class computes energy sharing in the Atar. It has methods to build a spline
        and return atar data with energysharing computed. You must initialized the class
        with a path to a histogram to build the spline. This spline, once built, can be
        passed into the compute energysharing method to return updated data. 
    '''

    #Here is a given path to a spline './BNL_Signal_Response.root'

    def __init__(self, path_to_spline, pixel_pitch = 200, pixels_per_plane = 100):
        
        '''
        Energysharing spline path must be provided by the user. There are given values for pixel pitch and pixels per plane (strips per plane) 
        based on current simulation parameters. These can be changed. 
        '''
        self.path_to_spline =  path_to_spline
        self.pixel_pitch = pixel_pitch
        self.pixels_per_plane = pixels_per_plane



    #For the paramaters in the atar simulation, the number of slits accessible on either side is two.

    def get_adjacent_strips(self, pixel, n = 2):
        '''
        Returns n adjacent strips within the same plane from a given strip (pixel). 
        N is determined by the spline width and can be adjusted as necessary. 
        Assumed pixel pitch of 200 microns and a spline width of 600 microns.
        '''
        return [(pixel + x) for x in range(-n, n+1) if int(np.ceil(float((pixel + x)/self.pixels_per_plane))) == int(np.ceil(float(pixel)/self.pixels_per_plane))]
        


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


    def compute_energy_sharing(self, root_path, spline, entry = "All"):

        '''
        Returns a root file in the current working directory of the form 'digitized <root_path>'. 
        Must be called with a path to the root file you want to digitize and a spline for energy sharing. The new
        root file has a digitized atar branch with energy sharing data computed. If no event is specified, it will compute energy sharing for all entries. 
        TODO: delete old atar branch and copy over relevant data to new file. 
        '''

        dest = os.getcwd() + '/digitized' + os.path.basename(root_path)

        shutil.copy(root_path, dest)


        
        #TTree of updated data

        tree = {
            'pixel_pdg':[],
            'pixel_edep':[],
            'pixel_time':[],
            'pixel_hits':[]
        }

        
        with uproot.open(root_path + ':atar') as infile:

            data = infile.arrays(['pixel_pdg', 'pixel_edep', 'pixel_time', 'pixel_hits'])
            pdg = data['pixel_pdg']
            edep = data['pixel_edep']
            time = data['pixel_time']
            hits = data['pixel_hits']
        
            if entry != 'All':
                pdg = [data['pixel_pdg'][entry]]
                edep = [data['pixel_edep'][entry]]
                time = [data['pixel_time'][entry]]
                hits = [data['pixel_hits'][entry]]
    
            for ent, event in enumerate(hits):
                #Go over all events 
                edep_list = []
                hits_list = []
                time_list = []
                pdg_list = []
                    
                for i, hiti in enumerate(event):
                    #Loop over all hits for each event and calculate energy sharing

                    these_energies = []
                    
                    for j, pos in enumerate(self.get_adjacent_strips(hiti)):
                        energy_position = (pos - hiti) * self.pixel_pitch
                        energyi = spline(energy_position)*edep[ent][i]
                        these_energies.append(energyi)
                    
                    #Normalize the energies to conserve energy
                    these_energies = np.array(these_energies) * (edep[ent][i]/np.sum(these_energies))

                    # print("total Energy", np.sum(these_energies), edep[i])

                    edep_list += list(these_energies)
                    hits_list += self.get_adjacent_strips(hiti)
                    time_list += list(np.full_like(these_energies, time[ent][i]))
                    pdg_list += list(np.full_like(these_energies, pdg[ent][i], dtype=int))
                
                #Add the updated data to the new tree branches

                tree['pixel_pdg'].append(pdg_list)
                tree['pixel_edep'].append(edep_list)
                tree['pixel_time'].append(time_list)
                tree['pixel_hits'].append(hits_list)

        # Add on TTree
        with uproot.update(dest) as f:
            f['digitized_atar'] = tree
            

    


    


