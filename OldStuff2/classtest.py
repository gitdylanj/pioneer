
import matplotlib.pyplot as plt
import numpy as np
import uproot  
from energy_sharingv2 import digitizeatar

path = "./BNL_Signal_Response.root"

f = uproot.open(path)

h = f['pmax_histogram'].to_hist()

x = digitizeatar(path)

out = x.build_spline()
xs = np.linspace(h.axes[0].centers[0], h.axes[0].centers[-1], 1000)



ys = x.build_spline()(xs)

fig, ax = plt.subplots()

plt.plot(xs, ys)

strips = x.get_adjacent_strips(200)
print(strips)

