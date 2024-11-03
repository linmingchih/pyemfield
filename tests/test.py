import sys, os
sys.path.insert(0, r'D:\OneDrive - ANSYS, Inc\GitHub\pyemfield')

from pyemfield import hfss_design, get_ffds, Beam

#%%

hd = hfss_design()

folder = r'D:\OneDrive - ANSYS, Inc\GitHub\pyemfield\tests\ffds'
hd.export_ffds(folder, 'Setup1 : Sweep', '30.0GHz')

#%%
folder = r'D:\OneDrive - ANSYS, Inc\GitHub\pyemfield\tests\ffds'
ffds = get_ffds(folder)

ffds.keys()

x = {j:(1,0) for i, j in ffds.items()}
b1 = Beam(x)

b1.ffd_excitation
b2 = b1.optimize_gain(60, 60)
b2.plot_realized_gain_contour()
b2.ffd_excitation

hd.update_excitation(b2.ffd_excitation)

