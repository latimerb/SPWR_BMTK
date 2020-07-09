import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import pdb
from mpl_toolkits.mplot3d import Axes3D
import h5py
import matplotlib.pyplot as plt
from scipy.signal import welch

# Load data
config_file = "simulation_config.json"
lfp_file = "./output/ecp.h5"
mem_pot_file = './output/v_report.h5'
raster_file = './output/spikes.h5'
conns_file = './network/SPWR_biophysical_SPWR_biophysical_edges.h5'

pyr = [0,21599]
axo = [21600,22139]
bask = [22140,26999]


# load 
#f = h5py.File(mem_pot_file,'r')
#mem_potential = f['report']['SPWR_biophysical']['data']

#f = h5py.File(lfp_file,'r')
#lfp = list(f['ecp']['data'])
#lfp_arr = np.asarray(lfp)

f = h5py.File(raster_file,'r')
gids = f['spikes']['SPWR_biophysical']['node_ids']
timestamps = f['spikes']['SPWR_biophysical']['timestamps']

f = h5py.File(conns_file,'r')
src = f['edges']['SPWR_biophysical_SPWR_biophysical']['source_node_id']
tgt = f['edges']['SPWR_biophysical_SPWR_biophysical']['target_node_id']
wgt = f['edges']['SPWR_biophysical_SPWR_biophysical']['0']['syn_weight']

id, pyr2baskconns = np.unique(src[(src[:]<=pyr[1]) & (tgt[:]>=bask[0]) \
			& (tgt[:]<=bask[1])],return_counts=True)


pyr2baskwgts = wgt[(src[:]<=pyr[1]) & (tgt[:]>=bask[0]) \
			& (tgt[:]<=bask[1])]


plt.figure()
plt.subplot(2,1,1)
plt.hist(pyr2baskconns,20)
plt.title('pyr2baskconns')
plt.subplot(2,1,2)
plt.hist(pyr2baskwgts,20)
plt.title('pyr2baskwgts')


# Plot data

#plt.figure()
#plt.plot(mem_potential[:,0])
#plt.plot(mem_potential[:,1])
#plt.plot(mem_potential[:,2])
#plt.plot(mem_potential[:,3])

#plt.figure()
#plt.plot(lfp_arr[:,0])
#plt.plot(lfp_arr[:,1])
#plt.plot(lfp_arr[:,2])


plt.figure()
plt.plot(timestamps,gids,'.')

#plt.figure()
#x = lfp_arr[:,0]
#fs = 10000
#f, Pxx_den = welch(x, fs, nperseg=2000)
#plt.semilogy(f, Pxx_den)
#plt.xlim(0,500)

plt.show()

