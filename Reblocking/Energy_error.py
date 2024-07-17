import numpy as np
import pandas as pd
import pyblock 
import sys
import matplotlib.pyplot as plt

stdoutOrigin=sys.stdout 
sys.stdout = open("Reblocking_log.txt", "w") #print the output in a file

data = np.loadtxt('energie.txt')
print('Steps in energie.txt = ', len(data) )

Tequil = 15 #au
ene = data[ data[:,1]>Tequil ,2]
print('Steps di statistica = ', len(ene) )

print( 'Media Et = ', ene.mean() )
print( 'std Et = ', ene.std() )

reblock_data = pyblock.blocking.reblock(ene)
opt = pyblock.blocking.find_optimal_block(len(ene), reblock_data)

print('Optimal block =', opt )
print(reblock_data[opt[0]])  #stampa solo l'iterazione di reblocking associata alla lunghezza di blocco (cioè circa len(ene)/2^opt)

pds_ene = pd.Series(ene)
(data_length, reblock_data, covariance) = pyblock.pd_utils.reblock(pds_ene) #Note some reblocking iterations discard a data point if there were an odd number of data points in the previous iteration.
print(data_length)
print(reblock_data)

plt.rcParams['figure.dpi'] = 300
fig = pyblock.plot.plot_reblocking(reblock_data)
ax = fig.gca()      #per estrarre gli assi dalla figura
ax.set_ylim(0.00006, 0.00013)
ax.set_xlabel('Reblock iteration', fontsize=14)
ax.set_ylabel(r'Standard error $\sigma_{\langle E_T \rangle}$', fontsize=14)
ax.tick_params(axis='both', labelsize=11)



fig.savefig("fig.png", dpi=300)

sys.stdout.close()
sys.stdout=stdoutOrigin