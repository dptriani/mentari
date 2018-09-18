import numpy as np
import matplotlib as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import SED
BoxSize = 62.5
Hubble_h = 0.73

num_of_file = 1
directory = '../millennium-full/' #change this to the output directory of sage
filename = 'model'
redshift = [127.000, 79.998, 50.000, 30.000, 19.916, 18.244, 16.725, 15.343, 14.086, 12.941, 11.897, 10.944, 10.073, 9.278, 8.550, 7.883, 7.272, 6.712, 6.197, 5.724, 5.289, 4.888, 4.520, 4.179, 3.866, 3.576, 3.308, 3.060, 2.831, 2.619, 2.422, 2.239, 2.070, 1.913, 1.766, 1.630, 1.504, 1.386, 1.276, 1.173, 1.078, 0.989, 0.905, 0.828, 0.755, 0.687, 0.624, 0.564, 0.509, 0.457, 0.408, 0.362, 0.320, 0.280, 0.242, 0.208, 0.175, 0.144, 0.116, 0.089, 0.064, 0.041, 0.020, 0.000]

print('reading redshift list')
G = SED.read_redshift_list(redshift, num_of_file, directory, filename)
print('finish reading, constructing histories')
mass, metal = SED.mass_metal_history(redshift, num_of_file,directory, filename)
'''
#take the recycle fraction in account
rec_frac = 0.43
c_mass = mass[0] * (1. - rec_frac) / Hubble_h #final computed mass

LastSnap = len(redshift) - 1
s_mass = G[LastSnap].StellarMass * 1.e10 / Hubble_h

#Plot stellar mass
binwidth = 0.3
ax0 = plt.subplot2grid((1,1), (0,0))
divider = make_axes_locatable(ax0) 
ax1 = divider.append_axes("bottom", size="50%", pad=0.3)

# Baldry+ 2008 modified data used for the MCMC fitting
Baldry = np.array([
            [7.05, 1.3531e-01, 6.0741e-02],
            [7.15, 1.3474e-01, 6.0109e-02],
            [7.25, 2.0971e-01, 7.7965e-02],
            [7.35, 1.7161e-01, 3.1841e-02],
            [7.45, 2.1648e-01, 5.7832e-02],
            [7.55, 2.1645e-01, 3.9988e-02],
            [7.65, 2.0837e-01, 4.8713e-02],
            [7.75, 2.0402e-01, 7.0061e-02],
            [7.85, 1.5536e-01, 3.9182e-02],
            [7.95, 1.5232e-01, 2.6824e-02],
            [8.05, 1.5067e-01, 4.8824e-02],
            [8.15, 1.3032e-01, 2.1892e-02],
            [8.25, 1.2545e-01, 3.5526e-02],
            [8.35, 9.8472e-02, 2.7181e-02],
            [8.45, 8.7194e-02, 2.8345e-02],
            [8.55, 7.0758e-02, 2.0808e-02],
            [8.65, 5.8190e-02, 1.3359e-02],
            [8.75, 5.6057e-02, 1.3512e-02],
            [8.85, 5.1380e-02, 1.2815e-02],
            [8.95, 4.4206e-02, 9.6866e-03],
            [9.05, 4.1149e-02, 1.0169e-02],
            [9.15, 3.4959e-02, 6.7898e-03],
            [9.25, 3.3111e-02, 8.3704e-03],
            [9.35, 3.0138e-02, 4.7741e-03],
            [9.45, 2.6692e-02, 5.5029e-03],
            [9.55, 2.4656e-02, 4.4359e-03],
            [9.65, 2.2885e-02, 3.7915e-03],
            [9.75, 2.1849e-02, 3.9812e-03],
            [9.85, 2.0383e-02, 3.2930e-03],
            [9.95, 1.9929e-02, 2.9370e-03],
            [10.05, 1.8865e-02, 2.4624e-03],
            [10.15, 1.8136e-02, 2.5208e-03],
            [10.25, 1.7657e-02, 2.4217e-03],
            [10.35, 1.6616e-02, 2.2784e-03],
            [10.45, 1.6114e-02, 2.1783e-03],
            [10.55, 1.4366e-02, 1.8819e-03],
            [10.65, 1.2588e-02, 1.8249e-03],
            [10.75, 1.1372e-02, 1.4436e-03],
            [10.85, 9.1213e-03, 1.5816e-03],
            [10.95, 6.1125e-03, 9.6735e-04],
            [11.05, 4.3923e-03, 9.6254e-04],
            [11.15, 2.5463e-03, 5.0038e-04],
            [11.25, 1.4298e-03, 4.2816e-04],
            [11.35, 6.4867e-04, 1.6439e-04],
            [11.45, 2.8294e-04, 9.9799e-05],
            [11.55, 1.0617e-04, 4.9085e-05],
            [11.65, 3.2702e-05, 2.4546e-05],
            [11.75, 1.2571e-05, 1.2571e-05],
            [11.85, 8.4589e-06, 8.4589e-06],
            [11.95, 7.4764e-06, 7.4764e-06],
            ], dtype=np.float32)
        
        # Finally plot the data
        # plt.errorbar(
        #     Baldry[:, 0],
        #     Baldry[:, 1],
        #     yerr=Baldry[:, 2],
        #     color='g',
        #     linestyle=':',
        #     lw = 1.5,
        #     label='Baldry et al. 2008',
        #     )
Baldry_xval = np.log10(10 ** Baldry[:, 0]  /Hubble_h/Hubble_h)
Baldry_xval = Baldry_xval - 0.26  # convert back to Chabrier IMF
Baldry_yvalU = (Baldry[:, 1]+Baldry[:, 2]) * Hubble_h*Hubble_h*Hubble_h
Baldry_yvalL = (Baldry[:, 1]-Baldry[:, 2]) * Hubble_h*Hubble_h*Hubble_h

ax0.fill_between(Baldry_xval, Baldry_yvalU, Baldry_yvalL, 
    facecolor='purple', alpha=0.25, label='Baldry et al. 2008 (z=0.1)')

mi = 8
ma = 12
NB = int((ma - mi) / binwidth)

(counts, binedges) = np.histogram(np.log10(c_mass), range=(mi, ma), bins=NB)
(Counts, Binedges) = np.histogram(np.log10(s_mass), range=(mi, ma), bins=NB)


# Set the x-axis values to be the centre of the bins
xaxeshisto = binedges[:-1] + 0.5 * binwidth
Xaxeshisto = Binedges[:-1] + 0.5 * binwidth

# Overplot the model histograms
ax0.plot(xaxeshisto, counts / (BoxSize/Hubble_h)**3 / binwidth, c ='r', label='computed mass')
ax0.plot(Xaxeshisto, Counts / (BoxSize/Hubble_h)**3 / binwidth, c='b', label='mass from SAGE')
ax0.set_yscale('log', nonposy='clip')
ax0.set_xlim(8.0, 12.0)


# Set the x-axis minor ticks
ax0.xaxis.set_minor_locator(plt.MultipleLocator(0.1))
ax0.set_ylabel(r'$\phi\ (\mathrm{Mpc}^{-3}\ \mathrm{dex}^{-1}$)')  # Set the y...

ax1.scatter(np.log10(s_mass), (np.log10(c_mass) - np.log10(s_mass)), marker='.', s=4, linewidths=0.5)
ax1.set_ylabel(r'$\Delta \log_{10}\ M_{\mathrm{stars}}\ (M_{\odot})$')
ax1.set_xlabel(r'$\log_{10} M_{\mathrm{stars}}\ SAGE (M_{\odot})$')
ax1.set_ylim(-0.1, 0.1)
ax1.set_xlim(8, 12)

leg = ax0.legend(loc=0, numpoints=1, labelspacing=0.1)
leg.draw_frame(False)  # Don't want a box frame
for t in leg.get_texts():  # Reduce the size of the text
    t.set_fontsize('medium')
    
plt.tight_layout()
plt.show()
'''
