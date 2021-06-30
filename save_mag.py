import h5py
import mentari_v2 as mtr
import os
import numpy as np

filter_list = ['GALEX_FUV', 'GALEX_NUV', 'TwoMass_Ks', 'VIRCAM_K', 'Sdss_u', 
              'Sdss_g', 'Sdss_r', 'Sdss_i', 'Sdss_z', 'IRAC_1', 'IRAC_2',
              'IRAC_3', 'IRAC_4', 'MIPS_24um', 'PACS_70um', 'PACS_160um',
               'SPIRE_250um', 'SPIRE_350um', 'SPIRE_500um', 'SCUBA_850WB']

#filter_list = ['GALEX_FUV', 'TwoMass_Ks', 'VIRCAM_K','IRAC_4','SPIRE_250um']
z = 0 

dirname_out = 'output/'
dirname_in = 'output/'
name_input = 'mentari_output_z0.0-'
name_output = 'mentari_mag_z0_all_0'
ext = '.hdf5'
firstfile = 0
lastfile = 3

for i in range(firstfile, lastfile+1):
    file_input = dirname_in + name_input + str(i) + ext
    with h5py.File(file_input, 'r') as f:
        
        wavelength_m1 = np.array(f['Wavelength_m1'])
        spectra_m1 = np.array(f['Spectra_m1'])
        wavelength_m2 = np.array(f['Wavelength_m2'])
        spectra_m2 = np.array(f['Spectra_m2'])
        wavelength_m3 = np.array(f['Wavelength_m3'])
        spectra_m3 = np.array(f['Spectra_m3'])
        wavelength_m4 = np.array(f['Wavelength_m4'])
        spectra_m4 = np.array(f['Spectra_m4'])
        
    file_output = dirname_out + name_output + str(i) + ext
    print(file_output)
    if os.path.isfile(file_output) == 0:
        mab1 = mtr.compute_mab(wavelength_m1, spectra_m1, filter_list, z)
        mab2 = mtr.compute_mab(wavelength_m2, spectra_m2, filter_list, z)
        mab3 = mtr.compute_mab(wavelength_m3, spectra_m3, filter_list, z)
        mab4 = mtr.compute_mab(wavelength_m4, spectra_m4, filter_list, z)

        with h5py.File(file_output, 'w') as f:
            f.create_dataset('default', data = mab1)
            f.create_dataset('SUNRISE', data = mab2)
            f.create_dataset('Somerville', data = mab3)
            f.create_dataset('CF00', data = mab4)

