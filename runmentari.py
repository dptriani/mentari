import numpy as np
import h5py
import mentari_v2 as mtr

directory_dusty = "../output/test-bulge4/"
snap_limit = 63
#BoxSize = ((62.5**3) * (6/8))**(1/3) #mini-millennium
firstfile = 0
lastfile = 0
Hubble_h = 0.73

mass_dusty, metals_dusty = mtr.build_mass_and_metallicity_history(1, directory_dusty, firstfile, lastfile, snap_limit)
dust, gas_metals, gas, rad  = mtr.build_dust_history(1, directory_dusty, firstfile, lastfile, snap_limit)

#Compute attenuation parameters
w = np.where((mass_dusty[:,snap_limit] > 0) & (dust[:,snap_limit] > 0))[0]
Mass = mass_dusty[w] / Hubble_h 
Metals = metals_dusty[w]

Dust = dust[w,snap_limit] / Hubble_h
Gas = gas[w,snap_limit] / Hubble_h
Rad = rad[w,snap_limit] / Hubble_h
eta_BC = [-0.7] * len(Dust)
eta_ISM_v2 = [-1.3] * len(Dust)

prescription = 0 #0 for Lagos+ 19; 1 for Somerville+ 12
tau_BC, eta_BC, tau_ISM, eta_ISM = mtr.compute_attenuation_parameters (prescription, Dust, Gas, Rad)

#Generate SED
Age = np.asarray([0.0124, 0.0246, 0.0491, 0.1037, 0.1871, 0.2120, 0.2399, 0.2709, 0.3054, 0.3438, 0.3864, 0.4335, 0.4856, 0.5430, 0.6062, 0.6756, 0.7517, 0.8349, 0.9259, 1.0249, 1.1327, 1.2496, 1.3763, 1.5131, 1.6606, 1.8192, 1.9895, 2.1717, 2.3662, 2.5734, 2.7934, 3.0265, 3.2726, 3.5318, 3.8038, 4.0886, 4.3856, 4.6944, 5.0144, 5.3488, 5.6849, 6.0337, 6.3901, 6.7531, 7.1215, 7.4940, 7.8694, 8.2464, 8.6238, 9.0004, 9.3750, 9.7463, 10.1133, 10.4750, 10.8303, 11.1783, 11.5181, 11.8490, 12.1702, 12.4811, 12.7810, 13.0695, 13.3459, 13.6098])
time_BC = 10**7
SSP = 0 #0 for BC03
wavelength, spectra, spectra_dusty = mtr.generate_SED(0, Age, Mass, Metals, 
             tau_BC, tau_ISM, eta_BC, eta_ISM, time_BC)

wavelength_IR, spectra_IR = mtr.add_IR_Dale(wavelength, spectra, spectra_dusty)
wave_FIR, spec_FIR = mtr.compute_IR_SUNRISE (Dust, wavelength, spectra, spectra_dusty)
wavelength_IR_b, spectra_IR_b = mtr.combine_Dale_SUNRISE(wave_FIR, spec_FIR, wavelength_IR, spectra_IR)

with h5py.File('mentari_output.h5', 'a') as f:
    f.create_dataset('StellarMass', data=Mass, maxshape=(None,snap_limit+1))

    f.create_dataset('Metallicity', data=Metals, maxshape=(None,snap_limit+1))
    f.create_dataset('DustMass', data=Dust, maxshape=(None,))
    f.create_dataset('GasMass', data=Gas, maxshape=(None,))
    f.create_dataset('Radius', data=Rad, maxshape=(None,))
    f.create_dataset('Wavelength_UVIR', data=wavelength_IR, maxshape=(None,))
    f.create_dataset('Spectra_UVIR', data=spectra_IR, maxshape=(None,len(wavelength_IR)))
    f.create_dataset('Wavelength_stellar', data=wavelength, maxshape=(None,))
    f.create_dataset('Spectra_stellar', data=spectra, maxshape=(None,len(wavelength)))
    f.create_dataset('Wavelength_SUNRISE', data=wavelength_IR_b, maxshape=(None,))
    f.create_dataset('Spectra_SUNRISE', data=spectra_IR_b, maxshape=(None,len(wavelength_IR_b)))
'''
with h5py.File('mentari_output.h5', 'a') as f:
    f['StellarMass'].resize((f['StellarMass'].shape[0] + Mass.shape[0]), axis=0)
    f['StellarMass'][-Mass.shape[0]:] = Mass
    f['Metallicity'].resize((f['Metallicity'].shape[0] + Metals.shape[0]), axis=0)
    f['Metallicity'][-Metals.shape[0]:] = Metals
    f['DustMass'].resize((f['DustMass'].shape[0] + Dust.shape[0]), axis=0)
    f['DustMass'][-Dust.shape[0]:] = Dust
    f['GasMass'].resize((f['GasMass'].shape[0] + Gas.shape[0]), axis=0)
    f['GasMass'][-Gas.shape[0]:] = Gas
    f['Radius'].resize((f['Radius'].shape[0] + Rad.shape[0]), axis=0)
    f['Radius'][-Rad.shape[0]:] = Rad
    f['Spectra_UVIR'].resize((f['Spectra_UVIR'].shape[0] + spectra_IR.shape[0]), axis=0)
    f['Spectra_UVIR'][-spectra_IR.shape[0]:] = spectra_IR
    f['Spectra_stellar'].resize((f['Spectra_stellar'].shape[0] + spectra.shape[0]), axis=0)
    f['Spectra_stellar'][-spectra.shape[0]:] = spectra
    f['Spectra_SUNRISE'].resize((f['Spectra_SUNRISE'].shape[0] + spectra_IR_b.shape[0]), axis=0)
    f['Spectra_SUNRISE'][-spectra_IR_b.shape[0]:] =spectra_IR_b
'''
