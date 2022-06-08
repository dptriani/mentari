
def save_mag(filter_list, input_filename, output_filename, z):
    import h5py
    import mentari_v2 as mtr
    import os
    import numpy as np

    file_input = input_filename
    print(file_input)
    with h5py.File(file_input, 'r') as f:
        
        wavelength_m1 = np.array(f['Wavelength_m1'])
        spectra_m1 = np.array(f['Spectra_m1'])
        wavelength_m2 = np.array(f['Wavelength_m2'])
        spectra_m2 = np.array(f['Spectra_m2'])
        wavelength_m3 = np.array(f['Wavelength_m3'])
        spectra_m3 = np.array(f['Spectra_m3'])
        wavelength_m4 = np.array(f['Wavelength_m4'])
        spectra_m4 = np.array(f['Spectra_m4'])
        wavelength_s = np.array(f['Wavelength_stellar'])
        spectra_s = np.array(f['Spectra_stellar'])
        
    
    file_output = output_filename
    print(file_output)
    if os.path.isfile(file_output) == 0:
        mab1 = mtr.compute_mab(wavelength_m1, spectra_m1, filter_list, z)
        mab2 = mtr.compute_mab(wavelength_m2, spectra_m2, filter_list, z)
        mab3 = mtr.compute_mab(wavelength_m3, spectra_m3, filter_list, z)
        mab4 = mtr.compute_mab(wavelength_m4, spectra_m4, filter_list, z)
        mabs = mtr.compute_mab(wavelength_s, spectra_s, filter_list, z)

        with h5py.File(file_output, 'w') as f:
            f.create_dataset('default', data = mab1)
            f.create_dataset('SUNRISE', data = mab2)
            f.create_dataset('Somerville', data = mab3)
            f.create_dataset('CF00', data = mab4)
            f.create_dataset('unattenuated', data=mabs)

    return
#-----------------------------------------------------------------------------------	

def distributed_processing(filter_list, input_file, output_file, z):
    
    import sys
    import os
    import time
    
    tstart = time.perf_counter()
    rank = 0
    ntasks = 1
    comm = None
    
    try:
        from mpi4py import MPI
        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()
        ntasks = comm.Get_size()
    except ImportError:
        pass
    
    sys.stdout.flush()
    nfiles = len(input_file)
    if nfiles < ntasks:
        print(f"[Rank={rank}]: Nfiles = {nfiles} < total tasks = {ntasks}. "
            "Some tasks will not have any work assigned (and will be idle)")
    
    if rank == 0:
        print(f"[Rank={rank}]: Running nfiles = {nfiles} over "\
              f"ntasks = {ntasks}...")
        
    for filenum in range(rank, nfiles, ntasks):
        save_mag( filter_list, input_file[filenum], output_file[filenum], z)
        
    # The barrier is only essential so that the total time printed
    # out on rank==0 is correct.
    if comm:
        comm.Barrier()

    if rank == 0:
        t1 = time.perf_counter()
        print(f"[Rank={rank}]: Running nfiles = {nfiles} over "\
              f"ntasks = {ntasks}...done. Time taken = {t1-tstart:0.3f} seconds")
        
    return True


if __name__ == '__main__':
    
#    filter_list = ['GALEX_FUV', 'GALEX_NUV', 'TwoMass_Ks', 'VIRCAM_K', 'Sdss_u', 
#                  'Sdss_g', 'Sdss_r', 'Sdss_i', 'Sdss_z', 'IRAC_1', 'IRAC_2',
#                  'IRAC_3', 'IRAC_4', 'MIPS_24um', 'PACS_70um', 'PACS_160um',
#                   'SPIRE_250um', 'SPIRE_350um', 'SPIRE_500um', 'SCUBA_850WB']

    filter_list = ['GALEX_FUV', 'TwoMass_Ks', 'VIRCAM_K','IRAC_4','SPIRE_250um']
    z = 0 

    dirname_out = 'output_mag/'
    dirname_in = 'output_v2/'
    z_in = [0.509, 1.078, 2.07, 3.06] 
    z_out = [0.5, 1.0, 2.0, 3.0]
    #z_in = [0.509, 3.06]
    #z_out = [0.5, 1.0]
    #z_in = [1.078, 2.07, 3.06]
    #z_out = [1.0, 2.0, 3.0]
    firstfile = 0
    lastfile = 1
    name_input = 'mentari_output_z'
    name_output = 'mentari_mag_'
    ext = '.hdf5'
    file_input = []
    file_output = []
    
    for j in range(len(z_in)):    
        for i in range(firstfile, lastfile+1):
            file_input.append(dirname_in + name_input + str(z_in[j]) + '-' + str(i) + ext)
            file_output.append(dirname_out + name_output + str(z_out[j]) + '_' + str(i) + ext)
#            print(file_input)
#    for i in range(firstfile, lastfile+1):
#        file_input.append(dirname_in + name_input + '0.0-'+ str(i) + ext)
#        file_output.append(dirname_out + name_output + '0.0_' + str(i) + ext)
        
    distributed_processing(filter_list, file_input, file_output,z)
