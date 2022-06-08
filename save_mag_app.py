
def save_mag(filter_list, input_filename, output_filename, z):
    import h5py
    import mentari_v2 as mtr
    import os
    import numpy as np

    file_input = input_filename
    print(file_input)
    with h5py.File(file_input, 'r') as f:
        mass = list(f['StellarMass'])
        wavelength_m1 = np.array(f['Wavelength_m1'])
        spectra_m1 = np.array(f['Spectra_m1'])
        wavelength_s = np.array(f['Wavelength_stellar'])
        spectra_s = np.array(f['Spectra_stellar'])
        
    
    file_output = output_filename
    print(file_output)
    if os.path.isfile(file_output) == 0:
        mab_app = mtr.compute_mab(wavelength_m1, spectra_m1, filter_list, z)
        mab_abs = mtr.compute_mab(wavelength_m1, spectra_m1, filter_list, 0)

        with h5py.File(file_output, 'w') as f:
            f.create_dataset('stellarmass', data = mass)
            f.create_dataset('apparent', data = mab_app)
            f.create_dataset('absolut', data = mab_abs)
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

#    filter_list = ['GALEX_FUV', 'TwoMass_Ks', 'VIRCAM_K','IRAC_4','SPIRE_250um']
    filter_list = ['WISE_W1'] 
    dirname_out = 'output_app/'
    dirname_in = 'output_app/'
    z_in = [0.755, 0.624, 0.564, 0.509]
    z_out = [0.755, 0.624, 0.564, 0.509]
    firstfile = 0
    lastfile = 31
    name_input = 'mentari_output_z'
    name_output = 'mag_'
    ext = '.hdf5'
    
    for j in range(len(z_in)): 
        file_input = []
        file_output = []
        for i in range(firstfile, lastfile+1):
            file_input.append(dirname_in + name_input + str(z_in[j]) + '-' + str(i) + ext)
            file_output.append(dirname_out + name_output + str(z_out[j]) + '_' + str(i) + ext)
        distributed_processing(filter_list, file_input, file_output,z_in[j])
#            print(file_input)
#    for i in range(firstfile, lastfile+1):
#        file_input.append(dirname_in + name_input + '0.0-'+ str(i) + ext)
#        file_output.append(dirname_out + name_output + '0.0_' + str(i) + ext)
        
       
