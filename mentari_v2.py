from __future__ import print_function
from random import sample, seed
import copy
from astropy.cosmology import FlatLambdaCDM
cosmo = FlatLambdaCDM(H0=73, Om0=0.25)
from os.path import dirname, abspath, join as pjoin
import re, os
import numpy as np

def galdtype_dusty(align=True):
     
    '''Define the data-type for the public version of Dusty-SAGE'''
    
    Galdesc_full = [
        ('SnapNum'                      , np.int32),
        ('Type'                         , np.int32),
        ('GalaxyIndex'                  , np.int64),
        ('CentralGalaxyIndex'           , np.int64),
        ('SAGEHaloIndex'                , np.int32),
        ('SAGETreeIndex'                , np.int32),
        ('SimulationHaloIndex'          , np.int64),
        ('mergeType'                    , np.int32),
        ('mergeIntoID'                  , np.int32),
        ('mergeIntoSnapNum'             , np.int32),
        ('dT'                           , np.float32),
        ('Pos'                          , (np.float32, 3)),
        ('Vel'                          , (np.float32, 3)),
        ('Spin'                         , (np.float32, 3)),
        ('Len'                          , np.int32),
        ('Mvir'                         , np.float32),
        ('CentralMvir'                  , np.float32),
        ('Rvir'                         , np.float32),
        ('Vvir'                         , np.float32),
        ('Vmax'                         , np.float32),
        ('VelDisp'                      , np.float32),
        ('ColdGas'                      , np.float32),
        ('f_H2'                         , np.float32),
        ('f_HI'                         , np.float32),
        ('cf'                           , np.float32),
        ('Zp'                           , np.float32),
        ('Pressure'                     , np.float32),
        ('StellarMass'                  , np.float32),
        ('BulgeMass'                    , np.float32),
        ('BulgeInstability'             , np.float32),
        ('HotGas'                       , np.float32),
        ('EjectedMass'                  , np.float32),
        ('BlackHoleMass'                , np.float32),
        ('IntraClusterStars'            , np.float32),
        ('MetalsColdGas'                , np.float32),
        ('MetalsStellarMass'            , np.float32),
        ('MetalsBulgeMass'              , np.float32),
        ('MetalsHotGas'                 , np.float32),
        ('MetalsEjectedMass'            , np.float32),
        ('MetalsIntraClusterStars'      , np.float32),
        ('ColdDust'                     , np.float32),
        ('HotDust'                      , np.float32),
        ('EjectedDust'                     , np.float32),
        ('SfrDisk'                      , np.float32),
        ('SfrBulge'                     , np.float32),
        ('SfrDiskZ'                     , np.float32),
        ('SfrBulgeZ'                    , np.float32),
        ('SfrDiskDTG'                     , np.float32),
        ('SfrBulgeDTG'                    , np.float32),
        ('dustdotform'                  , np.float32),
        ('dustdotgrowth'                    , np.float32),
        ('dustdotdestruct'                    , np.float32),
        ('DiskRadius'                   , np.float32),
        ('Cooling'                      , np.float32),
        ('Heating'                      , np.float32),
        ('QuasarModeBHaccretionMass'    , np.float32),
        ('TimeOfLastMajorMerger'         , np.float32),
        ('TimeOfLastMinorMerger'         , np.float32),
        ('OutflowRate'                  , np.float32),
        ('infallMvir'                   , np.float32),
        ('infallVvir'                   , np.float32),
        ('infallVmax'                   , np.float32)
        ]
    names = [Galdesc_full[i][0] for i in range(len(Galdesc_full))]
    formats = [Galdesc_full[i][1] for i in range(len(Galdesc_full))]
    if(align==True):
        Galdesc = np.dtype({'names':names, 'formats':formats}, align=True)
    else:
        Galdesc = np.dtype({'names':names, 'formats':formats})
    return Galdesc

#-----------------------------------------------------------------------------------	

def galdtype(align=True):
    
    '''Define the data-type for the public version of SAGE'''
    
    Galdesc_full = [
        ('SnapNum'                      , np.int32),
        ('Type'                         , np.int32),
        ('GalaxyIndex'                  , np.int64),
        ('CentralGalaxyIndex'           , np.int64),
        ('SAGEHaloIndex'                , np.int32),
        ('SAGETreeIndex'                , np.int32),
        ('SimulationHaloIndex'          , np.int64),
        ('mergeType'                    , np.int32),
        ('mergeIntoID'                  , np.int32),
        ('mergeIntoSnapNum'             , np.int32),
        ('dT'                           , np.float32),
        ('Pos'                          , (np.float32, 3)),
        ('Vel'                          , (np.float32, 3)),
        ('Spin'                         , (np.float32, 3)),
        ('Len'                          , np.int32),
        ('Mvir'                         , np.float32),
        ('CentralMvir'                  , np.float32),
        ('Rvir'                         , np.float32),
        ('Vvir'                         , np.float32),
        ('Vmax'                         , np.float32),
        ('VelDisp'                      , np.float32),
        ('ColdGas'                      , np.float32),
        ('StellarMass'                  , np.float32),
        ('BulgeMass'                    , np.float32),
        ('HotGas'                       , np.float32),
        ('EjectedMass'                  , np.float32),
        ('BlackHoleMass'                , np.float32),
        ('IntraClusterStars'            , np.float32),
        ('MetalsColdGas'                , np.float32),
        ('MetalsStellarMass'            , np.float32),
        ('MetalsBulgeMass'              , np.float32),
        ('MetalsHotGas'                 , np.float32),
        ('MetalsEjectedMass'            , np.float32),
        ('MetalsIntraClusterStars'      , np.float32),
        ('SfrDisk'                      , np.float32),
        ('SfrBulge'                     , np.float32),
        ('SfrDiskZ'                     , np.float32),
        ('SfrBulgeZ'                    , np.float32),
        ('DiskRadius'                   , np.float32),
        ('Cooling'                      , np.float32),
        ('Heating'                      , np.float32),
        ('QuasarModeBHaccretionMass'    , np.float32),
        ('TimeOfLastMajorMerger'         , np.float32),
        ('TimeOfLastMinorMerger'         , np.float32),
        ('OutflowRate'                  , np.float32),
        ('infallMvir'                   , np.float32),
        ('infallVvir'                   , np.float32),
        ('infallVmax'                   , np.float32)
        ]
    names = [Galdesc_full[i][0] for i in range(len(Galdesc_full))]
    formats = [Galdesc_full[i][1] for i in range(len(Galdesc_full))]
    if(align==True):
        Galdesc = np.dtype({'names':names, 'formats':formats}, align=True)
    else:
        Galdesc = np.dtype({'names':names, 'formats':formats})
    return Galdesc

#-----------------------------------------------------------------------------------	
def iterate_trees(SAM_option, directory, firstfile, lastfile):
    '''
    Iterating trees from the simulation output.
    Currently, it can only read trees from SAGE (Croton et al. 2006, 2016)
    and dusty-sage (Triani et al. 2020)
    
    Input:  - SAM option (int): (0) SAGE (1) Dusty-SAGE
            - path of the directory containing simulation output (string).
              Format of the simulation output: model_zX.XXX_Y
              X.XXX : redshift of the snapshot
              Y : file number
    
    Output: a tree, consist of properties of galaxies listed in galdtype_dusty() for Dusty SAGE or galdtype() for SAGE
    '''
    
    #define variables
    entries = [e for e in os.listdir(directory) 
               if os.path.isfile(os.path.join(directory, e))]
    entries = [e for e in entries if e.startswith('model_z')]
    redshift_strings = list(set([re.match(r'model_z(\d+\.?\d*)_\d+', e).group(1)
                                for e in entries]))
    #group_strings = list(set([re.match(r'model_z\d+\.?\d*_(\d+)', e).group(1)
    #                        for e in entries]))
    
    #group_strings.sort(key=lambda x: int(x))
    redshift_strings.sort(key=lambda x: float(x), reverse=True)
    
    if SAM_option == 0:
        Galdesc_false = galdtype(align=False)
        Galdesc=galdtype(align=True)
    elif SAM_option == 1:
        Galdesc_false = galdtype_dusty(align=False)
        Galdesc=galdtype_dusty(align=True)
    else:
        print("Choose a SAM: 0 - for SAGE, 1 - for Dusty-SAGE")
    
    #open files
    #for group in group_strings:
    for group in range(firstfile, lastfile+1):
        files = []
        for redshift in redshift_strings:
            fn = 'model_z%s_%s' % (redshift, group)
            files.append(open(os.path.join(directory, fn), 'rb'))
            
        n_trees = [np.fromfile(f, np.uint32, 1)[0] for f in files][0]
        n_gals = [np.fromfile(f, np.uint32, 1)[0] for f in files]
        chunk_sizes = [np.fromfile(f, np.uint32, n_trees) for f in files]
        tree_sizes = np.sum(chunk_sizes, axis=0)
        
        for ii in range(n_trees):
            tree_size = tree_sizes[ii]
            tree = np.empty(tree_size, dtype=Galdesc_false)
            offs=0
            for jj in range(len(chunk_sizes)):
                chunk_size = chunk_sizes[jj][ii]
                if chunk_size <= 0: continue

                data = np.fromfile(files[jj], Galdesc, chunk_size)

                for _v in data.dtype.names:
                    tree[_v][offs:offs+chunk_size] = data[_v]
                offs += chunk_size

            # First validate ID's.
            for f in ['Type', 'GalaxyIndex', 'CentralGalaxyIndex']:
                if min(tree[f]) < 0:
                    print("ERROR; min(tree[{0}]) = {1} should be non-zero"
                            .format(f, min(tree[f])))
                    raise ValueError()

            # Validate central galaxy index (unique id, generated by sage)
            ind = (np.where(tree['Type'] == 0))[0]
            if not bool(np.all(tree['GalaxyIndex'][ind] == tree['CentralGalaxyIndex'][ind])):
                print("tree[GalaxyIndex][ind] = {0}".format(tree['GalaxyIndex'][ind]))
                print("tree[CentralGalaxyIndex][ind] = {0}".format(tree['CentralGalaxyIndex'][ind]))

            assert bool(np.all(tree['GalaxyIndex'][ind] ==
                            tree['CentralGalaxyIndex'][ind])), \
                "Central Galaxy Index must equal Galaxy Index for centrals"
            
            yield tree
        
            
        for file in files:
            file.close()
            

#-----------------------------------------------------------------------------------	
            
def calculate_mass_and_metals(SAM_choice, tree, snap_limit):
    
    """Calculate mass history from Dusty-SAGE tree.
    In one fly, it will calculate the mass and metals history of a tree while mapping
    the descendant.
    
    Input:  - SAM_choice (int): (0) SAGE (1) Dusty-SAGE
            - a tree yielded by iterate_tree(directory)
            - snap_limit (int) -- last snapshot of the tree
    Output: 3-dimensions array.
            1st array (1xM)*: Stellar mass history of the tree (in Msun/h)
            2nd array (1xM)*: Stellar metallicity history (no unit)

    *M is the number of snapshot, ascending with increasing age of Universe.
    """
    
    recycle_fraction = 0.43
    sorted_idx = np.argsort(tree, order=('GalaxyIndex', 'SnapNum'))
    all_gal_ID = tree['GalaxyIndex']
    snapshot_nr = tree['SnapNum']
    merge_idx = tree['mergeIntoID']
    merge_snapshot = tree['mergeIntoSnapNum']
    merge_type = tree['mergeType']

    delta_bulge = tree['SfrBulge'] * tree['dT'] * 1.e6 * (1.0 - recycle_fraction)
    delta_disk = tree['SfrDisk'] * tree['dT'] * 1.e6 * (1.0 - recycle_fraction)
    delta_mass = delta_bulge + delta_disk

    
    if SAM_choice == 0:
        delta_bulge_metals = tree['SfrBulgeZ']  * tree['SfrBulge'] * tree['dT'] * 1.e6 * (1.0 - recycle_fraction)
        delta_disk_metals = tree['SfrDiskZ'] * tree['SfrDisk'] * tree['dT'] * 1.e6 * (1.0 - recycle_fraction)
        delta_metals = delta_bulge_metals + delta_disk_metals 

    elif SAM_choice == 1:
        delta_bulge_metals = tree['SfrBulgeZ']  * tree['SfrBulge'] * tree['dT'] * 1.e6 * (1.0 - recycle_fraction)
        delta_disk_metals = tree['SfrDiskZ'] * tree['SfrDisk'] * tree['dT'] * 1.e6 * (1.0 - recycle_fraction)
        delta_bulge_dust = tree['SfrBulgeDTG'] * tree['SfrBulge'] * tree['dT'] * 1.e6 * (1.0 - recycle_fraction)
        delta_disk_dust = tree['SfrDiskDTG'] * tree['SfrDisk'] * tree['dT'] * 1.e6 * (1.0 - recycle_fraction)
        delta_metals = delta_bulge_metals + delta_disk_metals + delta_bulge_dust + delta_disk_dust
    else:
        print("Choose a SAM: 0 - for SAGE, 1 - for Dusty-SAGE")
    
    unique_ID = np.unique(all_gal_ID)
    mass = np.zeros((len(unique_ID), max(snapshot_nr)+1))
    metals = np.zeros((len(unique_ID), max(snapshot_nr)+1))
    dust = np.zeros((len(unique_ID), max(snapshot_nr)+1))
    
    #map descendant and build mass and metal history
    for kk, gal_ID in enumerate(unique_ID):
        instant_mass = 0.0
        instant_metals = 0.0
        for ii, ID in enumerate(all_gal_ID[sorted_idx]):
            if(gal_ID == ID):
                instant_mass += delta_mass[sorted_idx[ii]] 
                mass[kk][snapshot_nr[sorted_idx[ii]]] = instant_mass
                assert mass[kk][snapshot_nr[sorted_idx[ii]]] >= mass[kk][snapshot_nr[sorted_idx[ii-1]]]
                
                instant_metals += delta_metals[sorted_idx[ii]]
                metals[kk][snapshot_nr[sorted_idx[ii]]] = instant_metals
                assert metals[kk][snapshot_nr[sorted_idx[ii]]] >= metals[kk][snapshot_nr[sorted_idx[ii-1]]]
                
    
    #make sure the mass and metals are increasing with snapshot_nr   
    for i in range(len(unique_ID)):
        for j in range(max(snapshot_nr)):
            if (mass[i][j+1] < mass[i][j]):
                mass[i][j+1] = mass[i][j]
                
            if (metals[i][j+1] < metals[i][j]):
                metals[i][j+1] = metals[i][j]
            
    #identify merger and add mass
    for snap in range(max(snapshot_nr)+1):
        wsnap = np.where(snapshot_nr == snap)[0]
        wmerge = np.where((merge_idx[wsnap] != -1) & (merge_snapshot[wsnap] < snap_limit) & (merge_type[wsnap] < 3))[0] #only include major (merge_type=1) and minor (merge_type=2) merger
        merger_snap = merge_snapshot[wsnap][wmerge]
        merger_id = merge_idx[wsnap][wmerge]
        
        if len(merger_id) > 0:
            for i, idx in enumerate(merger_id):
                
                wmergesnap = np.where(snapshot_nr == merger_snap[i])[0]
                central_ID = all_gal_ID[wmergesnap][idx]
                central_idx = np.where(unique_ID[:,None] == central_ID)[0]
                satellite_ID = all_gal_ID[wsnap][wmerge][i]
                satellite_idx = np.where(unique_ID[:,None] == satellite_ID)[0]
                
                #added satellite mass to central mass
                mass[central_idx] = mass[central_idx] + mass[satellite_idx]
            
                #eliminate the mass of satellite galaxies
                mass[satellite_idx] = np.zeros(max(snapshot_nr)+1)
                
                #added satellite metals to central
                metals[central_idx] = metals[central_idx] + metals[satellite_idx]
                
                #eliminate the metals of satellite galaxies
                metals[satellite_idx] = np.zeros(max(snapshot_nr)+1)
               
        #null more satellite (from mergetype 3 and 4)
        wmerge = np.where((merge_idx[wsnap] != -1) & (merge_snapshot[wsnap] < snap_limit) & (merge_type[wsnap] > 2))[0]
        merger_snap = merge_snapshot[wsnap][wmerge]
        merger_id = merge_idx[wsnap][wmerge]
        

        if len(merger_id) > 0:
            for i, idx in enumerate(merger_id):
                wmergesnap = np.where(snapshot_nr == merger_snap[i])[0]
                satellite_ID = all_gal_ID[wsnap][wmerge][i]
                satellite_idx = np.where(unique_ID[:,None] == satellite_ID)[0]
                
                #eliminate the mass of satellite galaxies but don't add it to the central
                mass[satellite_idx] = np.zeros(max(snapshot_nr)+1)
                
                #the metals as well
                metals[satellite_idx] = np.zeros(max(snapshot_nr)+1)
                

        #Finally, divide total metals to total mass:
        w = np.where((metals[:,snap] !=0) & (mass[:,snap] != 0))[0]
        metals[w,snap] = metals[w,snap] / mass[w,snap]
                
    return mass, metals
#-----------------------------------------------------------------------------------	
def calculate_dust_density (tree, snap_limit):
    """Calculate mass history from Dusty-SAGE tree.
    In one fly, it will calculate the mass and metals history of a tree while mapping
    the descendant.
    
    Input:  - SAM_choice (int): (0) SAGE (1) Dusty-SAGE
            - a tree yielded by iterate_tree(directory)
            - snap_limit (int) -- last snapshot of the tree
    Output: 4 variable.
            1st variable: Dust mass history (Msun/h)
            2nd variable: Gas phase metal mass history (Msun/h)
            3rd variable: Cold gas mass history (Msun/h)
            4th variable: Disk scale radius history (Mpc/h)
    """
    
    sorted_idx = np.argsort(tree, order=('GalaxyIndex', 'SnapNum'))
    all_gal_ID = tree['GalaxyIndex']
    snapshot_nr = tree['SnapNum']
    merge_idx = tree['mergeIntoID']
    merge_snapshot = tree['mergeIntoSnapNum']
    merge_type = tree['mergeType']
    
    instant_dust = tree['ColdDust'] * 1e10
    instant_metals = tree['MetalsColdGas'] * 1e10
    instant_gas = tree['ColdGas'] * 1e10
    instant_rad = tree['DiskRadius'] #in Mpc
    
    unique_ID = np.unique(all_gal_ID)
    dust = np.zeros((len(unique_ID), max(snapshot_nr)+1))
    gas_metals = np.zeros((len(unique_ID), max(snapshot_nr)+1))
    gas = np.zeros((len(unique_ID), max(snapshot_nr)+1))
    rad = np.zeros((len(unique_ID), max(snapshot_nr)+1))
    
    #map descendant and build mass and metal history
    for kk, gal_ID in enumerate(unique_ID):
        for ii, ID in enumerate(all_gal_ID[sorted_idx]):
            if(gal_ID == ID):
                dust[kk][snapshot_nr[sorted_idx[ii]]] = instant_dust[sorted_idx[ii]]
                gas_metals[kk][snapshot_nr[sorted_idx[ii]]] = instant_metals[sorted_idx[ii]]
                gas[kk][snapshot_nr[sorted_idx[ii]]] = instant_gas[sorted_idx[ii]]
                rad[kk][snapshot_nr[sorted_idx[ii]]] = instant_rad[sorted_idx[ii]]
    
    #identify merger and add mass
    for snap in range(max(snapshot_nr)+1):
        wsnap = np.where(snapshot_nr == snap)[0]
        wmerge = np.where((merge_idx[wsnap] != -1) & (merge_snapshot[wsnap] < snap_limit) & (merge_type[wsnap] < 3))[0] #only include major (merge_type=1) and minor (merge_type=2) merger
        merger_snap = merge_snapshot[wsnap][wmerge]
        merger_id = merge_idx[wsnap][wmerge]
        
        if len(merger_id) > 0:
            for i, idx in enumerate(merger_id):
                
                wmergesnap = np.where(snapshot_nr == merger_snap[i])[0]
                central_ID = all_gal_ID[wmergesnap][idx]
                central_idx = np.where(unique_ID[:,None] == central_ID)[0]
                satellite_ID = all_gal_ID[wsnap][wmerge][i]
                satellite_idx = np.where(unique_ID[:,None] == satellite_ID)[0]
                
                #eliminate the dust of satellite galaxies
                dust[satellite_idx] = np.zeros(max(snapshot_nr)+1)
                gas_metals[satellite_idx] = np.zeros(max(snapshot_nr)+1)
                gas[satellite_idx] = np.zeros(max(snapshot_nr)+1)
                rad[satellite_idx] = np.zeros(max(snapshot_nr)+1)
               
        #null more satellite (from mergetype 3 and 4)
        wmerge = np.where((merge_idx[wsnap] != -1) & (merge_snapshot[wsnap] < snap_limit) & (merge_type[wsnap] > 2))[0]
        merger_snap = merge_snapshot[wsnap][wmerge]
        merger_id = merge_idx[wsnap][wmerge]
        

        if len(merger_id) > 0:
            for i, idx in enumerate(merger_id):
                wmergesnap = np.where(snapshot_nr == merger_snap[i])[0]
                satellite_ID = all_gal_ID[wsnap][wmerge][i]
                satellite_idx = np.where(unique_ID[:,None] == satellite_ID)[0]
                
                #also the dust
                dust[satellite_idx] = np.zeros(max(snapshot_nr)+1)
                gas_metals[satellite_idx] = np.zeros(max(snapshot_nr)+1)
                gas[satellite_idx] = np.zeros(max(snapshot_nr)+1)
                rad[satellite_idx] = np.zeros(max(snapshot_nr)+1)
                
    return dust, gas_metals, gas, rad

#-----------------------------------------------------------------------------------	

def build_mass_and_metallicity_history(SAM_choice, directory, firstfile, lastfile, snap_limit):
    '''
    Build mass and metallicity history from the output directory of dusty-sage
    Input:  - SAM_choice (int): (0) SAGE (1) Dusty-SAGE
            - directory (string) -- path to the directory containing dusty-sage output tree
            - snap_limit (integer) -- number of last snapshot
            
    Output: - Mass (array(list(float))) -- an array containing a number of galaxy, each containing mass (in Msun/h) of each snapshot
            - Metallicity (array(list(float))) -- an array containing a number of galaxy, each containing stellar metallicity of each snapshot

    '''
    Mass = []
    Metals = []

    for tree in iterate_trees(SAM_choice, directory, firstfile, lastfile):
        mass, metal = calculate_mass_and_metals(SAM_choice, tree, snap_limit)
        Mass.extend(mass)
        Metals.extend(metal)

    Mass = np.array(Mass)
    Metals = np.array(Metals)
    
    return(Mass, Metals)
#-----------------------------------------------------------------------------------	

def build_dust_history(SAM_choice, directory, firstfile, lastfile, snap_limit):
    '''
    Build mass and metallicity history from the output directory of dusty-sage
    Input:  - SAM_choice (int): (0) SAGE (1) Dusty-SAGE
            - directory (string) -- path to the directory containing dusty-sage output tree
            - snap_limit (integer) -- index of the last snapshot
            
    Output: - Dust (array(list(float))) -- an array containing a number of galaxy, each containing mass (in Msun/h) of each snapshot
            - Gas metals (array(list(float))) -- an array containing a number of galaxy, each containing gas phase metal mass (in Msun/h) of each snapshot
            - Gas (array(list(float))) -- an array containing a number of galaxy, each containing gas mass (in Msun/h) of each snapshot
            - Rad (array(list(float))) -- an array containing a number of galaxy, each containing Disk Scale Radius (in Mpc/h) of each snapshot
    '''
    Dust = []
    GasMetals = []
    Gas = []
    Rad = []
    
    for tree in iterate_trees(SAM_choice, directory, firstfile, lastfile):
        dust, gas_metal, gas, rad = calculate_dust_density(tree, snap_limit)
        Dust.extend(dust)
        GasMetals.extend(gas_metal)
        Gas.extend(gas)
        Rad.extend(rad)

    Dust = np.array(Dust)
    GasMetals = np.array(GasMetals)
    Gas = np.array(Gas)
    Rad = np.array(Rad)
    
    return(Dust, GasMetals, Gas, Rad)

#-----------------------------------------------------------------------------------	

def open_file(filename):
    
    '''
    Open file, read each line, split each line into each float number.
    Create an list consists of all the number (float) in the file.
    
    Input: - filename (str) -- name of file to be opened
    Output: - M (list(float)) -- a list of float consists of all number in the file named filename
    '''
    
    f = open(filename, "r")
    M = []
    for elem in f.read().split():
        try:
            M.append(float(elem))
        except ValueError:
            pass
    f.close()
    
    return M
#-----------------------------------------------------------------------------------	
def generate_SED(SSP, Age, MassHist, MetalHist, tau_head_BC, tau_head_ISM, eta_BC, eta_ISM, time_BC):
    
    '''
    Generate intrinsice (stellar) SED by assembling SSP from BC03.
    
    Input:  - Choice_of_SSP (int) : 0 - BC03
            - Age : 1-dimension array consists of age of universe in Gyr
            - MassHist: N-dimension array, with N=number of galaxy.
                        Each array consists of stellar mass (in Msun) of each galaxy at corresponding age of Universe.
            - MetalHist: N-dimension array, with N=number of galaxy.
                         Each array consists of stellar metallicity (metals/stellar mass) of each galaxy at corresponding age of Universe.
                         
    Output: - Wavelength: 1-dimension array with 6900 wavelength in micrometer.
            - Luminosity: N-dimension array, with N=number of galaxy.
                          Each array consists of luminosity of galaxy at corresponding wavelength.
    
    '''
    
    #SSP = 0 (Bruzual & Charlot 2003 -- BC03)
    FileNames = ["files/bc2003_hr_m22_chab_ssp.ised_ASCII", "files/bc2003_hr_m32_chab_ssp.ised_ASCII",
            "files/bc2003_hr_m42_chab_ssp.ised_ASCII", "files/bc2003_hr_m52_chab_ssp.ised_ASCII", 
            "files/bc2003_hr_m62_chab_ssp.ised_ASCII", "files/bc2003_hr_m72_chab_ssp.ised_ASCII"]

    AllFiles = []
    for i in range(len(FileNames)):
        AllFiles.append(open_file(FileNames[i]))

    File1 = AllFiles[0]
    lookback = np.array(File1[1:222])
    wavelength = np.array(File1[236: 7136])
    metallicity = [0.0001, 0.0004, 0.004, 0.008, 0.02, 0.05] #metallicity grid in BC03
    time_grid = 221
    wave_grid = 6900
    lum = np.zeros((len(AllFiles), time_grid, wave_grid))

    for j in range(len(metallicity)):
        File = AllFiles[j]
        for i in range(time_grid):
            lum[j][i] = File[7137 + (wave_grid+54)*i: 7137 + (wave_grid+54)*i + 6900]

    #Check if all mass and metal history have the same number of timestep
    if(len(MassHist) > 1):
        for i in range(len(MassHist)):
            if len(MassHist[i]) != len(MassHist[0]):
                print("Not all galaxies have mass history at snapshot=",i)
            if len(MetalHist[i]) != len(MetalHist[0]):
                print("Not all galaxies have metal history at snapshot=", i)

    w = np.where(lookback < 10**7)[0]

    if(len(MassHist) > 1):
        gal_number = len(MassHist)
    else:
        gal_number = 1

    new_mass_hist = np.zeros((time_grid, gal_number))
    new_metal_hist = np.zeros((time_grid, gal_number))

    #   Build new mass and metal history based on the lookback time of BC03
    for i in range(gal_number):
        lookbacktime = age_to_lookback(Age)
        sorted_lbtime = sorted(lookbacktime)

        temp_mass_list = list(MassHist[i])
        temp_metal_list = list(MetalHist[i])
        temp_mass_list.reverse()
        temp_metal_list.reverse()
        new_mass_hist[:,i] = np.interp(lookback, sorted_lbtime, temp_mass_list) 
        new_metal_hist[:,i] = np.interp(lookback, sorted_lbtime, temp_metal_list) 

    #Count half metallicity of metallicity grids from BC03
    half_metal = [0] * (len(metallicity) - 1)
    for i in range(len(half_metal)):
        half_metal[i] = (metallicity[i] + metallicity[i+1]) / 2

    #print('Building SED')
    total_lum = np.zeros((gal_number, wave_grid))
    total_lum_dusty = np.zeros((gal_number, wave_grid))

    tau_ISM = compute_tau(tau_head_ISM, eta_ISM, wavelength)
    tau_BC = tau_ISM + compute_tau(tau_head_BC, eta_BC, wavelength)
    
    
    attenuation_factor_BC = np.e**(-tau_BC)
    attenuation_factor_ISM = np.e**(-tau_ISM)
    
    for i in range(len(lookback) - 1):
        
        if i%22 == 0:
            print(int(i*100/219),'%',end = '')
        else: 
            print('.',end = '')
        
        #print("Timestep", i, "/", len(lookback) - 2)
        delta_mass = new_mass_hist[i] - new_mass_hist[i+1]
        deltamass = np.reshape(delta_mass, (-1, 1))

        w1 = np.where(new_metal_hist[i] < half_metal[0])[0]
        total_lum[w1] += deltamass[w1] * lum[0][i]

        for j in range(len(half_metal)-1):
            w2 = np.where((new_metal_hist[i] > half_metal[j]) & (new_metal_hist[i] <= half_metal[j+1]))[0]
            total_lum[w2] += deltamass[w2] * lum[j+1][i]

        w3 = np.where(new_metal_hist[i] > half_metal[-1])[0]
        total_lum[w3] += deltamass[w3] * lum[-1][i]


        #birthcloud extinction
        if (lookback[i] < time_BC):
            total_lum_dusty = total_lum * attenuation_factor_BC
        else:
            total_lum_dusty = total_lum * attenuation_factor_ISM
        
    return wavelength, total_lum, total_lum_dusty
#-----------------------------------------------------------------------------------	
def compute_tau(tau_head, eta, wavelength):
    
    '''
    Compute optical depth as a function of wavelength using Charlot and Fall (2000) model
    
    Input:  - tau_head (float or array): optical depth at 5500 Angstorm
            - eta (float or array): power law index
            - wavelength (array): in Angstorn           
    Output: - tau (array): the computed optical depth have the same dimension with wavelength. 
    '''

    
    if type(tau_head) != float:
        tau = np.zeros((len(tau_head), len(wavelength)))
        for i in range(len(tau)):
            tau[i] = tau_head[i] * (wavelength/5500)**eta[i]
    else:
        tau = tau_head * (wavelength/5500)**eta

    return(tau)
#-----------------------------------------------------------------------------------	
def compute_area(rad):
    '''
    Compute area of a disk (2 * area of a circle)
    
    Input:  - rad (float or array): radius
    Output: - area (float or array): area of a disk for given radius
    '''

    return (2 * np.pi * rad**2)
#-----------------------------------------------------------------------------------	
def compute_tauBC_Trayford(ColdDust, ColdGas, rad):
    '''
    Compute optical depth at 5500 A for birth clouds using relation in Trayford+ 19
    (relation between dust surface density and optical depth)
    
    Input:  - ColdDust (float or array): total cold dust mass in the ISM in Msun
            - ColdGas (float or array): total cold gas mass in the ISM in Msun
            - rad (float or array): disk radius in Mpc
    Output: - Sigma_BC (float or array): dust surface density in the birth clouds (logscale Msun/kpc2)
            - tau_BC (float or array): optical depth at 5500 A for birth clouds
    '''

    
    fdust = ColdDust / ColdGas
    ScaleRad = rad * 1e6 #convert to pc
    halfrad = 1.68 * ScaleRad 
    threerad = 0.4 * ScaleRad
    
    fdust_MW = 0.33
    Zsun = 0.0189
    Sigma_MW = 85 #Msun/pc2 
    
    #area = compute_area(halfrad)
    area = compute_area(threerad)
    Sigma_gas = ColdGas / area
    w = np.where(Sigma_gas < Sigma_MW)[0]
    if len(w) > 1:
        Sigma_gas[w] = Sigma_MW
    elif len(w) == 1:
        Sigma_gas = Sigma_MW
    tau_BC = (fdust * Sigma_gas) / (fdust_MW * Zsun * Sigma_MW)
    Sigma_BC = np.log10(fdust * Sigma_gas * 1e6) #kpc
    
    return Sigma_BC, tau_BC

#-----------------------------------------------------------------------------------	
def compute_tauISM_Trayford(ColdDust, ColdGas, rad):
    '''
    Compute ISM optical depth from relation in Trayford+ 19
    
    Input:  - ColdDust (float or array): total cold dust mass in the ISM (Msun)
            - ColdGas (float or array): total cold gas mass in the ISM (Msun)
            - rad (float or array): disk radius in Mpc
    Output: - Sigma_dust (float or array): dust surface density in the diffuse ISM (in logscale - Msun/kpc2)
            - tau_ISM (float or array): optical depth at 5500 A for diffuse ISM
    '''

    ScaleRad = rad * 1e3 #convert to kpc
    halfrad = 1.68 * ScaleRad
    threerad = 0.4 * ScaleRad
    Sigma_MW = 85 * 1e6 #Msun/kpc2
    
    Sigma_dust_Trayford = [4.088, 4.351, 4.579, 4.823, 5.057, 5.292, 5.528, 5.765, 6.001, 6.234, 6.470, 6.704, 6.941, 7.177, 7.416]
    tau_head_Trayford = [0.031, 0.059, 0.078, 0.129, 0.203, 0.308, 0.467, 0.647, 0.838, 1.065, 1.235, 1.475, 1.571, 1.645, 1.806]

    #area = compute_area(halfrad)
    area = compute_area(threerad)
        
    Sigma_dust = np.log10(ColdDust / area)
    #Sigma_dust = np.log10(fdust * Sigma_gas)
    tau_ISM = np.interp(Sigma_dust, Sigma_dust_Trayford, tau_head_Trayford)
    
    return Sigma_dust, tau_ISM

#-----------------------------------------------------------------------------------	

def compute_etaISM_Trayford(ColdDust, ColdGas, rad):
    '''
    Compute powerlaw index for the diffuse ISM component from relation in Trayford+ 19
    
    Input:  - ColdDust (float or array): total cold dust mass in the ISM in Msun
            - ColdGas (float or array): total cold gas mass in the ISM in Msun
            - rad (float or array): disk radius in Mpc
    Output: - Sigma_dust (float or array): dust surface density in the diffuse ISM (logscale Msun/kpc2)
            - eta_ISM (float or array): powerlaw index for diffuse ISM
    '''


    ScaleRad = rad * 1e3 #convert to kpc
    halfrad = 1.68 * ScaleRad
    threerad = 0.4 * ScaleRad
    
    Sigma_dust_Trayford = [4.116, 4.354, 4.588, 4.822, 5.061, 5.293, 5.533, 5.763, 6.003, 6.236, 6.469, 6.706, 6.943, 7.181, 7.411]
    eta_ISM_Trayford = [-1.379, -1.357, -1.334, -1.243, -1.170, -1.062, -0.922, -0.778, -0.668, -0.570, -0.506, -0.453, -0.381, -0.318, -0.307]
    #area = compute_area(halfrad)
    area = compute_area(threerad)
    
    Sigma_dust = np.log10(ColdDust / area)
    #Sigma_dust = np.log10(fdust * Sigma_gas)
    eta_ISM = np.interp(Sigma_dust, Sigma_dust_Trayford, eta_ISM_Trayford)
    
    return Sigma_dust, eta_ISM

#-----------------------------------------------------------------------------------	

def compute_tauISM_Somerville (Dust, Rad):
    '''
    Compute ISM optical depth from prescription in Somerville+ 2012
    
    Input:  - ColdDust (float or array): total cold dust mass in the ISM in Msun
            - rad (float or array): disk radius in Mpc
    Output: - Sigma_dust (float or array): dust surface density in the diffuse ISM (Msun/pc2)
            - tau_ISM (float or array): optical depth at 5500 A for diffuse ISM
    '''

    
    Chi_gas = 0.42
    rad_gas = Chi_gas * Rad * 1e6
    
    tau_dust_0 = 0.3
    Sigma = Dust / (rad_gas**2) 
    tau_ISM = tau_dust_0 * Sigma
    
    return(Sigma, tau_ISM)

#-----------------------------------------------------------------------------------	

def compute_tauBC_Somerville (Dust, Rad):
    
    
    '''
    Compute birth clouds' optical depth from prescription in Somerville+ 2012
    
    Input:  - ColdDust (float or array): total cold dust mass in the ISM in Msun
            - rad (float or array): disk radius in Mpc
    Output: - Sigma_dust (float or array): dust surface density in the Birth Clouds (Msun/pc2)
            - tau_ISM (float or array): optical depth at 5500 A for Birth Clouds
    '''

    
    Chi_gas = 0.42
    rad_gas = Chi_gas * Rad * 1e6
    
    tau_dust_0 = 0.3
    Sigma = Dust / (rad_gas**2)
    tau_ISM = tau_dust_0 * Sigma
    
    mu_BC = 6
    tau_BC = mu_BC * tau_ISM
    
    return (Sigma, tau_BC)

#-----------------------------------------------------------------------------------	

def compute_attenuation_parameters (prescription_choice, DustMass, GasMass, Radius):
    '''
    Compute attenuation parameters based on the Charlot & Fall (2000) model.
    We adopted two prescriptions to compute the parameters:
    Input:  - prescription_choice (int): (0) Lagos+ 19 (1) Somerville+ 12
            - DustMass (float or array): Dust mass in Msun
            - GasMass (float or array): Gas mass in Msun
            
    Output: - Mass (array(list(float))) -- an array containing a number of galaxy, each containing mass (in Msun/h) of each snapshot
            - Metallicity (array(list(float))) -- an array containing a number of galaxy, each containing stellar metallicity of each snapshot

    '''

    
    eta_BC = [-0.7] * len(DustMass)
    
    if prescription_choice == 0:
        Sigma_BC, tau_BC = compute_tauBC_Trayford(DustMass, GasMass, Radius)
        Sigma_tau_ISM, tau_ISM = compute_tauISM_Trayford(DustMass, GasMass, Radius)
        Sigma_eta_ISM, eta_ISM = compute_etaISM_Trayford(DustMass, GasMass, Radius)

    elif prescription_choice == 1:
        eta_ISM = [-1.3] * len(DustMass)
        Sigma_ISM, tau_ISM = compute_tauISM_Somerville (DustMass, Radius)
        Sigma_BC, tau_BC = compute_tauBC_Somerville (DustMass, Radius)
    
    else:
        print("Choose 0 for attenuation prescriptions from Lagos+19 and 1 for Somerville+12")
        
    return tau_BC, eta_BC, tau_ISM, eta_ISM

#-----------------------------------------------------------------------------------	
def determine_idx_Rieke(LIR):
    
    '''
    Determine the index of the Dale+ 14 IR template to be used based on the total IR luminosity
    from the prescription in Rieke+ 09
    
    Input:  - LIR (float or array): total IR luminosity
    Output: - idx (float or array): index of the spectra in the Dale+ 14 IR template
    '''

    alpha_SF, log_fnu_SF = np.loadtxt('files/alpha.dat', unpack=True)
    
    if (LIR > 10**11.6):  
        LIR = 10**11.6

    alpha = 10.096 - 0.741 * np.log10(LIR)

    delta_alpha = abs(alpha_SF - alpha)
    idx = np.where(delta_alpha==min(delta_alpha))[0]    
    return idx

#-----------------------------------------------------------------------------------	

def determine_idx_Marcillac(LIR):
    '''
    Determine the index of the Dale+ 14 IR template to be used based on the total IR luminosity
    from the prescription in Marcillac+ 06
    
    Input:  - LIR (float or array): total IR luminosity
    Output: - idx (float or array): index of the spectra in the Dale+ 14 IR template
    '''

    alpha_SF, log_fnu_SF = np.loadtxt('files/alpha.dat', unpack=True)
    log_fnu = 0.128 * np.log10(LIR) - 1.611
    delta_fnu = abs(log_fnu_SF - log_fnu)
    idx = np.where(delta_fnu==min(delta_fnu))[0]    
    return idx

#-----------------------------------------------------------------------------------	

def add_IR_Dale (wavelength, spectra, spectra_dusty):

    '''
    Add the NIR-FIR spectra from Dale+ 14 IR template to the UV-NIR spectra from BC03/
    
    Input:  - wavelength (N-dimensional array)
            - spectra (N-dimensional array) - intrinsic stellar spectra corresponding to each wavelength
            - spectra_dusty (N-dimensional array) - attenuated spectra corresponding to each wavelength
    Output: - wavelength (M-dimensional array) - M > N
            - spectra (M-dimensional array) - attenuated spectra with IR addition
    '''
    
    Dale_template = np.loadtxt('files/spectra.0.00AGN.dat', unpack=True)
    lambda_IR = Dale_template[0] * 1e4 #convert from micron to Angstrom

    Ldust = (spectra - spectra_dusty)
    w = np.where(wavelength < 912)[0]
    idx_912 = w[-1]

    all_wave = np.unique(np.concatenate((wavelength, lambda_IR)))
    all_wave.sort(kind='mergesort')
    UVIR = np.zeros((len(Ldust), len(all_wave)))

    for i in range(len(Ldust)):
        LIR_mentari = np.trapz(Ldust[i][idx_912:-1], wavelength[idx_912:-1])
        #---------------------------------------------------
        '''
        #Compute alpha based on Rieke+ 2009
        if (LIR_mentari > 10**11.6):  
            LIR_mentari = 10**11.6

        alpha = 10.096 - 0.741 * np.log10(LIR_mentari)

        delta_alpha = abs(alpha_SF - alpha)
        idx = np.where(delta_alpha==min(delta_alpha))[0]
        
        #----------------------------------------------------
        #Compute log_fnu based on Marcillac+ 2006
        log_fnu = 0.128 * np.log10(LIR_mentari) - 1.611
        delta_fnu = abs(log_fnu_SF - log_fnu)
        idx = np.where(delta_fnu==min(delta_fnu))[0]
        #----------------------------------------------------
        '''
        idx = determine_idx_Marcillac(LIR_mentari)
        spectra_IR = 10 ** Dale_template[idx[0]+1] 

        LIR_dale = np.trapz(spectra_IR, lambda_IR)
        scaling = LIR_mentari / LIR_dale
        spectra_IR_dale = spectra_IR * scaling 

        new_spectra = np.interp(all_wave, wavelength, spectra_dusty[i])
        new_IR = np.interp(all_wave, lambda_IR, spectra_IR_dale)
        all_spec = new_spectra + new_IR
        UVIR[i] = all_spec

    return (all_wave, UVIR)

#-----------------------------------------------------------------------------------	

"""
Contains some simple functions for working with the Safarzadeh et al. (2015)
FIR SED templates.

The find_template function returns the SED template that
is the best match for the log (L_IR/L_sun) and log (M_dust/M_sun) values
specified by the user.

The load_templates function reads the ASCII file that
contains the template data and returns arrays that contain the L_IR and M_dust
values for each template, a 2-D array that contains the normalized SEDs
(lambda*L_lambda/L_IR), and the wavelength array.

"""
def find_template_SUNRISE(
    lir, #"""logarithm of the IR luminosity in solar units"""
    mdust, #"""logarithm of the dust mass in solar units"""
    tol=0.5, #"""maximum allowed discrepancy between template and requested values, delta(L_IR, M_dust); see function docstring"""
    rtol=0.2, #"""maximum allowed difference between the template and requested L_IR/M_dust ratio"""
    sed_file="safarzadeh_et_al_2015.txt", #"""full path to file that contains the templates"""
    verbose=True #"""prints messages if set to True"""
    ):
    """Return template with (L_IR, M_dust) values closest to those requested.

    Similarity is defined by the Euclidean distance between the requested and template 
    (L_IR, M_dust) values,

        delta(L_IR, M_dust) = [(L_IR,requested - L_IR, template)^2
            + (M_dust,requested - M_dust,template)^2]^(0.5)

    If no template with delta(L_IR, M_dust) < tol is found, then search for a template with an
    L_IR/M_dust ratio that differs from that requested by at most rtol.

    Parameters:
        - lir : log (L_IR/L_sun) value for which a template is desired
        - mdust : log (M_dust/M_sun) value for which a template is desired
        - tol : maximum allowed delta(L_IR, M_dust) value
        - rtol : maximum allowed difference between requested and template L_IR/M_dust ratio;
        only used if delta(L_IR, M_dust) < tol cannot be satisfied
        - sed_file : full path to file that contains the template
        - verbose : print helpful messages if set to True

    Returns:
        - lambda_array: array of wavelengths (in micron)
        - sed: template SED that is most appropriate for the (L_IR, M_dust) requested; the
            template is normalized such that its total IR luminosity equals  the *requested*
            L_IR value (i.e. the normalized SED from the template file is multipled by L_IR

    Example:
        Find an SED for a galaxy with L_IR = 10^11 L_sun and M_dust = 10^8 M_sun:

            lam, sed = find_template(11.,8.)
"""


    lir_array, mdust_array, sed_array, lambda_array = load_templates_SUNRISE()

    dist = ((lir_array - lir)**2 + (mdust_array - mdust)**2)**(0.5)
    
    # first try to find a template with delta(L_IR,M_dust) < tol (see function docstring)
    if dist.min() < tol :
        '''
        if verbose :
            print ("Template sufficiently close in the (L_IR, M_dust) plane found")
            print ("Requested: log L_IR = ", lir, ", log M_dust = ", mdust)
            print ("Template: log L_IR = ", lir_array[dist.argmin()], ", log M_dust = ", \
                mdust_array[dist.argmin()])
        '''
        return lambda_array, sed_array[dist.argmin(),:]*10.**lir
    else : # not found, so search for one with delta(L_IR/M_dust) < rtol (see function docstring)
        '''
        if verbose :
            print ("Template with delta(L_IR, M_dust) < ",tol," not found. Returning template with nearest L_IR/M_dust value instead")
        '''
        ratio=lir-mdust
        ratio_array=lir_array-mdust_array
        dist_ratio = np.abs(ratio - ratio_array)
        if dist_ratio.min() < rtol :
            '''
            if verbose :
                print ("Found template with sufficiently similar L_IR/M_dust ratio")
                print ("Requested: log L_IR = ", lir, ", log M_dust = ", mdust, \
                    ", log L_IR/M_dust = ",ratio)
                print ("Template: log L_IR = ", lir_array[dist_ratio.argmin()], ", log M_dust = ", \
                    mdust_array[dist_ratio.argmin()],", log L_IR/M_dust = ", \
                    ratio_array[dist_ratio.argmin()])
            '''
            sed = sed_array[dist_ratio.argmin(),:]*10.**lir
            return lambda_array, sed
        else : # didn't find a suitable template
            print ("ERROR: acceptable template not found")
            return

#-----------------------------------------------------------------------------------	

def load_templates_SUNRISE(
    sed_file = "safarzadeh_et_al_2015.txt" #"""full path to file that contains the templates"""
    ):
    """Read the ASCII file that contains the template data.

    Parameters:
        - sed_file: the full path to the template file.

    Returns:
        - lir_array : log (L_IR/L_sun) for each template

        - mdust_array : log (M_dust/M_sun) for each template

        - sed_array: the actual SEDs (lambda*L_lambda) normalized by dividing by L_IR.
        The first dimension is the template number, and the second is wavelength.
    
        - lambda_array: the wavelengths (in microns) at which the SEDs are provided

    Example:
        import template_functions as tf
        lir_array, mdust_array, sed_array, lambda_array = tf.load_templates()
    """

    lir_array = np.loadtxt(sed_file,skiprows=21,usecols=[0])
    mdust_array = np.loadtxt(sed_file,skiprows=21,usecols=[1])
    sed_array = np.loadtxt(sed_file,skiprows=21,usecols=(np.arange(19)+3))
    lambda_array = np.loadtxt(sed_file,skiprows=20,usecols=(np.arange(19)+3))[0,:]*1.e6

    return lir_array, mdust_array, sed_array, lambda_array

#-----------------------------------------------------------------------------------	

def compute_IR_SUNRISE (Dust, wavelength, spectra, spectra_dusty):
    
    '''
    Compute the FIR spectra from SUNRISE model Safarzadeh+ 15    
    Input:  - Dust (array): dust mass in Msun/h
            - wavelength (array): in UV-NIR
            - spectra (array): intrinsic spectra corresponding to each wavelength
            - spectra_dusty (array): attenuated spectra corresponding to each wavelength
    Output: - wave_IR (array): wavelength in FIR
            - lum_IR (array): IR spectra at each corresponding FIR wavelength
    '''
    w = np.where(Dust > 5e5)[0]
    if len(w) > 0:
        new_dust = Dust[w]
        Ldust = (spectra[w] - spectra_dusty[w])
        w = np.where(wavelength < 912)[0]
        idx_912 = w[-1]

        LIR_mentari = np.trapz(Ldust[0][idx_912:-1], wavelength[idx_912:-1])
        lam, sed  = find_template_SUNRISE(np.log10(LIR_mentari), np.log10(new_dust[0]))
        wave_IR = lam * 1e4
        lum_IR = np.zeros((len(Ldust), len(wave_IR)))

        for i in range(len(Ldust)):
            LIR_mentari = np.trapz(Ldust[i][idx_912:-1], wavelength[idx_912:-1])
            lam, sed = find_template_SUNRISE(np.log10(LIR_mentari), np.log10(new_dust[i]))
            lum_IR[i] = sed / (lam * 1e4)   

    else:
	
        wave_IR = 0
        lum_IR = 0
    return (wave_IR, lum_IR)

#-----------------------------------------------------------------------------------	

def combine_Dale_SUNRISE(wavelength_sunrise, spectra_sunrise, wavelength_dale, spectra_dale):
    '''
    Combine UV-MIR spectra from BC03 and Dale+ 14 with the FIR spectra from SUNRISE (Safarzadeh+15)    
    Input:  - wavelength_sunrise (array): in FIR
            - spectra_sunrise (array): FIR spectra from Safarzadeh+14
            - wavelength_dale (array): in UV-NIR
            - spectra_dale (array): combination of BC03 and Dale+14 spectra in UV-NIR
    Output: - all_wave (array): wavelength in UV-IR
            - new_spec (array): combination of BC03 (UV-NIR), Dale+14 (MIR), and Safarzadeh+ 15 (FIR)
    '''

    w = np.where(wavelength_dale < wavelength_sunrise[0])
    all_wave = np.unique(np.concatenate((wavelength_dale[w], wavelength_sunrise)))
    all_wave.sort(kind='mergesort')
    
    new_spec = np.zeros((len(spectra_sunrise), len(all_wave)))
    for i in range(len(spectra_sunrise)):
        w = np.where((all_wave < wavelength_sunrise[0]) | (all_wave == wavelength_sunrise[0]))[0]
        UVIR = np.interp(all_wave[w], wavelength_dale, spectra_dale[i])
        new_spec[i][w] = UVIR
        w = np.where(all_wave > wavelength_sunrise[0])[0]
        new_spec[i,w] = np.interp(all_wave[w], wavelength_sunrise, spectra_sunrise[i])

    return (all_wave, new_spec)

#-----------------------------------------------------------------------------------	

def age_to_lookback(age): #age in Gyr, lookback in yr
    
    '''
    Convert age of Universe to lookback time
    Input: age (list(float)) -- age of universe, in Gyr
    Output: lookback time (list(float)) -- corresponding lookback time, in yr
    '''
    lookback = (np.array([13.6098]*len(age)) - age) * 1.e9
    return lookback

#-----------------------------------------------------------------------------------	

def read_filters():
    
    '''
    Reading filters wavelength and response listed in 'files/allfilters.dat'
    '''
    
    F = type('', (), {})
    F.wavelength, F.response = np.loadtxt('files/allfilters.dat', unpack=True)
    F.Johnson_V_wave = F.wavelength[0:24]
    F.Johnson_V = F.response[0:24]
    F.Johnson_U_wave = F.wavelength[24:49]
    F.Johnson_U = F.response[24:49]
    F.Johnson_B_wave = F.wavelength[49:70]
    F.Johnson_B = F.response[49:70]
    F.Buser_B2_wave = F.wavelength[70:110]
    F.Buser_B2 = F.response[70:110]
    F.Cousins_R_wave = F.wavelength[110:175]
    F.Cousins_R = F.response[110:175]
    F.Cousins_I_wave = F.wavelength[175:214]
    F.Cousins_I = F.response[175:214]
    F.Deep_B_wave = F.wavelength[214:584]
    F.Deep_B = F.response[214:584]
    F.Deep_R_wave = F.wavelength[584:750]
    F.Deep_R = F.response[584:750]
    F.Deep_I_wave = F.wavelength[750:1106]
    F.Deep_I = F.response[750:1106]
    F.TwoMass_J_wave = F.wavelength[1106:1214]
    F.TwoMass_J = F.response[1106:1214]
    F.TwoMass_H_wave = F.wavelength[1214:1272]
    F.TwoMass_H = F.response[1214:1272]
    F.TwoMass_Ks_wave = F.wavelength[1272:1347]
    F.TwoMass_Ks = F.response[1272:1347]
    F.Sdss_u_wave = F.wavelength[1347:1394]
    F.Sdss_u = F.response[1347:1394]
    F.Sdss_g_wave = F.wavelength[1394:1483]
    F.Sdss_g = F.response[1394:1483]
    F.Sdss_r_wave = F.wavelength[1483:1558]
    F.Sdss_r = F.response[1483:1558]
    F.Sdss_i_wave = F.wavelength[1558:1647]
    F.Sdss_i = F.response[1558:1647]
    F.Sdss_z_wave = F.wavelength[1647:1788]
    F.Sdss_z = F.response[1647:1788]
    F.WFPC2_F255W_wave = F.wavelength[1788:11788]
    F.WFPC2_F255W = F.response[1788:11788]
    F.WPC2_F300W_wave = F.wavelength[11788:21788]
    F.WFPC2_F300W = F.response[11788:21788]
    F.WFPC2_F336W_wave = F.wavelength[21788:31788]
    F.WFPC2_F336W = F.response[21788:31788]
    F.WFPC2_F439W_wave = F.wavelength[31788:41788]
    F.WFPC2_F439W = F.response[31788:41788]
    F.WFPC2_F450W_wave = F.wavelength[41788:51788]
    F.WFPC2_F450W = F.response[41788:51788]
    F.WFPC2_F555W_wave = F.wavelength[51788:61788]
    F.WFPC2_F555W = F.response[51788:61788]	
    F.WFPC2_F606W_wave = F.wavelength[61788:71788]
    F.WFPC2_F606W = F.response[61788:71788]
    F.WFPC2_F814W_wave = F.wavelength[71788:81788]
    F.WFPC2_F814W = F.response[71788:81788]
    F.WFPC2_F850W_wave = F.wavelength[81788:91788]
    F.WFPC2_F850W = F.response[81788:91788]
    F.WFCACS_F435W_wave = F.wavelength[91788:101788]
    F.WFCACS_F435W = F.response[91788:101788]
    F.WFCACS_F475W_wave = F.wavelength[101788:111788]
    F.WFCACS_F475W = F.response[101788:111788]
    F.WFCACS_F555W_wave = F.wavelength[111788:121788]
    F.WFCACS_F555W = F.response[111788:121788]
    F.WFCACS_F606W_wave = F.wavelength[121788:131788]
    F.WFCACS_F606W = F.response[121788:131788]
    F.WFCACS_F625W_wave = F.wavelength[131788:141788]
    F.WFCACS_F625W = F.response[131788:141788]
    F.WFCACS_F775W_wave = F.wavelength[141788:151788]
    F.WFCACS_F775W = F.response[141788:151788]
    F.WFCACS_F814W_wave = F.wavelength[151788:161788]
    F.WFCACS_F814W = F.response[151788:161788]
    F.WFCACS_F850W_wave = F.wavelength[161788:171788]
    F.WFCACS_F850W = F.response[161788:171788]
    F.WFC3UVIS_F218W_wave = F.wavelength[171788:180742]
    F.WFC3UVIS_F218W = F.response[171788:180742]
    F.WFC3UVIS_F225W_wave = F.wavelength[180742:189757]
    F.WFC3UVIS_F225W = F.response[180742:189757]
    F.WFC3UVIS_F275W_wave = F.wavelength[189757:198762]
    F.WFC3UVIS_F275W = F.response[189757:198762]
    F.WFC3UVIS_F336W_wave = F.wavelength[198762:207777]
    F.WFC3UVIS_F336W = F.response[198762:207777]
    F.WFC3UVIS_F390W_wave = F.wavelength[207777:216792]
    F.WFC3UVIS_F390W = F.response[207777:216792]	
    F.WFC3UVIS_F438W_wave = F.wavelength[216792:225807]
    F.WFC3UVIS_F438W = F.response[216792:225807]
    F.WFC3UVIS_F475W_wave = F.wavelength[225807:234822]
    F.WFC3UVIS_F475W = F.response[225807:234822]
    F.WFC3UVIS_F555W_wave = F.wavelength[234822:243837]
    F.WFC3UVIS_F555W = F.response[234822:243837]	
    F.WFC3UVIS_F606W_wave = F.wavelength[243837:252792]
    F.WFC3UVIS_F606W = F.response[243837:252792]
    F.WFC3UVIS_F775W_wave = F.wavelength[252792:261807]
    F.WFC3UVIS_F775W = F.response[252792:261807]
    F.WFC3UVIS_F814W_wave = F.wavelength[261807:270822]
    F.WFC3UVIS_F814W = F.response[261807:270822]	
    F.WFC3UVIS_F850W_wave = F.wavelength[270822:279837]
    F.WFC3UVIS_F850W = F.response[270822:279837]
    F.WFC3IR_F098M_wave = F.wavelength[279837:284338]
    F.WFC3IR_F098M = F.response[279837:284338]
    F.WFC3IR_F105W_wave = F.wavelength[284338:293339]
    F.WFC3IR_F105W = F.response[284338:293339]
    F.WFC3IR_F110W_wave = F.wavelength[284338:302340]
    F.WFC3IR_F110W = F.response[284338:302340]
    F.WFC3IR_F125W_wave = F.wavelength[302340:311340]
    F.WFC3IR_F125W = F.response[302340:311340]
    F.WFC3IR_F140W_wave = F.wavelength[311340:320341]
    F.WFC3IR_F140W = F.response[311340:320341]
    F.WFC3IR_F160W_wave = F.wavelength[320341:329342]
    F.WFC3IR_F160W = F.response[320341:329342]
    F.IRAC_1_wave = F.wavelength[329342:329847]
    F.IRAC_1 = F.response[329342:329847]
    F.IRAC_2_wave = F.wavelength[329847:330274]
    F.IRAC_2 = F.response[329847:330274]
    F.IRAC_3_wave = F.wavelength[330274:330644]
    F.IRAC_3 = F.response[330274:330644]
    F.IRAC_4_wave = F.wavelength[330644:331065]
    F.IRAC_4 = F.response[330644:331065]
    F.ISAAC_Ks_wave = F.wavelength[331065:331265]
    F.ISAAC_Ks = F.response[331065:331265]
    F.FORS_V_wave = F.wavelength[331265:331765]
    F.FORS_V = F.response[331265:331765]
    F.FORS_R_wave = F.wavelength[331765:332265]
    F.FORS_R = F.response[331765:332265]	
    F.NIC_F110W_wave = F.wavelength[332265:334264]
    F.NIC_F110W = F.response[332265:334264]	
    F.NIC_F160W_wave = F.wavelength[334264:335868]
    F.NIC_F160W = F.response[334264:335868]	
    F.GALEX_FUV_wave = F.wavelength[335868:336369]
    F.GALEX_FUV = F.response[335868:336369]	
    F.GALEX_NUV_wave = F.wavelength[335868:337710]
    F.GALEX_NUV = F.response[335868:337710]
    F.DES_g_wave = F.wavelength[337710:337900]
    F.DES_g = F.response[337710:337900]	
    F.DES_r_wave = F.wavelength[337900:338100]
    F.DES_r = F.response[337900:338100]	
    F.DES_i_wave = F.wavelength[338100:338290]
    F.DES_i = F.response[338100:338290]
    F.DES_z_wave = F.wavelength[338290:338480]
    F.DES_z = F.response[338290:338480]
    F.DES_Y_wave = F.wavelength[338480:338570]
    F.DES_Y = F.response[338480:338570]
    F.WFCAM_Z_wave = F.wavelength[338570:338723]
    F.WFCAM_Z = F.response[338570:338723]
    F.WFCAM_Y_wave = F.wavelength[338723:338890]
    F.WFCAM_Y = F.response[338723:338890]
    F.WFCAM_J_wave = F.wavelength[338890:339139]
    F.WFCAM_J = F.response[338890:339139]
    F.WFCAM_H_wave = F.wavelength[339139:339642]
    F.WFCAM_H = F.response[339139:339642]
    F.WFCAM_K_wave = F.wavelength[339642:340216]
    F.WFCAM_K = F.response[339642:340216]
    F.Steidel_Un_wave = F.wavelength[340216:340259]
    F.Steidel_Un = F.response[340216:340259]
    F.Steidel_G_wave = F.wavelength[340259:340430]
    F.Steidel_G = F.response[340259:340430]
    F.Steidel_Rs_wave = F.wavelength[340430:341239]
    F.Steidel_Rs = F.response[340430:341239]
    F.Steidel_I_wave = F.wavelength[341239:341636]
    F.Steidel_I = F.response[341239:341636]
    F.MegaCam_u_wave = F.wavelength[341636:341768]
    F.MegaCam_u = F.response[341636:341768]
    F.MegaCam_g_wave = F.wavelength[341768:342009]
    F.MegaCam_g = F.response[341768:342009]
    F.MegaCam_r_wave = F.wavelength[342009:342239]
    F.MegaCam_r = F.response[342009:342239]
    F.MegaCam_i_wave = F.wavelength[342239:342378]
    F.MegaCam_i = F.response[342239:342378]
    F.MegaCam_z_wave = F.wavelength[342378:342530]
    F.MegaCam_z = F.response[342378:342530]
    F.WISE_W1_wave = F.wavelength[342530:342717]
    F.WISE_W1 = F.response[342530:342717]
    F.WISE_W2_wave = F.wavelength[342717:342967]
    F.WISE_W2 = F.response[342717:342967]
    F.WISE_W3_wave = F.wavelength[342967:344467]
    F.WISE_W3 = F.response[342967:344467]
    F.WISE_W4_wave = F.wavelength[344467:345679]
    F.WISE_W4 = F.response[344467:345679]
    F.UVOT_w2_wave = F.wavelength[345679:346320]
    F.UVOT_w2 = F.response[345679:346320]
    F.UVOT_m2_wave = F.wavelength[346320:346636]
    F.UVOT_m2 = F.response[346320:346636]
    F.UVOT_w1_wave = F.wavelength[346636:347177]
    F.UVOT_w1 = F.response[346636:347177]
    F.MIPS_24um_wave = F.wavelength[347177:347305]
    F.MIPS_24um = F.response[347177:347305]
    F.MIPS_70um_wave = F.wavelength[347305:347416]
    F.MIPS_70um = F.response[347305:347416]
    F.MIPS_160um_wave = F.wavelength[347416:347815]
    F.MIPS_160um = F.response[347416:347815]
    F.SCUBA_450WB_wave = F.wavelength[347815:348511]
    F.SCUBA_450WB = F.response[347815:348511]
    F.SCUBA_850WB_wave = F.wavelength[348511:348994]
    F.SCUBA_850WB = F.response[348511:348994]
    F.PACS_70um_wave = F.wavelength[348994:349208]
    F.PACS_70um = F.response[348994:349208]
    F.PACS_100um_wave = F.wavelength[349208:349447]
    F.PACS_100um = F.response[349208:349447]
    F.PACS_160um_wave = F.wavelength[349447:349680]
    F.PACS_160um = F.response[349447:349680]
    F.SPIRE_250um_wave = F.wavelength[349680:349810]
    F.SPIRE_250um = F.response[349680:349810]
    F.SPIRE_350um_wave = F.wavelength[349810:349901]
    F.SPIRE_350um = F.response[349810:349901]
    F.SPIRE_500um_wave = F.wavelength[349901:349999]
    F.SPIRE_500um = F.response[349901:349999]
    F.IRAS_12um_wave = F.wavelength[349999:350017]
    F.IRAS_12um = F.response[349999:350017]
    F.IRAS_25um_wave = F.wavelength[350017:350049]
    F.IRAS_25um = F.response[350017:350049]
    F.IRAS_60um_wave = F.wavelength[350049:350070]
    F.IRAS_60um = F.response[350049:350070]
    F.IRAS_100um_wave = F.wavelength[350070:350086]
    F.IRAS_100um = F.response[350070:350086]
    F.Bessel_L_wave = F.wavelength[350086:350107]
    F.Bessel_L = F.response[350086:350107]
    F.Bessel_Lprime_wave = F.wavelength[350107:350127]
    F.Bessel_Lprime = F.response[350107:350127]
    F.Bessel_M_wave = F.wavelength[350127:350144]
    F.Bessel_M = F.response[350127:350144]
    F.Stromgren_u_wave = F.wavelength[350144:350173]
    F.Stromgren_u = F.response[350144:350173]
    F.Stromgren_v_wave = F.wavelength[350173:350202]
    F.Stromgren_v = F.response[350173:350202]
    F.Stromgren_b_wave = F.wavelength[350202:350231]
    F.Stromgren_b = F.response[350202:350231]
    F.Stromgren_y_wave = F.wavelength[350231:350260]
    F.Stromgren_y = F.response[350231:350260]
    F.Idealized_1500A_wave = F.wavelength[350260:350301]
    F.Idealized_1500A = F.response[350260:350301]
    F.Idealized_2300A_wave = F.wavelength[350301:350362]
    F.Idealized_2300A = F.response[350301:350362]
    F.Idealized_2800A_wave = F.wavelength[350362:350437]
    F.Idealized_2800A = F.response[350362:350437]
    F.JWST_F070W_wave = F.wavelength[350437:350837]
    F.JWST_F070W = F.response[350437:350837]
    F.JWST_F090W_wave = F.wavelength[350837:351139]
    F.JWST_F090W = F.response[350837:351139]
    F.JWST_F115W_wave = F.wavelength[351139:351555]
    F.JWST_F115W = F.response[351139:351555]
    F.JWST_F150W_wave = F.wavelength[351555:352221]
    F.JWST_F150W = F.response[35155:352221]
    F.JWST_F200W_wave = F.wavelength[352221:353128]
    F.JWST_F200W = F.response[352221:353128]
    F.JWST_F277W_wave = F.wavelength[353128:354553]
    F.JWST_F277W = F.response[353128:354553]
    F.JWST_F356W_wave = F.wavelength[354553:355899]
    F.JWST_F356W = F.response[354553:355899]
    F.JWST_F444W_wave = F.wavelength[355899:357351]
    F.JWST_F444W = F.response[355899:357351]
    F.NEWFIRM_J1_wave = F.wavelength[357351:357447]
    F.NEWFIRM_J1 = F.response[357351:357447]
    F.NEWFIRM_J2_wave = F.wavelength[357447:357526]
    F.NEWFIRM_J2 = F.response[357447:357526]
    F.NEWFIRM_J3_wave = F.wavelength[357526:357599]
    F.NEWFIRM_J3 = F.response[357526:357599]
    F.NEWFIRM_H1_wave = F.wavelength[357599:357669]
    F.NEWFIRM_H1 = F.response[357599:357669]
    F.NEWFIRM_H2_wave = F.wavelength[357669:357733]
    F.NEWFIRM_H2 = F.response[357669:357733]
    F.NEWFIRM_K_wave = F.wavelength[357733:357798]
    F.NEWFIRM_K = F.response[357733:357798]
    F.VIRCAM_Y_wave = F.wavelength[357798:357914]
    F.VIRCAM_Y = F.response[357798:357914]
    F.VIRCAM_J_wave = F.wavelength[357914:358062]
    F.VIRCAM_J = F.response[357914:358062]
    F.VIRCAM_H_wave = F.wavelength[358062:358286]
    F.VIRCAM_H = F.response[358062:358286]
    F.VIRCAM_K_wave = F.wavelength[358286:358545]
    F.VIRCAM_K = F.response[358286:358545]
    F.SuprimeCam_B_wave = F.wavelength[358545:358735]
    F.SuprimeCam_B = F.response[358545:358735]
    F.SuprimeCam_gplus_wave = F.wavelength[358735:358925]
    F.SuprimeCam_gplus = F.response[358735:358925]
    F.SuprimeCam_V_wave = F.wavelength[358925:359111]
    F.SuprimeCam_V = F.response[358925:359111]
    F.SuprimeCam_rplus_wave = F.wavelength[359111:359300]
    F.SuprimeCam_rplus = F.response[359111:359300]
    F.SuprimeCam_iplus_wave = F.wavelength[359300:359518]
    F.SuprimeCam_iplus = F.response[359300:359518]
    F.SuprimeCam_zplus_wave = F.wavelength[359518:359703]
    F.SuprimeCam_zplus = F.response[359518:359703]
    F.PanSTARRS1_g_wave = F.wavelength[359703:359882]
    F.PanSTARRS1_g = F.response[359703:359882]
    F.PanSTARRS1_r_wave = F.wavelength[359882:360069]
    F.PanSTARRS1_r = F.response[359882:360069]
    F.PanSTARRS1_i_wave = F.wavelength[360069:360250]
    F.PanSTARRS1_i = F.response[360069:360250]
    F.PanSTARRS1_z_wave = F.wavelength[360250:360418]
    F.PanSTARRS1_z = F.response[360250:360418]
    F.PanSTARRS1_y_wave = F.wavelength[360418:360624]
    F.PanSTARRS1_y = F.response[360418:360624]

    return F
#-----------------------------------------------------------------------------------	
def luminosity_distance(z, h0=73., omega_m=0.27, omega_l=0.73):
    
    '''
    Computing luminosity distance
    Input:  - z (float) -- redshift
            - h0 (float) (optional) -- hubble constant (in km/pc)
            - omega_m (float) (optional) -- matter density parameter
            - omega_l (float) (optional) -- dark energy density parameter
            
    Output: - luminosity distance (float) -- in parsec
    '''
    
    c = 2.9979e18 #velocity of lights
    omega_k = 1. - omega_m - omega_l
    dh = c/1.e13/h0 * 1.e6 #in pc
    
    if z > 0.:
        dc, edc = integrate.quad(lambda x: (omega_m * (1.+x)** 3 + omega_k * (1+x)**2 + omega_l)**(-.5), 0., z, epsrel=1e-4)
        dc = dh * dc
    else:
    # Bad idea as there is something *wrong* going on
        print('LumDist: z <= 0 -> Assume z = 0!')
        z = 0.
        #dist = 0.
        return 0
    
    if omega_k > 0.:
    	dm = dh * np.sinh(dc/dh * np.sqrt(omega_k)) / np.sqrt(omega_k)
    elif omega_k < 0.:
    	dm = dh * np.sin(dc/dh * np.sqrt(-omega_k)) / np.sqrt(-omega_k)
    else:
    	dm = dc
    return dm * (1+z)
#-----------------------------------------------------------------------------------	
def doppler_shift(wavelength, luminosity, z): #wavelength in micrometer
    '''
    Shifting intrinsic spectrum to observed spectrum using doppler shift formula
    Input:  - wavelength (list(float)) -- wavelength (in Angstorm)
            - luminosity (list(list(float))) -- intrinsic luminosity of each galaxy in each wavelength
            - z -- intrinsic redshift
            
    Output: - wavelength (list(float)) -- redshifted wavelength (in Angstorm)
            - luminosity (list(list(float))) -- observed luminosity (in erg/cm2/s/AA)
    '''
    
    pc2cm = 3.0856e18
    solar_lum = 3.839e33 # in cgs
    if z == 0:
    	distance = 10 * pc2cm #distance in cm: 1pc = 3.0856e18 cm
    else:
    	wavelength = wavelength * (1. + z)
    	distance = luminosity_distance(z) * pc2cm #distance in cm: 1pc = 3.0856e18 cm
    spectrum = luminosity * solar_lum / (4*np.pi*distance**2) #spec in erg/cm2/s/AA 
    return (wavelength, spectrum)
#-----------------------------------------------------------------------------------	
def compute_individual_mab(wavelength, luminosity, filt_wave, filt, z):
    
    '''
    Compute AB magnitude (mAB) using a single filter
    Input:  - wavelength (list(float)) -- intrinsic wavelength (in Angstorm)
            - luminosity (list(list(float))) -- intrinsic luminosity (in erg/cm2/s/AA)
            - filt_wave (str) -- name of the variable of filter wavelength (from the list)
            - filt (str) -- name of the filter (from the list)
            - z (float) -- redshift
            
    Output: - AB magnitude (float) of the input filter
    '''
    
    from scipy.integrate import simps
    c = 2.9979e18
    wavelength, spectrum = doppler_shift(wavelength, luminosity, z)
    filt_int  = np.interp(wavelength, filt_wave, filt)
    filtSpec = filt_int * spectrum
    flux = simps(filtSpec, wavelength)
    I1 = simps(spectrum*filt_int*wavelength,wavelength)
    I2 = simps(filt_int/wavelength, wavelength) 
    fnu = I1/I2/c
    mAB = -2.5*np.log10(fnu) - 48.6
    return(mAB)
#-----------------------------------------------------------------------------------	
def compute_mab(wavelength, luminosity, filter_list, z):
    '''
    Compute mab from a list of filters
    Input : - wavelength (list(float)) -- intrinsic wavelength (in Angstorm)
            - luminosity (list(list(float))) -- intrinsic luminosity (in erg/cm2/s/AA)
            - filter_list (list(str)) -- list of filter name
            - z (float) -- redshift
            
    Output: - AB magnitude (list(float)) -- computed AB magnitude of the input filters
    
    '''
    F = read_filters()
    mab_list = []
    for i in range(len(filter_list)):
    	filters_wave = eval('F.' + filter_list[i] + '_wave')
    	filters = eval('F.' + filter_list[i])
    	mab = compute_individual_mab(wavelength, luminosity, filters_wave, filters, z)
    	mab_list.append(mab)
    return(mab_list)
#-----------------------------------------------------------------------------------	
