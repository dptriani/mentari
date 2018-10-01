import numpy as np
Hubble_h = 0.73

def galdtype():
	# Define the data-type for the public version of SAGE
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
	Galdesc = np.dtype({'names':names, 'formats':formats}, align=True)
	return Galdesc

#-----------------------------------------------------------------------------------	

def read_one_file(name):
	
    Galdesc = galdtype()
    fin = open(name, 'rb')
    Ntrees = np.fromfile(fin,np.dtype(np.int32),1)[0]
    NtotGals = np.fromfile(fin,np.dtype(np.int32),1)[0]
    GalsPerTree = np.fromfile(fin, np.dtype((np.int32, Ntrees)),1)[0]
    G = []
    for i in range(Ntrees):
        G1 = np.fromfile(fin, Galdesc, GalsPerTree[i])
        G1 = G1.view(np.recarray)
        G.append(G1)
    return G

#-----------------------------------------------------------------------------------
    
def read_redshift_list(redshift, firstfile, lastfile, directory, filename):

    FirstSnap = 0
    LastSnap = len(redshift) - 1
    snapshot = list(range(len(redshift)-1, -1, -1))
    Galdesc = galdtype()
    GHist = []
    GalsPerTree = []

    for i in range(len(redshift)):
        GList = []
        GalsTree = []
        for j in range(firstfile, lastfile+1):
            name = (directory+filename+'_z'+f'{redshift[i]:.3f}'+'_'+f'{j}')
            G1 = read_one_file(name)
            GList.append(G1)
        GHist.append(GList)
        
    All_G = []
    for snap in range (FirstSnap, LastSnap+1):
        G_snap = []
        for filenr in range(len(GHist[snap])):
        #for filenr in ([0]):
            for tree in range(len(GHist[snap][filenr])):
            #for tree in ([0]):
                G_snap.extend(GHist[snap][filenr][tree])
        G_snap = np.array(G_snap)
        G_snap = G_snap.view(np.recarray)
        All_G.append(G_snap)
    
    return All_G
#-----------------------------------------------------------------------------------

def read_redshift_list_in(redshift, firstfile, lastfile, directory, filename):

    FirstSnap = 0
    LastSnap = len(redshift) - 1
    snapshot = list(range(len(redshift)-1, -1, -1))
    Galdesc = galdtype()
    GHist = []
    GalsPerTree = []

    for i in range(len(redshift)):
        GList = []
        GalsTree = []
        for j in range(firstfile, lastfile+1):
            name = (directory+filename+'_z'+f'{redshift[i]:.3f}'+'_'+f'{j}')
            G1 = read_one_file(name)
            GList.append(G1)
        GHist.append(GList)
    return GHist
    
#-----------------------------------------------------------------------------------
    
def mass_metal_history(redshift, firstfile, lastfile, directory, filename):

    FirstSnap = 0
    LastSnap = len(redshift) - 1
    snapshot = list(range(len(redshift)-1, -1, -1))
        
    GHist = read_redshift_list_in(redshift, firstfile, lastfile, directory, filename)
    ID_central_bulge = []
    ID_satellite_bulge = []
    ID_central_disk = []
    ID_satellite_disk = []
    ID0 = []
    

    All_G = []
    for snap in range (FirstSnap, LastSnap+1):
        G_snap = []
        for filenr in range(len(GHist[snap])):
        #for filenr in ([0]):
            for tree in range(len(GHist[snap][filenr])):
            #for tree in ([0]):
                G_snap.extend(GHist[snap][filenr][tree])
        G_snap = np.array(G_snap)
        G_snap = G_snap.view(np.recarray)
        All_G.append(G_snap)
    
    #Mass from SAGE
    Mass = []
    MassBulge = []
    MassDisk = []

    #Look for merger
    for filenr in range(len(GHist[LastSnap])):
        for tree in range(len(GHist[LastSnap][filenr])):
            w = np.where((GHist[LastSnap][filenr][tree].StellarMass * 1.e10 / Hubble_h) > 1.e8)[0]
            ID_z0 = GHist[LastSnap][filenr][tree].GalaxyIndex
            Mass.extend(GHist[LastSnap][filenr][tree].StellarMass * 1.e10)
            MassBulge.extend(GHist[LastSnap][filenr][tree].StellarMass * 1.e10)

            ID0.extend(ID_z0)
            ID_central_temp = list(ID_z0)
            ID_central = list(ID_z0)

            for snap in (snapshot):
                if GHist[snap][filenr][tree].mergeIntoID != []:
                    #major merger
                    w1 = np.where((GHist[snap][filenr][tree].mergeIntoID != -1) & (GHist[snap][filenr][tree].mergeType == 2))[0]
                    MergeSnap1 = GHist[snap][filenr][tree].mergeIntoSnapNum[w1]
                    MergeIndex1 = GHist[snap][filenr][tree].mergeIntoID[w1]
                    GalaxyID1 = GHist[snap][filenr][tree].GalaxyIndex[w1]

                    for i in range(len(MergeSnap1)):
                        if GHist[MergeSnap1[i]][filenr][tree].GalaxyIndex[MergeIndex1[i]] in ID_central_temp:
                            loc1 = np.where(ID_central_temp == GHist[MergeSnap1[i]][filenr][tree].GalaxyIndex[MergeIndex1[i]])[0]
                            ID_central_bulge.append(ID_central[loc1[0]])
                            ID_satellite_bulge.append(GalaxyID1[i])
                            ID_central_temp.append(GalaxyID1[i])
                            ID_central.append(ID_central[loc1[0]])

                    #minor_merger
                    w2 = np.where((GHist[snap][filenr][tree].mergeIntoID != -1) & (GHist[snap][filenr][tree].mergeType == 1))[0]
                    MergeSnap2 = GHist[snap][filenr][tree].mergeIntoSnapNum[w2]
                    MergeIndex2 = GHist[snap][filenr][tree].mergeIntoID[w2]
                    GalaxyID2 = GHist[snap][filenr][tree].GalaxyIndex[w2]

                    for i in range(len(MergeSnap2)):
                        if GHist[MergeSnap2[i]][filenr][tree].GalaxyIndex[MergeIndex2[i]] in ID_central_temp:
                            loc2 = np.where(ID_central_temp == GHist[MergeSnap2[i]][filenr][tree].GalaxyIndex[MergeIndex2[i]])[0]
                            ID_central_disk.append(ID_central[loc2[0]])
                            ID_satellite_disk.append(GalaxyID2[i])
                            ID_central_temp.append(GalaxyID2[i])
                            ID_central.append(ID_central[loc2[0]])

    index = np.array(ID0)

    #Mass from SAGE
    Mass = np.array(Mass)
    MassBulge = np.array(MassBulge)
    MassDisk = Mass - MassBulge
    MassHist = np.zeros((LastSnap+1, len(index)))
    MetalHist = np.zeros((LastSnap+1, len(index)))

    #Array for the computed mass
    mass = np.zeros(len(index))
    massbulge = np.zeros(len(index))
    metalbulge = np.zeros(len(index))
    massdisk = np.zeros(len(index))
    metaldisk = np.zeros(len(index))
    masshist = np.zeros((LastSnap+1, len(index)))
    metalhist = np.zeros((LastSnap+1, len(index)))
    bulgemasshist = np.zeros((LastSnap+1, len(index)))
    bulgemetalhist = np.zeros((LastSnap+1, len(index)))
    diskmasshist = np.zeros((LastSnap+1, len(index)))
    diskmetalhist = np.zeros((LastSnap+1, len(index)))

    #Set things to be array
    ID_central_bulge = np.array(ID_central_bulge)
    ID_satellite_bulge = np.array(ID_satellite_bulge)
    ID_central_disk = np.array(ID_central_disk)
    ID_satellite_disk = np.array(ID_satellite_disk)
    
    print("Constructing mass history of snapshot:")
    for snap in range(len(redshift)):
        print(snap, "/", len(redshift)-1)
        metal = np.zeros(len(index))
        allmass = np.zeros(len(index))
        if len(All_G[snap]) != 0:
            mass_loc, index_loc = np.nonzero(All_G[snap].GalaxyIndex[:,None] == index)
            massdisk[index_loc] = massdisk[index_loc] + All_G[snap].SfrDisk[mass_loc] * All_G[snap].dT[mass_loc] * 1e6
            massbulge[index_loc] = massbulge[index_loc] + All_G[snap].SfrBulge[mass_loc] * All_G[snap].dT[mass_loc] * 1e6
            metal[index_loc] = ((All_G[snap].SfrDiskZ[mass_loc] * All_G[snap].SfrDisk[mass_loc] * All_G[snap].dT[mass_loc] * 1e6) + (All_G[snap].SfrBulgeZ[mass_loc] * All_G[snap].SfrBulge[mass_loc] * All_G[snap].dT[mass_loc] * 1e6))
            allmass[index_loc] = (All_G[snap].SfrDisk[mass_loc] + All_G[snap].SfrBulge[mass_loc]) * All_G[snap].dT[mass_loc] * 1.e6

            #take care of major merger
            #--------------------------
            sfh_loc_b, gal_loc_b = np.nonzero(All_G[snap].GalaxyIndex[:,None] == ID_satellite_bulge)
            unq_ori_index_b = np.unique(ID_central_bulge[gal_loc_b])
            idx_loc_b = np.nonzero(unq_ori_index_b[:,None] == index)[1]
            temporary_ori_index_b = ID_central_bulge[gal_loc_b]

            #group the index based on the duplication in ori_index
            idx_sort_b = np.argsort(temporary_ori_index_b)
            sorted_array_b = temporary_ori_index_b[idx_sort_b]
            vals_b, idx_start_b, count_b = np.unique(sorted_array_b, return_counts=True, return_index=True)
            #sets of indices
            res_b = np.split(idx_sort_b, idx_start_b[1:])
            #filter them with respect to their size, keeping only items occurring more than once
            vals_b = vals_b[count_b > 1]
            res_b = list(filter(lambda x: x.size > 0, res_b)) #x.size shows number of minimum occurence

            temp_mass_b = [sum((All_G[snap].SfrDisk[sfh_loc_b[i]] + All_G[snap].SfrBulge[sfh_loc_b[i]]) * All_G[snap].dT[sfh_loc_b[i]] * 1.e6) for i in res_b]
            temp_metal_b = [sum((All_G[snap].SfrDiskZ[sfh_loc_b[j]] * All_G[snap].SfrDisk[sfh_loc_b[j]] + All_G[snap].SfrBulgeZ[sfh_loc_b[j]] * All_G[snap].SfrBulge[sfh_loc_b[j]]) * All_G[snap].dT[sfh_loc_b[j]] * 1.e6) for j in res_b]

            massbulge[idx_loc_b] = massbulge[idx_loc_b] + temp_mass_b
            allmass[idx_loc_b] = allmass[idx_loc_b] + temp_mass_b
            metal[idx_loc_b] = metal[idx_loc_b] + temp_metal_b

            #take care of minor merger
            #--------------------------
            sfh_loc_d, gal_loc_d = np.nonzero(All_G[snap].GalaxyIndex[:,None] == ID_satellite_disk)
            unq_ori_index_d = np.unique(ID_central_disk[gal_loc_d])
            idx_loc_d = np.nonzero(unq_ori_index_d[:,None] == index)[1]
            temporary_ori_index_d = ID_central_disk[gal_loc_d]

            #group the index based on the duplication in ori_index
            idx_sort_d = np.argsort(temporary_ori_index_d)
            sorted_array_d = temporary_ori_index_d[idx_sort_d]
            vals_d, idx_start_d, count_d = np.unique(sorted_array_d, return_counts=True, return_index=True)
            #sets of indices
            res_d = np.split(idx_sort_d, idx_start_d[1:])
            #filter them with respect to their size, keeping only items occurring more than once
            vals_d = vals_d[count_d > 1]
            res_d = list(filter(lambda x: x.size > 0, res_d)) #x.size shows number of minimum occurence

            temp_mass_d = [sum((All_G[snap].SfrDisk[sfh_loc_d[i]] + All_G[snap].SfrBulge[sfh_loc_d[i]])* All_G[snap].dT[sfh_loc_d[i]] * 1.e6) for i in res_d]
            temp_metal_d = [sum((All_G[snap].SfrDiskZ[sfh_loc_d[j]] * All_G[snap].SfrDisk[sfh_loc_d[j]] + All_G[snap].SfrBulgeZ[sfh_loc_d[j]] * All_G[snap].SfrBulge[sfh_loc_d[j]]) * All_G[snap].dT[sfh_loc_d[j]] * 1.e6) for j in res_d]

            massbulge[idx_loc_d] = massbulge[idx_loc_d] + temp_mass_d
            allmass[idx_loc_d] = allmass[idx_loc_d] + temp_mass_d
            metal[idx_loc_d] = metal[idx_loc_d] + temp_metal_d

            w_metal = np.where(allmass > 0)[0]
            metal[w_metal] = metal[w_metal] / allmass[w_metal]
            mass[index_loc] = massdisk[index_loc] + massbulge[index_loc]
            masshist[LastSnap-snap][index_loc] = mass[index_loc]
            diskmasshist[LastSnap-snap][index_loc] = massdisk[index_loc]
            bulgemasshist[LastSnap-snap][index_loc] = massbulge[index_loc]

            metalhist[LastSnap-snap][index_loc] = metal[index_loc]

    return masshist, metalhist
#=========================================================================================
    
def SED(lookbacktime, MassHist, MetalHist):

    f = open("files/bc2003_hr_m22_chab_ssp.ised_ASCII", "r")
    M1 = []
    for elem in f.read().split():
        try:
            M1.append(float(elem))
        except ValueError:
            pass
    f.close()

    f = open("files/bc2003_hr_m32_chab_ssp.ised_ASCII", "r")
    M2 = []
    for elem in f.read().split():
        try:
            M2.append(float(elem))
        except ValueError:
            pass
    f.close()

    f = open("files/bc2003_hr_m42_chab_ssp.ised_ASCII", "r")
    M3 = []
    for elem in f.read().split():
        try:
            M3.append(float(elem))
        except ValueError:
            pass
    f.close()

    f = open("files/bc2003_hr_m52_chab_ssp.ised_ASCII", "r")
    M4 = []
    for elem in f.read().split():
        try:
            M4.append(float(elem))
        except ValueError:
            pass
    f.close()

    f = open("files/bc2003_hr_m62_chab_ssp.ised_ASCII", "r")
    M5 = []
    for elem in f.read().split():
        try:
            M5.append(float(elem))
        except ValueError:
            pass
    f.close()

    f = open("files/bc2003_hr_m72_chab_ssp.ised_ASCII", "r")
    M6 = []
    for elem in f.read().split():
        try:
            M6.append(float(elem))
        except ValueError:
            pass
    f.close()

    lookback = M1[1:222]
    wavelength = M1[236: 7136]
    metallicity = [0.0001, 0.0004, 0.004, 0.008, 0.02, 0.05]
    ZZ = [M1, M2, M3, M4, M5, M6]
    time_grid = 221
    wave_grid = 6900
    lum = np.zeros((len(ZZ), time_grid, wave_grid))

    for j in range(len(metallicity)):
        M = ZZ[j]
        for i in range(time_grid):
            lum[j][i] = M[7137 + (wave_grid+54)*i: 7137 + (wave_grid+54)*i + 6900]

    gal_number = len(MassHist[0])
    computed_SFH = np.zeros((time_grid, gal_number))
    total_lum = np.zeros((gal_number, wave_grid))
    total_lum_0 = np.zeros((gal_number, wave_grid))
    new_mass_hist = np.zeros((time_grid, gal_number))
    new_metal_hist = np.zeros((time_grid, gal_number))

    for i in range(gal_number):
        new_mass_hist[:,i] = np.interp(lookback, lookbacktime, MassHist[:,i]) 
        new_metal_hist[:,i] = np.interp(lookback, lookbacktime, MetalHist[:,i]) 

    Z_ID = np.zeros_like(new_metal_hist)
    t1, m1 = np.where((new_metal_hist <= metallicity[0]) & (new_metal_hist < (metallicity[0] + metallicity[1])/2))
    t2, m2 = np.where((new_metal_hist > (metallicity[0] + metallicity[1])/2) & (new_metal_hist <= (metallicity[1] + metallicity[2])/2))
    t3, m3 = np.where((new_metal_hist > (metallicity[1] + metallicity[2])/2) & (new_metal_hist <= (metallicity[2] + metallicity[3])/2))
    t4, m4 = np.where((new_metal_hist > (metallicity[2] + metallicity[3])/2) & (new_metal_hist <= (metallicity[3] + metallicity[4])/2))
    t5, m5 = np.where((new_metal_hist > (metallicity[3] + metallicity[4])/2) & (new_metal_hist <= (metallicity[4] + metallicity[5])/2))
    t6, m6 = np.where(new_metal_hist > (metallicity[4] + metallicity[5])/2)

    Z_ID[t1, m1] = 0
    Z_ID[t2, m2] = 1
    Z_ID[t3, m3] = 2
    Z_ID[t4, m4] = 3
    Z_ID[t5, m5] = 4
    Z_ID[t6, m6] = 5

    n_metal_hist = np.zeros_like(new_metal_hist)
    n_metal_hist[t1, m1] = metallicity[0]
    n_metal_hist[t2, m2] = metallicity[1]
    n_metal_hist[t3, m3] = metallicity[2]
    n_metal_hist[t4, m4] = metallicity[3]
    n_metal_hist[t5, m5] = metallicity[4]
    n_metal_hist[t6, m6] = metallicity[5]
    
    print('Building SED')
    for i in range(len(lookback) - 1):
        print(i, "/", len(lookback) - 2)
        delta_mass = new_mass_hist[i] - new_mass_hist[i+1]
        deltamass = np.reshape(delta_mass, (-1, 1))

        temp_lum = [lum[int(j)][i] for j in Z_ID[i]]
        delta_lum = deltamass * temp_lum
        total_lum = total_lum + delta_lum

    return wavelength, total_lum

#=========================================================================================
    
def read_filters():
     
    F = type('', (), {})
    F.wavelength, F.response = np.loadtxt('filters/allfilters.dat', unpack=True)
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
     
    omega_k = 1. - omega_m - omega_l
    dh = c/1.e13/h0 * 1.e6 #in pc
    
    if z > 0.:
    	dc, edc = integrate.quad(lambda x: (omega_m * (1. + x) ** 3 + omega_k * (1 + x) ** 2 + omega_l) ** (-.5), 0., z, epsrel=1e-4)
    	dc = dh * dc
    else:
    # Bad idea as there is something *wrong* going on
    	print('LumDist: z <= 0 -> Assume z = 0!')
    	z = 0.
    	dist = 0.
    	
    if omega_k > 0.:
    	dm = dh * np.sinh(dc/dh * np.sqrt(omega_k)) / np.sqrt(omega_k)
    elif omega_k < 0.:
    	dm = dh * np.sin(dc/dh * np.sqrt(-omega_k)) / np.sqrt(-omega_k)
    else:
    	dm = dc
    return dm * (1+z)

#-----------------------------------------------------------------------------------
    
def doppler_shift(wavelength, luminosity, z): #wavelength in micrometer
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
    
def compute_mab(wavelength, luminosity, filt_wave, filt, z):
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
    
def mab(wavelength, luminosity, filter_list, z):
    
    F = read_filters()
    mab_list = []
    for i in range(len(filter_list)):
    	filters_wave = eval('F.' + filter_list[i] + '_wave')
    	filters = eval('F.' + filter_list[i])
    	mab = compute_mab(wavelength, luminosity, filters_wave, filters, z)
    	mab_list.append(mab)
    return(mab_list)
    
#=========================================================================================
