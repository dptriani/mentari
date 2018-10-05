# sage-sed
## Installation
If you have git installed, sage-sed can be obtained with the following commands: 
```
cd /path/to/desired/location/
git clone https://github.com/dptriani/sage-sed
```

Otherwise download it as a zip file [here](https://github.com/dptriani/sage-sed/archive/master.zip) and then unzip it. 

## Contents
* `files`: contains the filters wavelength and response file obtained from [FSPS](https://github.com/cconroy20/fsps/blob/master/data/allfilters.dat)
and the SSP tables from [Bruzual and Charlot (2003)](http://www.bruzual.org/bc03/).
* `list of filters`: list of filters provided. When specifying filters to compute the AB magnitude, please refer to the
name format in this list.
* `SED.py`: the main module.
* `tutorial.ipynb`: a tutorial on how to use the module.

## Functions in SED.py
This is a python module to compute SED and AB magnitude of galaxies from tabulated mass and metal histories, 
specifically from simulated galaxies from SAGE https://github.com/darrencroton/sage. 

### To use the module simply type:

`import SED` in a python script. 

### List of Functions:

#### * `read_redshift_list(redshift, firstfile, lastfile, directory, filename)` : 
read the properties of all galaxies from the output of SAGE
##### Input:
* `redshift` : a list of redshift of SAGE output.
* `firstfile` : first file of SAGE output
* `lastfile` : last file of SAGE output
* `directory` : specific directory of SAGE output
* `filename` : specific filename of SAGE output, default: model
##### Output:
* `G` : a `m * n` dimension array containing the properties of galaxies. `m` is the number of redshift listed, `n is the number of galaxies in that redshift`

#### * `mass_metal_history(redshift, firstfile, lastfile, directory, filename)` :
build mass history and metal history from the output of SAGE

##### Input:
* `redshift` : a list of redshift of SAGE output
* `firstfile` : first file of SAGE output
* `lastfile` : last file of SAGE output
* `directory` : specific directory of SAGE output
* `filename` : specific filename of SAGE output, default: model

##### Output:
* `mass` : mass of each galaxies in each redshift listed in 10^10 Msun/h. The order of the redshift will be reversed (the first element of this array is for the last redshift value in redshift list)
* `metal` : metallicity of each galaxies in each redshift listed. The order of the redshift will be reversed (the first element of this array is for the last redshift value in redshift list)

#### * `SED(lookbacktime, mass, metal)`:
building SED from a tabulated mass and metal history of galaxy(es) in each lookback time
##### input :
* `lookbacktime` : a list lookback time in year
* `mass` : an array containing mass(es) of galaxy(es) in Msun
* `metal` : metallicities of galaxy(es)
#### output:
* `wavelength` : a list of wavelength in Angstorm
* `spectra` :  luminosity of galaxy(es) in each wavelength in Lsun/Angstorm

#### * `mab(wavelength, spectra, filter_list, z)`
##### input :
* `wavelength` : list of wavelength in Angstorm
* `spectra` : luminosity of galaxy(es) in each wavelength in Lsun/Angstorm
* `filter list` : a list containing the name of filters (please refer to `list of filters` for the naming convention)
* `z` : a redshift (0 for absolute magnitude, other for apparent)
##### output :
* `mab` : a list of magnitudes in each filters
