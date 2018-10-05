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
* `SEDpipeline.ipynb`: a tutorial on how to use the module.

## Functions in SED.py
This is a python module to compute SED and AB magnitude of galaxies from tabulated mass and metal histories, 
specifically from simulated galaxies from SAGE https://github.com/darrencroton/sage. 

### To use the module simply type:

`import SED` in a python script. 

### List of Functions:

#### `read_redshift_list(redshift, firstfile, lastfile, directory, filename)` : 
read the properties of all galaxies from the output of SAGE.
##### Input:
* `redshift` : a list of redshift of SAGE output.
* `firstfile` : first file of SAGE output.
* `lastfile` : last file of SAGE output.
* `directory` : specific directory of SAGE output.
* `filename` : specific filename of SAGE output, default: model.
##### Output:
* `G[redshift][properties]` : a 2 dimension array containing the properties of galaxies in each redshift listed. 

#### mass_metal_history(redshift, firstfile, lastfile, directory, filename) :
build mass history and metal history from the output of SAGE.

##### Input:
* `redshift` : a list of redshift of SAGE output.
* `firstfile` : first file of SAGE output.
* `lastfile` : last file of SAGE output.
* `directory` : specific directory of SAGE output.
* `filename` : specific filename of SAGE output, default: model.

##### Output:
* `mass` : mass history of each galaxies in each redshift listed. The order of the redshift will be reversed (the first element of this array is for the last redshift value in redshift list).
* `metal` : metal history of each galaxies in each redshift listed. The order of the redshift will be reversed (the first element of this array is for the last redshift value in redshift list).


