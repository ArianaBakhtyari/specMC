# Bayesian Inference using MCMC for ammonia spectral model

## How to setup by Cloning this Repo:

        1. git clone {github address} 
        2. pip install -r requirements.txt
        (make sure you have atleast python version 3.6 installed)
        3. python setup.py install

## General Instructions
There are two ways to run the code 
        - Using the built in prompts
        - Using the command line
When running the script you must run _python runScript.py {firstDataFile} {fittype} {guesses}_ for the most general fit. If you'd like to compare two Datasets then you have the option of adding the _{secondDataFile}_. You also have the choice between a one-component and a two-component fit. When using two-components you must enter double the guesses, p0,  standard deciation, prior, error required for the fit. 

The following is a description of each variable:

- fittype- The fittype you'd like for your model
- guesses- Initial guesses for your fit
- p0- Starting points for your walkers
- standard deviation of p0- std for your walkers
- prior- The prior for each component
- error of prior- The error for each component
- y or n for Fitting the model- If you'd like to fit the model
- x0- The x coordinate of the pixel you'd like to analyze
- y0- The y coordinate of the pixel you'd like to analyze
- number of walkers- Number of walkers for mcmc
- number of steps- Number of steps you'd like each walker to take
- number of steps to reject (Burn)- Number of steps you'd like to reject when analayzing the data (The initial steps are always noisy and inaccurate)

## Example to run:
###     1. Using the built in prompts:
                Run a one-component fixfortho fit (one dataset):
_python runScript.py {firstDataFile} {fittype} {guesses}_

`                python runScript.py ../notebook/data/nh3_11_cubeA.fits fixfortho 20,10,20,0.11,8`

                Run a two-component fixfortho fit (two datasets):
_python runScript.py {firstDataFile} {fittype} {guesses x2} {secondDataFile}_
 
`                python runScript.py ../notebook/data/nh3_11_cubeA.fits fixfortho 20,10,20,0.11,8,20,10,20,0.11,8 ../notebook/data/nh3_22_cubeA.fits`

                Run a one-component fixfortho fit (two datasets):
_python runScript.py {firstDataFile} {fittype} {guesses} {secondDataFile}_
 
`                python runScript.py ../notebook/data/nh3_11_cubeA.fits fixfortho 20,10,20,0.11,8 ../notebook/data/nh3_22_cubeA.fits`

###     2. Using the command line
        Run a one-component fixfortho fit (one dataset):
_python runScript.py {firstDataFile} {fittype} {guesses} {p0} {standard deviation of p0} {prior} {error of prior} {y or n for Fitting the model} {x0} {y0} {number of walkers} {number of steps} {number of steps to reject}_

`        python runScript.py ../notebook/data/nh3_11_cubeA.fits fixfortho 20,10,20,0.11,8 10,5.3,25,0.13,8.16 3,1,2,0.026,0.051 10,7,14,0.15,8.2 2,2,3,1,0.5 y 68 125 2000 5000 2500`

        Run a two-component fixfortho fit (two datasets):

_python runScript.py {firstDataFile} {fittype} {guesses} {secondDataFile} {p0 x2} {standard deviation of p0 x2} {prior x2} {error of prior x2} {y or n for Fitting the model} {x0} {y0} {number of walkers} {number of steps} {number of steps to reject}_

`                python runScript.py ../notebook/data/nh3_11_cubeA.fits fixfortho 20,10,20,0.11,8,20,10,20,0.11,8 ../notebook/data/nh3_22_cubeA.fits 10,5.3,25,0.13,8.16,10,5.3,25,0.13,8.16 3,1,2,0.026,0.051,3,1,2,0.026,0.051 10,7,14,0.15,8.2,10,7,14,0.15,8.2 2,2,3,1,0.5,2,2,3,1,0.5 y 68 125 2000 5000 2500 `

        Run a one-component fixfortho fit (two datasets):
_python runScript.py {firstDataFile} {fittype} {guesses} {secondDataFile} {p0} {standard deviation of p0} {prior} {error of prior} {y or n for Fitting the model} {x0} {y0} {number of walkers} {number of steps} {number of steps to reject}_

`        python runScript.py ../notebook/data/nh3_11_cubeA.fits fixfortho 20,10,20,0.11,8 ../notebook/data/nh3_22_cubeA.fits 10,5.3,25,0.13,8.16 3,1,2,0.026,0.051 10,7,14,0.15,8.2 2,2,3,1,0.5 y 68 125 2000 5000 2500`

## Description of Files: 

### Files Created:
 - runScript.py -> Runs the script for the fitting and model
 - plots.py -> Creates Subplots, Corner plots, and fitting plots for the specific fit
 - model.py -> Implements necessary fit, mcmc, burn-in, and prompts


### The two pyspeckit files:
- model_withPrior.py - incorporates the multivariate normal prior of the model.
- fitter_withPrior.py- Passes necessary parameters to the model_withPrior.py file

## How to read a .h5 file

import h5py
filename = "file.hdf5"

with h5py.File(filename, "r") as f:
    # Print all root level object names (aka keys) 
    # these can be group or dataset names 
    print("Keys: %s" % f.keys())
    # get first object name/key; may or may NOT be a group
    a_group_key = list(f.keys())[0]

    # get the object type for a_group_key: usually group or dataset
    print(type(f[a_group_key])) 

    # If a_group_key is a group name, 
    # this gets the object names in the group and returns as a list
    data = list(f[a_group_key])

    # If a_group_key is a dataset name, 
    # this gets the dataset values and returns as a list
    data = list(f[a_group_key])
    # preferred methods to get dataset values:
    ds_obj = f[a_group_key]      # returns as a h5py dataset object
    ds_arr = f[a_group_key][()]  # returns as a numpy array
    