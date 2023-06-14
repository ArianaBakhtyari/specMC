# Bayesian Inference using MCMC for ammonia spectral model

#How to Setup By Cloning this Repo:

        1. pip setup.py install
        2. pip install -r requirements.txt
        *make sure you have atleast python version 3.6 installed
        3. python setup.py install

Example to run:
       1. Using the built in prompts:
                Run a two-component fixfotho fit:

                python runScript.py ../notebook/data/nh3_11_cubeA.fits fixfortho 20,10,20,0.11,8,20,10,20,0.11,8 ../notebook/data/nh3_22_cubeA.fits

                Run a one-component fixfortho fit:

                python runScript.py ../notebook/data/nh3_11_cubeA.fits fixfortho 20,10,20,0.11,8 ../notebook/data/nh3_22_cubeA.fits

        2. Using the command line
                python runScript.py ../notebook/data/nh3_11_cubeA.fits fixfortho 20,10,20,0.11,8,20,10,20,0.11,8 ../notebook/data/nh3_22_cubeA.fits 10,5.3,25,0.13,8.16,10,5.3,25,0.13,8.16 3,1,2,0.026,0.051,3,1,2,0.026,0.051 10,7,14,0.15,8.2,10,7,14,0.15,8.2 2,2,3,1,0.5,2,2,3,1,0.5 y 68 125 2000 5000 2500

Description of Files: 

Files Created:
 - runScript.py -> Runs the script for the fitting and model
 - plots.py -> Creates Subplots, Corner plots, and fitting plots for the specific fit
 - model.py -> Implements necessary fit, mcmc, burn-in, and prompts


 The two pyspeckit files:
- models.py - incorporates the multivariate normal prior of the model.
- fitters.py- Passes necessary parameters to the models.py file