#This manipulates the data in runScript.py
import os
import sys
from typing import Any
#sys.path.append(f"{os.getenv('STEM')}/pyspeckit-IncorporatePrior")
import pyspeckit
from astropy import units as u
from astropy.io import fits
from spectral_cube import SpectralCube
import matplotlib.pyplot as plt
import emcee
import numpy as np
from models import ammonia, ammonia_fixfortho
from fitters_withPrior import Specfit 


class SpecModel:
    
    def __init__(self, inFiles, fittype):
        self.fittype= fittype
        if fittype=='gaussian':
            self.titles=['Amplitude', 'Shift', 'Width']
            self.plotTitles=self.titles
        if fittype=='ammonia':
            self.titles=['RotationTemp', 'ExcitTemp', 'Column Density', 'Line Width', 'Line offset', 'Fortho']
            self.plotTitles=self.titles
        if fittype=='cold_ammonia':
            self.titles=['Kinetic Temperature']
            self.plotTitles=self.titles
        if fittype=='multiv':
            self.titles=['trot', 'tex', 'ntot', 'width', 'xoff_v', 'fortho']
            self.plotTitles=['$T_{Rot}$', '$T_{Ex}$', '$N_{Tot}$', 'Width', '$X_{offV}$', 'Fortho']
        if fittype=='fixfortho':
            self.titles=['trot', 'tex', 'ntot', 'width', 'xoff_v']
            self.plotTitles=['$T_{Rot}$', '$T_{Ex}$', '$N_{Tot}$', 'Width', '$X_{offV}$']
        if len(inFiles) ==1:
            self.cube1, self.cube1a =self.createCube(inFiles[0])
            self.model=1
        elif len(inFiles) ==2:
            self.cube1, self.cube1a =self.createCube(inFiles[0])
            self.cube2, self.cube2a =self.createCube(inFiles[1])
            self.model=2
        
    def getSampleBall(self):
        sampleball=[]
        p0=[]
        std=[]
        for j in range(0,len(self.titles)):
            while True:
                try:
                    p=float(input("Please enter the {} value ".format((self.titles[j]))))
                    break
                except ValueError:
                    print("ERROR: Please enter a valid FLOAT")
                    pass
            while True:
                try:
                    s=float(input("and the standard deviation "))
                    break
                except ValueError:
                    print("ERROR: Please enter a valid FLOAT")
                    pass
            p0.append(p)
            std.append(s)
        sampleball.append(p0)
        sampleball.append(std)
        self.sampleball=sampleball
    
    def askForFit(self):
        Fit=None
        while Fit not in ("y","n"):
            Fit =input ("Would you like to fit your model? (y/n)")
        self.Fit=Fit

    def findComponents(self,guesses):
        self.ncomp=int(len(guesses)/len(self.titles))
        i=0
        actualTitles=[]
        while i is not self.ncomp:
            x=0
            for x in range(0, len(self.titles)):
                actualTitles.append(self.titles[x])
                x=x+1
            i=i+1
        self.titles=actualTitles
        self.cornerTitles=[]
        if self.ncomp != 1:
            component=1
            while component <= self.ncomp:
                j=0
                for j in range(0, len(self.plotTitles)):
                    self.cornerTitles.append('{}_{}'.format(self.plotTitles[j], component))
                    j=j+1
                component=component +1
            
    def askForPrior(self):
        thePrior=[]
        prior=[]
        error=[]
        for j in range(0,len(self.titles)):
            while True:
                try:
                    p=float(input("Please enter the {} PRIOR value ".format(self.titles[j])))
                    break
                except ValueError:
                    print("ERROR. Please enter a valid FLOAT ")
                    pass
            while True:
                try:
                    e=float(input("and the respective error "))
                    break
                except ValueError:
                    print("ERROR. Please enter a valid FLOAT ")
                    pass
            prior.append(p)
            error.append(e)
        thePrior.append(prior)
        thePrior.append(error)
        self.prior=thePrior

    def promptUser(self, arguments):
        if len(arguments) <= 5:
            self.getSampleBall()
            self.askForPrior()
            self.askForFit()
            self.askForPixel()
            self.askForWalkers()
            self.askForBurnIn()
        elif len(arguments) > 5:
            self.sampleball=[]
            self.prior=[]    
            guess=sys.argv[3].split(',')
            guess = [float(n) for n in guess]
            self.sampleball.append(self.toAnArrayOfInt(sys.argv[5])) #p0
            self.sampleball.append(self.toAnArrayOfInt(sys.argv[6])) #std
            self.prior.append(self.toAnArrayOfInt(sys.argv[7])) #prior
            self.prior.append(self.toAnArrayOfInt(sys.argv[8])) #error
            self.Fit=sys.argv[9] #fit
            self.x0=int(sys.argv[10]) #x0
            self.y0=int(sys.argv[11]) #y0
            self.Walkers=int(sys.argv[12]) #walkers
            self.steps=int(sys.argv[13]) #steps
            self.Burn=int(sys.argv[14]) #burnin

    def toAnArrayOfInt(self, inputString):
        array=inputString.split(',')
        newarray = [float(n) for n in array]
        return newarray
    
    def askForPixel(self):
        while True:
            try:
                x0=int(input("Please enter the x0 value "))
                break
            except ValueError:
                print("ERROR. Please enter a valid INTEGER ")
                pass
        while True:
            try:
                y0=int(input("and the respective y0 "))
                break
            except ValueError:
                print("ERROR. Please enter a valid INTEGER ")
                pass
        self.x0=x0
        self.y0=y0

    def askForBurnIn(self):
        Burn = "m"
        while Burn not in ("y","n"):
            Burn =input ("Would you like to select how many points to Burn-In? (y/n)")
        if Burn == "y":
            while True:
                try:
                    burnIn=int(input("How many steps would you like to reject?"))
                    break
                #add if burnIn is more than walker numbers
                except ValueError:
                    print("ERROR. Please enter a valid INTEGER ")
                    pass
        else:
            burnIn=2500
        self.Burn=burnIn
        
    def askForWalkers(self):
        while True:
            try:
                walkers=int(input("How many walkers would you like: "))
                break
            except ValueError:
                print("ERROR. Please enter a valid INTEGER ")
                pass
        while True:
            try:
                steps=int(input("How many steps would you like: "))
                break
            except ValueError:
                print("ERROR. Please enter a valid INTEGER ")
                pass
        self.Walkers=walkers
        self.steps=steps
                
    def createCube(self, inFile):
        fn1=fits.open(inFile)
        cube1 = SpectralCube.read(fn1)
        cube1a = cube1.with_spectral_unit(u.Hz)
        return (cube1, cube1a)

    def convertUnits(self, cubea):
        pcubea = pyspeckit.Cube(cube=cubea)
        pcubea.xarr.velocity_convention='radio'
        pcubea.xarr.convert_to_unit("km/s")
        print(pcubea)
        return cubea

    def fit(self, fittype, guesses):
        self.guesses = guesses
        if self.model ==1:
            pcube1 = pyspeckit.Cube(cube=self.cube1a, xO=self.x0, yO=self.y0)
            pcube1.xarr.velocity_convention='radio'
            pcube1.xarr.convert_to_unit("km/s")
            fitter= ammonia.nh3_multi_v_model_generator(n_comp=self.ncomp)
            fitter1= ammonia_fixfortho.nh3_fixfortho_model(n_comp =self.ncomp)
            sp = pyspeckit.Spectrum(data=pcube1.data, xarr=pcube1.xarr, xarrkwargs={'units':'km/s'})
            sp.specfit.Registry.add_fitter('multiv',fitter, fitter.npars)
            sp.specfit.Registry.add_fitter('fixfortho',fitter1, fitter1.npars)
            sp.specfit(fittype=fittype, guesses=guesses)
            sp.plotter()
            sp.specfit.plot_fit()
            self.sp= sp
            print(sp.specfit.plot_fit())
        else :
            pcube1a = pyspeckit.Cube(cube=self.cube1a, xO=self.x0, yO=self.y0)
            pcube1a.xarr.velocity_convention='radio'
            pcube1a.xarr.convert_to_unit("km/s")
            pcube2a = pyspeckit.Cube(cube=self.cube2a, xO=self.x0, yO=self.y0)
            pcube2a.xarr.velocity_convention='radio'
            pcube2a.xarr.convert_to_unit("km/s")
            fitter= ammonia.nh3_multi_v_model_generator(n_comp=self.ncomp)
            fitter1= ammonia_fixfortho.nh3_fixfortho_model(n_comp=self.ncomp)
            cubes = pyspeckit.CubeStack([pcube1a, pcube2a])
            cubes.specfit.Registry.add_fitter('multiv',fitter, fitter.npars)
            cubes.specfit.Registry.add_fitter('fixfortho',fitter1, fitter1.npars)
            self.sp=cubes
            cubes.fiteach(fittype=fittype, guesses=guesses)
            cubes.plotter()
            cubes.specfit.plot_fit()
        plt.savefig(f'./Results/{fittype}plot.png')

#this function breaks. Fit must always be True
    def noFit(self,guesses):
        self.guesses = guesses
        if self.model==1:
            pcube1 = pyspeckit.Cube(cube=self.cube1a, x0= self.x0, y0=self.y0)
            pcube1.xarr.velocity_convention='radio'
            pcube1.xarr.convert_to_unit("km/s")
            sp = pyspeckit.Spectrum(data=pcube1.data, xarr=pcube1.xarr, xarrkwargs={'units':'km/s'})
            self.sp=sp
        else: 
            pcube1a = pyspeckit.Cube(cube=self.cube1a, x0=self.x0, y0=self.y0)
            pcube1a.xarr.velocity_convention='radio'
            pcube1a.xarr.convert_to_unit("km/s")
            pcube2a = pyspeckit.Cube(cube=self.cube2a, x0=self.x0, y0=self.y0)
            pcube2a.xarr.velocity_convention='radio'
            pcube2a.xarr.convert_to_unit("km/s")
            cubes = pyspeckit.CubeStack([pcube1a, pcube2a])
            self.sp=cubes

    def runEmcee(self, array):
        self.ndim= len(array)
        self.nbins=self.ndim *self.ncomp
        self.nwalkers=self.Walkers
        self.nsteps=self.steps
        self.emcee_ensemble= Specfit.get_emcee(self.sp.specfit, self.proto_gauss_prior)
        #self.emcee_ensemble= self.sp.specfit.get_emcee(self.proto_gauss_prior()) 
        #HARD CODED EXAMPLE: 
        #self.p0 = emcee.utils.sample_ball((10, 5.3, 25,0.13, 8.16, 10, 5.3, 25, 0.13, 8.16),(3, 1, 2, 0.026, 0.051, 3, 1, 2, 0.026, 0.051), self.nbins*2)
        self.p0 = emcee.utils.sample_ball(self.sampleball[0],self.sampleball[1], self.nbins)
        print(self.p0.shape)
        print("Running emcee")
        self.emcee_ensemble.run_mcmc(self.p0,self.nsteps)

    def proto_gauss_prior(self):
        self.priorvals={}
        self.priorvals["mean"]=[]
        self.priorvals["cov"]=[]
        covariance =[]
        for prior in self.prior[1]:
            sigma=prior
            covariance.append(sigma **2)
        self.priorvals["mean"]=self.prior[0]
        self.priorvals["cov"]=covariance
        self.priorvals["steps"]=0
        return self.priorvals

    def rejectFirst(self,steps_to_burn):
        stb = steps_to_burn
        self.burned_chain = self.emcee_ensemble.chain[:,stb:,:] 
        self.burned_flatchain = np.reshape(self.burned_chain,(-1,self.ndim)) 
        self.burned_lnpost = self.emcee_ensemble.lnprobability[:,stb:] 

    def rejectLast(self,steps_to_burn):
        stb = steps_to_burn
        self.burned_chain = self.emcee_ensemble.chain[:,:stb,:] 
        self.burned_flatchain = np.reshape(self.burned_chain,(-1,self.ndim))
        self.burned_lnpost = self.emcee_ensemble.lnprobability[:,:stb] 
        