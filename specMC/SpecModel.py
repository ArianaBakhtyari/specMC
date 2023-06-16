import sys
from typing import Any
import pyspeckit
from astropy import units as u
from astropy.io import fits
from spectral_cube import SpectralCube
import matplotlib.pyplot as plt
import emcee
import numpy as np
from .models import ammonia, ammonia_fixfortho
from .fitters_withPrior import Specfit


class SpecModel:
    """
    This class creates the specModel object
    """
    def __init__(self, inFiles, fittype):
        """
        This function initializes the SpecModel object
        """
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
        """
        This function prompts the user for the initial points for the walker and its standard deviation
        """
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
        """ 
        This function asks if the user would like their function to be fit (Note: No fitting does not work)
        """
        Fit=None
        while Fit not in ("y","n"):
            Fit =input ("Would you like to fit your model? (y/n)")
        self.Fit=Fit

    def findComponents(self,guesses):
        """
        This function calculates the respective components and assigns the proper amount of titles for the number of components
        """
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
        """
        This function prompts the user for the prior values and their respective error
        """
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
        """
        This function prompts the user based on if the user wants to use the command line or if they choose to be prompted
        """
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
            self.sampleball.append(self.toAnArrayOfInt(sys.argv[len(arguments)-10])) #p0
            self.sampleball.append(self.toAnArrayOfInt(sys.argv[len(arguments)-9])) #std
            self.prior.append(self.toAnArrayOfInt(sys.argv[len(arguments)-8])) #prior
            self.prior.append(self.toAnArrayOfInt(sys.argv[len(arguments)-7])) #error
            self.Fit=sys.argv[len(arguments)-6] #fit
            self.x0=int(sys.argv[len(arguments)-5]) #x0
            self.y0=int(sys.argv[len(arguments)-4]) #y0
            self.Walkers=int(sys.argv[len(arguments)-3]) #walkers
            self.steps=int(sys.argv[len(arguments)-2]) #steps
            self.Burn=int(sys.argv[len(arguments)-1]) #burnin

    def toAnArrayOfInt(self, inputString):
        """
        This function returns an array when inputted a string
        """
        array=inputString.split(',')
        newarray = [float(n) for n in array]
        return newarray
    
    def askForPixel(self):
        """
        This function prompts the user for the desired pixel to be analyzed
        """
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
        """
        This function prompts the user asking if they'd like to select their own Burn-In value
        """
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
        """
        This function prompts the user for the number of walkers they'd like to use
        """
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
        """
        This function creates a cube to be analyzed
        """
        fn1=fits.open(inFile)
        cube1 = SpectralCube.read(fn1)
        cube1a = cube1.with_spectral_unit(u.Hz)
        return (cube1, cube1a)

    def convertUnits(self, cubea):
        """
        This function ensures the units are in the proper units
        """
        pcubea = pyspeckit.Cube(cube=cubea)
        pcubea.xarr.velocity_convention='radio'
        pcubea.xarr.convert_to_unit("km/s")
        print(pcubea)
        return cubea

    def get_sp(self, x, y, guesses):
        # need a way to
        pcube1a = pyspeckit.Cube(cube=self.cube1a)  # , xO=self.x0, yO=self.y0)
        pcube1a.xarr.velocity_convention = 'radio'
        pcube1a.xarr.convert_to_unit("km/s")
        pcube2a = pyspeckit.Cube(cube=self.cube2a)  # , xO=self.x0, yO=self.y0)
        pcube2a.xarr.velocity_convention = 'radio'
        pcube2a.xarr.convert_to_unit("km/s")
        cubes = pyspeckit.CubeStack([pcube1a, pcube2a])

        if self.fittype == 'fixfortho':
            print(self.fittype)
            fitter = ammonia_fixfortho.nh3_fixfortho_model(n_comp=self.ncomp)
        elif self.fittype == 'multiv':
            print(self.fittype)
            fitter = ammonia.nh3_multi_v_model_generator(n_comp=self.ncomp)

        self.sp = cubes
        self.sp.specfit.register_fitter(self.fittype, fitter, fitter.npars)

        # AttributeError: The 'specfit' object has no 'fitter' yet.  This means you haven't yet run a fit.
        # The fitter is not accessible until after a fit has been run.
        self.sp = self.sp.get_spectrum(x, y)

        #  replace the  pyspeckit Specfit object with the one fromt his package
        self.sp.specfit = Specfit(self.sp, self.sp.Registry)
        self.sp.specfit.register_fitter(self.fittype, fitter, fitter.npars)
        self.sp.specfit(fittype=self.fittype, guesses=guesses)


    def fit(self, fittype, guesses):
        """
        This function fits the model and creates the spectrum analyzed
        """
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

    def noFit(self,guesses):
        """
        This function runs if the user does not want to Fit their model (Note: this function breaks. Fit must always be True)
        """
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

    def runEmcee(self, array, progress=True):
        """
        This function runs emcee for the respective spectrum with the respective prior
        """
        self.ndim= len(array)
        self.nbins=self.ndim * 2
        self.nwalkers=self.Walkers
        self.nsteps=self.steps
        self.emcee_ensemble = Specfit.get_emcee(self.sp.specfit, self.proto_gauss_prior(), self.nwalkers)

        #HARD CODED EXAMPLE: 
        #self.p0 = emcee.utils.sample_ball((10, 5.3, 25,0.13, 8.16, 10, 5.3, 25, 0.13, 8.16),(3, 1, 2, 0.026, 0.051, 3, 1, 2, 0.026, 0.051), self.nbins*2)
        #self.p0 = emcee.utils.sample_ball(self.sampleball[0],self.sampleball[1], self.nbins) #commented out by mchen

        self.p0 = emcee.utils.sample_ball(self.sampleball[0],self.sampleball[1], self.nwalkers)
        print(self.p0.shape)
        print("Running emcee")
        self.emcee_ensemble.run_mcmc(self.p0, self.nsteps, progress=progress)

    def runEmcee_continue(self, nsteps, progress=True):
        '''
        Continue emcee sampling from where the last run left off
        Parameters
        ----------
        nsteps <int>
            the number of additional steps to run with the current emcee sampler
        '''
        print("Initial steps: {0}".format(self.emcee_ensemble.backend.iteration))
        self.emcee_ensemble.run_mcmc(None, nsteps, progress=progress)
        print("Final steps: {0}".format(self.emcee_ensemble.backend.iteration))
        # update nsteps
        self.nsteps = self.emcee_ensemble.backend.iteration

    def proto_gauss_prior(self):
        """
        This function sets up the respective mean and calculates the respective prior to be used in fitters_withPrior.py and model_withPrior.py
        """
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
        """
        This function removes the first respective number of steps identified
        """
        stb = steps_to_burn
        self.burned_chain = self.emcee_ensemble.chain[:,stb:,:] 
        self.burned_flatchain = np.reshape(self.burned_chain,(-1,self.ndim)) 
        self.burned_lnpost = self.emcee_ensemble.lnprobability[:,stb:] 

    def rejectLast(self,steps_to_burn):
        """
        This function removes the last respective number of steps identified
        """
        stb = steps_to_burn
        self.burned_chain = self.emcee_ensemble.chain[:,:stb,:] 
        self.burned_flatchain = np.reshape(self.burned_chain,(-1,self.ndim))
        self.burned_lnpost = self.emcee_ensemble.lnprobability[:,:stb] 
        