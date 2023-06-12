#This file creates all the different plots to analyze the models in runScript.py
import corner
import numpy as np
import matplotlib.pyplot as plt

def plotCorner(emcee_ensemble, nameOfFile, titles):
    corner.corner(emcee_ensemble.flatchain);
    plt.figure(corner.corner(emcee_ensemble.flatchain, labels = titles))
    plt.savefig(f'./Results/{nameOfFile}.png')

def plotBurnCorner(flatchain, nameOfFile, titles):
    corner.corner(flatchain);
    plt.figure(corner.corner(flatchain, labels = titles))
    plt.savefig(f'./Results/{nameOfFile}.png')

def plotMSP(emcee_ensemble,titles, nbins):
        values=[]
        medVals=[]
        stdVals=[]
        for x in range(len(titles)):
                values.append(emcee_ensemble.flatchain[:,x])
        median_values = np.median(emcee_ensemble.flatchain,axis=0)
        medVals = median_values
        stddev_values = np.std(emcee_ensemble.flatchain,axis=0)
        stdVals = stddev_values
        for i,title in enumerate(titles):
                print(f'Median and stddev for {title} are:'+'{:.3f} +/- {:.3f}'.format(medVals[i],stdVals[i]))
        plt.figure(figsize=(15,7))
        for j in range(len(titles)):
                plt.subplot(1,6,j+1)
                plt.hist(values[j],bins=nbins)
                plt.title(f'{titles[j]} samples')
        plt.tight_layout();
        plt.savefig("./Results/MedianStdPlots.png")

def plotSubplots(emcee_ensemble, titles, nameOfFile):
        walker_to_inspect=0 # index of walker
        median_lnpost = np.median(emcee_ensemble.lnprobability[walker_to_inspect,:])
        plt.figure(figsize=(15,7))
        x=0
        for x, title in enumerate(titles):
                plt.subplot(2,8,x+1)
                plt.plot(emcee_ensemble.chain[walker_to_inspect,:,x])
                plt.title(f'{title}')
                plt.xlabel('Step number')
        plt.subplot(2,8,6)
        plt.plot(emcee_ensemble.lnprobability[walker_to_inspect,:])
        plt.axhline(median_lnpost,color='red', linestyle='--')
        plt.title('Ln Posterior')
        plt.xlabel('Step number')
        plt.tight_layout()
        plt.savefig(f'./Results/{nameOfFile}.png')


def plotBurnSubplots(emcee_ensemble, steps_to_burn, burn_setting, titles, nameOfFile):
        stb = steps_to_burn
        walker_to_inspect=0 # index of walker
        if burn_setting is "before":
                median_lnpost = np.median(emcee_ensemble.lnprobability[walker_to_inspect,stb:])
        else:  
                median_lnpost = np.median(emcee_ensemble.lnprobability[walker_to_inspect,:stb])
        plt.figure(figsize=(15,7))
        x=0
        for x, title in enumerate(titles):
                plt.subplot(2,8,x+1)
                if burn_setting is "before":
                        plt.plot(emcee_ensemble.chain[walker_to_inspect,stb:,x])
                else:
                       plt.plot(emcee_ensemble.chain[walker_to_inspect,:stb,x])
                plt.title(f'{title}')
                plt.xlabel('Step number')
        plt.subplot(2,8,6)
        if burn_setting is "before":
                plt.plot(emcee_ensemble.lnprobability[walker_to_inspect,stb:])
        else: 
               plt.plot(emcee_ensemble.lnprobability[walker_to_inspect,:stb])
        plt.axhline(median_lnpost,color='red', linestyle='--')
        plt.title('Ln Posterior')
        plt.xlabel('Step number')
        plt.tight_layout()
        plt.savefig(f'./Results/{nameOfFile}.png')