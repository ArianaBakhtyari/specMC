
"""
This file creates all necessary plots to further aanalyze the data
"""
import corner
import numpy as np
import matplotlib.pyplot as plt

def plotCorner(emcee_ensemble, nameOfFile, titles, savename=None):
    """ 
    This function plots the original
    """
    #corner.corner(emcee_ensemble.flatchain);
    plt.figure(corner.corner(emcee_ensemble.flatchain, labels = titles))
    #savename = f'./Results/{nameOfFile}.png'
    if savename is not None:
        plt.savefig(savename)

def plotBurnCorner(flatchain, nameOfFile, titles, savename=None):
    """
    This function plots the corner plots with the corresponding burn-in values
    """
    #corner.corner(flatchain);

    plt.figure(corner.corner(flatchain, labels = titles))
    # f'./Results/{nameOfFile}.png'
    if savename is not None:
        plt.savefig(savename)


def plotMSP(emcee_ensemble,titles, nbins, savename=None):

        """
        This function plots the Median and Standard deviation histograms for each component
        """

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
        # "./Results/MedianStdPlots.png"
        if savename is not None:
                plt.savefig(savename)


def plotSubplots(emcee_ensemble, titles, nameOfFile, savename=None):
        """
        This function plots each components value with respect to each step
        """
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
        # f'./Results/{nameOfFile}.png'
        if savename is not None:
                plt.savefig(savename)


def plotBurnSubplots(emcee_ensemble, steps_to_burn, burn_setting, titles, nameOfFile, savename=None):
        """
        This function plots each components value with respect to the remaining steps after rejecting a certain number of steps
        """
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
        if savename is not None:
                plt.savefig(savename)
                #plt.savefig(f'./Results/{nameOfFile}.png')