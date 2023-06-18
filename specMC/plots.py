
"""
This file creates all necessary plots to further analyze the data

The following plots can be create
    - corner plot- showing the one and two dimesnional projections of each component required for the fit
    - historgram for each component required to fit the model- a type of bar plot that groups the data into bins
    - subplots for each component required to fit the model- plotting each value of each required component

With each plot you have the ability to create a Burn-in version where you can reject a certain amount of steps defined by the user. The functions are plotBurnCorner() and plotBurnSubplots().
"""
import corner
import numpy as np
import matplotlib.pyplot as plt

def plotCorner(emcee_ensemble, titles, savename=None):
    """ 
    This function plots the original corner plot

    Parameters
    -----------
    emcee_ensemble:
        results of mcmc on the sample provided by the user
    nameOfFile: string
        desired name of plot
    titles: array
        array of strings of components required for the desired fit
    """
    #corner.corner(emcee_ensemble.flatchain);
    plt.figure(corner.corner(emcee_ensemble.flatchain, labels = titles))
    #savename = f'./Results/{nameOfFile}.png'
    if savename is not None:
        plt.savefig(savename)

def plotBurnCorner(flatchain, titles, savename=None):
    """
    This function plots the corner plots with the corresponding burn-in values. 

    Parameters

    -----------
    flatchain: np.ndarray
    nameOfFile: string
        desired name of plot
    titles: array
        array of strings of components required for the desired fit
    """
    #corner.corner(flatchain);

    plt.figure(corner.corner(flatchain, labels = titles))
    # f'./Results/{nameOfFile}.png'
    if savename is not None:
        plt.savefig(savename)


def plotMSP(emcee_ensemble,titles, nbins, savename=None):

        """
        This function plots the Median and Standard deviation histograms for each component required for the desired fit

        Parameters

        -----------
        emcee_ensemble:
            results of mcmc on the sample provided by the user
        titles: array
            array of strings of components required for the desired fit
        nbins: int
            number of bins
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


def plotSubplots(emcee_ensemble, titles, i_walker=0, savename=None):
        """
        This function plots each components value with respect to each step

        Parameters

        -----------
        emcee_ensemble:
            results of mcmc on the sample provided by the user
        titles: array
            array of strings of components for the desired fit
        nameOfFile: string
            desired name of plot
        i_walker:
                index of the walker to inspect
        """
        median_lnpost = np.median(emcee_ensemble.lnprobability[i_walker,:])
        plt.figure(figsize=(15,7))
        x=0
        for x, title in enumerate(titles):
                plt.subplot(2,8,x+1)
                plt.plot(emcee_ensemble.chain[i_walker,:,x])
                plt.title(f'{title}')
                plt.xlabel('Step number')
        plt.subplot(2,8,6)
        plt.plot(emcee_ensemble.lnprobability[i_walker,:])
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

        Parameters

        -----------
        emcee_ensemble:
            results of mcmc on the sample provided by the user
        steps_to_burn: int
            number of steps to reject
        burn_setting: string
                choice of "before", first # of steps, or "after", last # of steps, to be rejected
        titles: array
            array of strings of components for the desired fit
        nameOfFile: string
            desired name of plot
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