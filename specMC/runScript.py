import sys
import numpy as np
from SpecModel import SpecModel
from os.path import exists
from plots import *

class Run():
    print("Starting")
    print("creating model")
    try:
        if exists(sys.argv[4]):
            M=SpecModel([sys.argv[1], sys.argv[4]], sys.argv[2])
        elif exists(sys.argv[4]) == False:
            M=SpecModel([sys.argv[1]], sys.argv[2])
    except IndexError as e:
        M=SpecModel([sys.argv[1]], sys.argv[2])


    guess=sys.argv[3].split(',')
    guess = [float(n) for n in guess]
    guesses = np.array(guess)

    M.findComponents(guesses)
    M.promptUser(sys.argv)

    if M.Fit =="y":
        print("fitting spectrum")
        M.fit(sys.argv[2], guesses)
    else:
        print("setting up model")
        M.noFit(guesses)

    print("Setting up emcee")
    guesss=guess[:len(guesses)]
    M.runEmcee(guesss)

    print("Plotting graphs")
    plotMSP(M.emcee_ensemble,M.plotTitles, M.nbins)

    print("Plotting corner plots")
    plotCorner(M.emcee_ensemble, "firstCorner", M.cornerTitles)

    print("Plotting subplots")
    plotSubplots(M.emcee_ensemble, M.plotTitles, "SubPlots1")

    M.rejectFirst(M.Burn)
    plotBurnCorner(M.burned_flatchain, "Cornerpost-burn-in" , M.cornerTitles)
    plotBurnSubplots(M.emcee_ensemble, M.Burn,"before", M.plotTitles, "SubPlots2")

    M.rejectLast(M.Burn)
    plotBurnCorner(M.burned_flatchain, "Cornerpre-burn-in", M.cornerTitles)
    plotBurnSubplots(M.emcee_ensemble, M.Burn, "after", M.plotTitles, "SubPlots3")
