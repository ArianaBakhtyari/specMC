"""
========================================
Multi-component Ammonia inversion transition TROT fitter
========================================

Ammonia inversion transition TROT fitter built from Adam Ginsburg's pyspeckit, which was translated from Erik
Rosolowsky's https://github.com/low-sky/nh3fit

.. moduleauthor:: Mike Chen <chen.m@queensu.ca>

Module API
^^^^^^^^^^

"""
from __future__ import division

import numpy as np

from pyspeckit.spectrum.models import model
from pyspeckit.spectrum.models.ammonia_constants import (line_names, freq_dict, aval_dict, ortho_dict,
                                                         voff_lines_dict, tau_wts_dict)
from pyspeckit.spectrum.models.ammonia_constants import (ckms, h, kb)



#from ...mpfit import mpfit
#from ...spectrum.parinfo import ParinfoList,Parinfo
#from . import fitter
#from . import model
import matplotlib.cbook as mpcb
import copy
from astropy import log
from six import iteritems
#from . import mpfit_messages
import operator
import string
import warnings

#from .ammonia_constants import (line_names, freq_dict, aval_dict, ortho_dict,
#                                voff_lines_dict, tau_wts_dict)

TCMB = 2.7315 # K


def ammonia(xarr, trot=20, tex=None, ntot=14, width=1, xoff_v=0.0, fortho=0.0,
            tau=None, fillingfraction=None, return_tau=False,
            return_tau_profile=False, background_ta=TCMB, verbose=False,
            return_components=False, debug=False, line_names=line_names,
            ignore_neg_models=False, bg_is_baseline=False):
    
    # Convert X-units to frequency in GHz
    if xarr.unit.to_string() != 'GHz':
        xarr = xarr.as_unit('GHz')
    
    CMB = T_antenna(TCMB, xarr.value)
    
    model = ammonia_slab(xarr, trot=trot, tex=tex, ntot=ntot, width=width, xoff_v=xoff_v, fortho=fortho,
                         tau=tau, fillingfraction=fillingfraction, return_tau=return_tau,
                         background_ta=CMB,
                         verbose=verbose, debug=debug, line_names=line_names, bg_is_baseline=bg_is_baseline)
    return model - CMB



def ammonia_multi_v(xarr, ncomp, *parameters, background_tb=TCMB, **kwargs):
    """
    :param xarr:
    :param ncomp:
    :param parameters:
        a list of trot, tex, ntot, width, xoff_v, fortho for each component
    :param background_tb:
    :param kwargs:
    :return:
    """

    # Convert X-units to frequency in GHz
    if xarr.unit.to_string() != 'GHz':
        xarr = xarr.as_unit('GHz')

    args = parameters
    # need to sanity check the number of parameters is consistent with the number of components
    if len(args)/ncomp != 6:
        print("[ERROR]: the number of parameters ({}) does not match the number of components ({})".format(len(args)                                                                                      ,ncomp))
        return None

    # Start with a constant background radiation field (e.g., Cosmic Microwave Background by default)
    BGRF = T_antenna(background_tb, xarr.value)
    model = BGRF
    #model = ammonia_slab(xarr, background_ta=BGRF, **kwargs)

    #args = parameters
    # iteratively move through each slabs towards the observer (i.e., radiative transfer)
    for trot, tex, ntot, width, xoff_v, fortho in zip(args[::6], args[1::6], args[2::6],
                                                      args[3::6], args[4::6], args[5::6]):
        model = ammonia_slab(xarr, trot=trot, tex=tex, ntot=ntot, width=width, xoff_v=xoff_v, fortho=fortho,
                                     background_ta=model, **kwargs)

    return model - BGRF



def ammonia_slab(xarr, trot=20, tex=None, ntot=14, width=1, xoff_v=0.0, fortho=0.0,
                 tau=None, fillingfraction=None, return_tau=False,
                 background_ta=TCMB,
                 verbose=False, debug=False, line_names=line_names, bg_is_baseline=False):
    """
    Generate a model Ammonia spectrum based on input temperatures, column, and
    gaussian parameters.  The returned model will be in Kelvin (brightness
    temperature) units.

    Note that astropy units are not used internally for performance reasons.  A
    wrapped version of this module including those units would be a good idea,
    as it is definitely possible to implement this with unit support and good
    performance.

    Parameters
    ----------
    xarr: `pyspeckit.spectrum.units.SpectroscopicAxis`
        Array of wavelength/frequency values
    trot: float
        The rotational temperature of the lines.  This is the excitation
        temperature that governs the relative populations of the rotational
        states.
    tex: float or None
        Excitation temperature. Assumed LTE if unspecified (``None``) or if
        tex>trot.  This is the excitation temperature for *all* of the modeled
        lines, which means we are explicitly assuming T_ex is the same for all
        lines.
    ntot: float
        Total log column density of NH3.  Can be specified as a float in the
        range 5-25
    width: float
        Line width (Gaussian sigma) in km/s
    xoff_v: float
        Line offset in km/s
    fortho: float
        Fraction of NH3 molecules in ortho state.  Default assumes all para
        (fortho=0).
    tau: None or float
        If tau (optical depth in the 1-1 line) is specified, ntot is NOT fit
        but is set to a fixed value.  The optical depths of the other lines are
        fixed relative to tau_oneone
    fillingfraction: None or float
        fillingfraction is an arbitrary scaling factor to apply to the model
    return_tau: bool
        Return a dictionary of the optical depths in each line instead of a
        synthetic spectrum
    return_tau_profile: bool
        Return a dictionary of the optical depth profiles in each line, i.e.,
        the optical depths that will be used in conjunction with T_ex to produce
        the synthetic spectrum
    return_components: bool
        Return a list of arrays, one for each hyperfine component, instead of
        just one array
    background_tb : float
        The background brightness temperature.  Defaults to TCMB.
    ignore_neg_models: bool
        Normally if background=TCMB and the model is negative, an exception
        will be raised.  This parameter will simply skip that exception.  Use
        with extreme caution: negative models (absorption spectra against the
        CMB) are not physical!  You may want to allow this in some cases
        because there can be numerical issues where the model goes negative
        when it shouldn't.
    verbose: bool
        More messages
    debug: bool
        For debugging.

    Returns
    -------
    spectrum: `numpy.ndarray`
        Synthetic spectrum with same shape as ``xarr``
    component_list: list
        List of `numpy.ndarray`'s, one for each hyperfine component
    tau_dict: dict
        Dictionary of optical depth values for the various lines
        (if ``return_tau`` is set)
    """

    from pyspeckit.spectrum.models.ammonia_constants import (ckms, ccms, h, kb,
                                                             Jortho, Jpara, Brot, Crot)

    # Convert X-units to frequency in GHz
    if xarr.unit.to_string() != 'GHz':
        xarr = xarr.as_unit('GHz')

    if tex is None:
        log.warning("Assuming tex=trot")
        tex = trot
    elif isinstance(tex, dict):
        for k in tex:
            assert k in line_names,"{0} not in line list".format(k)
        line_names = tex.keys()
    elif tex > trot:
        warnings.warn("tex > trot in the ammonia model.  "
                      "This is unphysical and "
                      "suggests that you may need to constrain tex.  See "
                      "ammonia_model_restricted_tex.")
    if width < 0:
        return np.zeros(xarr.size)*np.nan
    elif width == 0:
        return np.zeros(xarr.size)

    
    from pyspeckit.spectrum.models.ammonia_constants import line_name_indices, line_names as original_line_names

    # recreate line_names keeping only lines with a specified tex
    # using this loop instead of tex.keys() preserves the order & data type
    line_names = [k for k in original_line_names if k in line_names]

    if 5 <= ntot <= 25:
        # allow ntot to be specified as a logarithm.  This is
        # safe because ntot < 1e10 gives a spectrum of all zeros, and the
        # plausible range of columns is not outside the specified range
        lin_ntot = 10**ntot
    else:
        print("ntot: {}".format(ntot))
        raise ValueError("ntot, the logarithmic total column density,"
                         " must be in the range 5 - 25")

    tau_dict = {}

    """
    Column density is the free parameter.  It is used in conjunction with
    the full partition function to compute the optical depth in each band
    """
    Zpara = (2*Jpara+1)*np.exp(-h*(Brot*Jpara*(Jpara+1)+
                                   (Crot-Brot)*Jpara**2)/(kb*trot))
    Zortho = 2*(2*Jortho+1)*np.exp(-h*(Brot*Jortho*(Jortho+1)+
                                       (Crot-Brot)*Jortho**2)/(kb*trot))
    Qpara = Zpara.sum()
    Qortho = Zortho.sum()

    log.debug("Partition Function: Q_ortho={0}, Q_para={1}".format(Qortho, Qpara))

    for linename in line_names:
        if ortho_dict[linename]:
            # define variable "ortho_or_para_frac" that will be the ortho
            # fraction in the case of an ortho transition or the para
            # fraction for a para transition
            ortho_or_parafrac = fortho
            Z = Zortho
            Qtot = Qortho
        else:
            ortho_or_parafrac = 1.0-fortho
            Z = Zpara
            Qtot = Qpara

        # for a complete discussion of these equations, please see
        # https://github.com/keflavich/pyspeckit/blob/ammonia_equations/examples/AmmoniaLevelPopulation.ipynb
        # https://github.com/pyspeckit/pyspeckit/blob/master/examples/AmmoniaLevelPopulation.ipynb
        # and
        # http://low-sky.github.io/ammoniacolumn/
        # and
        # https://github.com/pyspeckit/pyspeckit/pull/136

        # short variable names for readability
        frq = freq_dict[linename]
        partition = Z[line_name_indices[linename]]
        aval = aval_dict[linename]

        # Total population of the higher energy inversion transition
        population_rotstate = lin_ntot * ortho_or_parafrac * partition/Qtot

        if isinstance(tex, dict):
            expterm = ((1-np.exp(-h*frq/(kb*tex[linename]))) /
                       (1+np.exp(-h*frq/(kb*tex[linename]))))
        else:
            expterm = ((1-np.exp(-h*frq/(kb*tex))) /
                       (1+np.exp(-h*frq/(kb*tex))))
        fracterm = (ccms**2 * aval / (8*np.pi*frq**2))
        widthterm = (ckms/(width*frq*(2*np.pi)**0.5))

        tau_i = population_rotstate * fracterm * expterm * widthterm
        tau_dict[linename] = tau_i

        log.debug("Line {0}: tau={1}, expterm={2}, pop={3},"
                  " partition={4}"
                  .format(linename, tau_i, expterm, population_rotstate,
                          partition))

    # allow tau(11) to be specified instead of ntot
    # in the thin case, this is not needed: ntot plays no role
    # this process allows you to specify tau without using the approximate equations specified
    # above.  It should remove ntot from the calculations anyway...
    if tau is not None:
        tau11_temp = tau_dict['oneone']
        # re-scale all optical depths so that tau is as specified, but the relative taus
        # are sest by the kinetic temperature and partition functions
        for linename,t in iteritems(tau_dict):
            tau_dict[linename] = t * tau/tau11_temp

    if return_tau:
        return tau_dict
    
    model_spectrum = _ammonia_spectrum(xarr, tex=tex, tau_dict=tau_dict, width=width, xoff_v=xoff_v,
                                       fortho=fortho, line_names=line_names,
                                       background_ta=background_ta,
                                       fillingfraction=fillingfraction
                                      )
                                      #return_components=return_components
                                      #return_tau_profile=return_tau_profile

    '''
    if not return_tau_profile and model_spectrum.min() < 0 and background_ta == TCMB and not ignore_neg_models:
        raise ValueError("Model dropped below zero.  That is not possible "
                         " normally.  Here are the input values: "+
                         ("tex: {0} ".format(tex)) +
                         ("trot: %f " % trot) +
                         ("ntot: %f " % ntot) +
                         ("width: %f " % width) +
                         ("xoff_v: %f " % xoff_v) +
                         ("fortho: %f " % fortho)
                         )
    '''

    if verbose or debug:
        log.info("trot: %g  tex: %s  ntot: %g  width: %g  xoff_v: %g  "
                 "fortho: %g  fillingfraction: %g" % (trot, tex, ntot, width,
                                                      xoff_v, fortho,
                                                      fillingfraction))

    if bg_is_baseline:
    # return the model with the backgroud baseline subtracted
        return model_spectrum - background_ta
    else:
        return model_spectrum



def _ammonia_spectrum(xarr, tex, tau_dict, width, xoff_v, fortho, line_names, background_ta=0.0, fillingfraction=None,
                      return_components=False):
    """
    Helper function: given a dictionary of ammonia optical depths, an excitation tmeperature... etc, and produce a
    spectrum based on a one-slab (i.e., single velocity component model)
    Note: this is a modified version of the _ammonia_spectrum found in pyspeckit/spectrum/models/ammonia.py
    Note2: the final spectrum returned do not have the background emission subtracted
    Parameters
    ----------
    background_ta : float or ndarray
        "Antenna temperature" of the background emission
    Returns
    -------
    model : `model.SpectralModel`
        A SpectralModel class build from N different metastable inversion
        hyperfine models
    """

    # fillingfraction is an arbitrary scaling for the data; the model will be (normal model) * fillingfraction
    if fillingfraction is None:
        fillingfraction = 1.0

    # "runspec" means "running spectrum": it is accumulated over a loop
    runspec = np.zeros(len(xarr))

    # tau array that spans all the lines
    tau_arr = np.zeros(len(xarr))
    
    for linename in line_names:
        voff_lines = np.array(voff_lines_dict[linename])
        tau_wts = np.array(tau_wts_dict[linename])

        lines = (1-voff_lines/ckms)*freq_dict[linename]/1e9
        tau_wts = tau_wts / (tau_wts).sum()
        nuwidth = np.abs(width/ckms*lines)
        nuoff = xoff_v/ckms*lines

        # tau array for each (J,K) lines
        tauprof = np.zeros(len(xarr))
        for kk,nuo in enumerate(nuoff):
            tauprof_ = (tau_dict[linename] * tau_wts[kk] *
                        np.exp(-(xarr.value+nuo-lines[kk])**2 /
                               (2.0*nuwidth[kk]**2)))
            tauprof += tauprof_

        T0 = (h*xarr.value*1e9/kb) # "temperature" of wavelength

        if isinstance(tex, dict):
            runspec = ((T0/(np.exp(T0/tex[linename])-1)*(1-np.exp(-tauprof)))*fillingfraction
                       + runspec)

        else:
            runspec = ((T0/(np.exp(T0/tex)-1)*(1-np.exp(-tauprof)))*fillingfraction
                       + runspec)
            
        tau_arr = tau_arr + tauprof

    # note, the background emission is assumed to be beam-filling
    return runspec + background_ta*np.exp(-tau_arr)


def T_antenna(Tbright, nu):
    """
    Calculate antenna temperatures over nu (in GHz)
    """
    T0 = (h*nu*1e9/kb)
    return T0/(np.exp(T0/Tbright)-1)



#=======================================================================================================================
# for generating registry with pyspeckit fitters

def nh3_multi_v_model_generator(n_comp, linenames = None):
    """
    Works for up to 2 componet fits at the moment
    Parameters
    ----------
    n_comp : int
        The number of velocity componets to fit
    linenames : list
        A list of line names from the set ('oneone', ..., 'eighteight'); default is just 'oneone'
    Returns
    -------
    model : `model.SpectralModel`
        A SpectralModel class build from N different metastable inversion
        hyperfine models
    """
    n_para_slab = 6
    n_para = n_comp*n_para_slab
    idx_comp = np.arange(n_comp)

    if linenames is None:
        #linenames = ['oneone']
        linenames = line_names

    def nh3_vtau_multimodel(xarr, *args):
        assert len(args) == n_para
        return ammonia_multi_v(xarr, n_comp, *args, line_names=linenames)


    # define parameter limits
    limitedmin = (True, True, True, True, False, True)
    limitedmax = (False, False, True, False, False, True)
    minpars = (TCMB, TCMB, 5, 0, 0, 0)
    maxpars = (0, 0, 25, 0, 0, 1)

    if n_comp >= 2:
        for i in range(2, n_comp+1):
            limitedmin += limitedmin
            limitedmax += limitedmax
            minpars += minpars
            maxpars += maxpars

    parlimited = list(zip(limitedmin, limitedmax))
    parlimits = list(zip(minpars, maxpars))

    #print(parlimited)
    #print(parlimits)

    mod = model.SpectralModel(nh3_vtau_multimodel, n_para,
            parnames=[x
                      for ln in idx_comp
                      for x in ('trot{0}'.format(ln),
                                'tex{0}'.format(ln),
                                'ntot{0}'.format(ln),
                                'width{0}'.format(ln),
                                'xoff_v{0}'.format(ln),
                                'fortho{0}'.format(ln)
                                )],
            #parlimited=[(False,False), (True,False), (True,False), (True,False)]*n_para,
            #parlimits=[(0,0), ]*n_para,
            parlimited=parlimited,
            parlimits=parlimits,

            shortvarnames=[x
                           for ln in idx_comp
                           for x in ('T_{{rot,{0}}}'.format(ln),
                                     'T_{{ex,{0}}}'.format(ln),
                                     'N_{{tot,{0}}}'.format(ln),
                                     '\\sigma_{{v,{0}}}'.format(ln),
                                     'v_{{lsr,{0}}}'.format(ln),
                                     'f_{{orth,{0}}}'.format(ln)
                                     )],
            fitunit='GHz')
            # the keyword fitunits is now fitunit?
    return mod
