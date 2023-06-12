import numpy as np
from pyspeckit.spectrum.models import model
from pyspeckit.spectrum.models.ammonia_constants import line_names


from models import ammonia as nh3

'''
========================================
Multi-component Ammonia inversion transition TROT fitter
========================================
'''

TCMB = 2.7315 # K

def nh3_fixfortho(xarr, ncomp, *parameters, background_tb=TCMB, fortho=0.0, **kwargs):
    """
    A wrapper
    :param xarr:
    :param ncomp:
    :param parameters:
        a list of trot, tex, ntot, width, xoff_v, fortho for each component
    :param background_tb:
    :param fortho:
        <float> the fortho ration to be fixed at
    :param kwargs:
    :return:
    """

    # number of parameters per slab
    n_para_slab = 5
    nprs = n_para_slab

    # Convert X-units to frequency in GHz
    if xarr.unit.to_string() != 'GHz':
        xarr = xarr.as_unit('GHz')

    args = parameters
    # need to sanity check the number of parameters is consistent with the number of components
    if len(args)/ncomp != n_para_slab:
        print("[ERROR]: the number of parameters ({}) does not match the number of components ({})".format(len(args)                                                                                      ,ncomp))
        return None

    # Start with a constant background radiation field (e.g., Cosmic Microwave Background by default)
    BGRF = nh3.T_antenna(background_tb, xarr.value)
    model = BGRF

    #args = parameters
    # iteratively move through each slabs towards the observer (i.e., radiative transfer)
    for trot, tex, ntot, width, xoff_v in zip(args[::nprs], args[1::nprs], args[2::nprs],
                                                      args[3::nprs], args[4::nprs]):
        model = nh3.ammonia_slab(xarr, trot=trot, tex=tex, ntot=ntot, width=width, xoff_v=xoff_v, fortho=fortho,
                                     background_ta=model, **kwargs)

    return model - BGRF


def nh3_fixfortho_model(n_comp, linenames = None):
    """
    Generates a pyspeckit SpectralModel for multi-component ammonia model (i.e., fitter for the fitting)

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
    n_para_slab = 5
    n_para = n_comp*n_para_slab
    idx_comp = np.arange(n_comp)

    if linenames is None:
        linenames = line_names

    def nh3_vtau_multimodel(xarr, *args):
        assert len(args) == n_para
        return nh3_fixfortho(xarr, n_comp, *args, line_names=linenames)


    # define parameter limits
    limitedmin = (True, True, True, True, False)
    limitedmax = (False, False, True, False, False)
    minpars = (TCMB, TCMB, 5, 0, 0)
    maxpars = (0, 0, 25, 0, 0)

    if n_comp >= 2:
        # repeats the limit for each slab
        for i in range(2, n_comp+1):
            limitedmin += limitedmin
            limitedmax += limitedmax
            minpars += minpars
            maxpars += maxpars

    parlimited = list(zip(limitedmin, limitedmax))
    parlimits = list(zip(minpars, maxpars))


    mod = model.SpectralModel(nh3_vtau_multimodel, n_para,
            parnames=[x
                      for ln in idx_comp
                      for x in ('trot{0}'.format(ln),
                                'tex{0}'.format(ln),
                                'ntot{0}'.format(ln),
                                'width{0}'.format(ln),
                                'xoff_v{0}'.format(ln)
                                )],
            parlimited=parlimited,
            parlimits=parlimits,

            shortvarnames=[x
                           for ln in idx_comp
                           for x in ('T_{{rot,{0}}}'.format(ln),
                                     'T_{{ex,{0}}}'.format(ln),
                                     'N_{{tot,{0}}}'.format(ln),
                                     '\\sigma_{{v,{0}}}'.format(ln),
                                     'v_{{lsr,{0}}}'.format(ln)
                                     )],
            fitunit='GHz')
    return mod
