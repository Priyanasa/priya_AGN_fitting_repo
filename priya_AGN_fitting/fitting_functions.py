from astropy.wcs import WCS
import numpy as np
import matplotlib.pyplot as plt
from lmfit.models import PowerLawModel, GaussianModel, VoigtModel, SkewedGaussianModel, SkewedVoigtModel

emiss_line_dict = {'Ly_alpha':1216, 'NV':1240, 'SiIV':1397, 'CIV':1549, 'HeII':1640}

def get_coords_from_header(header):
    """Takes a spectral cube FITS header and returns RA, Dec, and
    Wavelength arrays
    
    Parameters
    ==========
    header : class 'astropy.io.fits.header.Header'
        The header from a spectral FITS cube as read in via `astropy.io.fits`
    
    Returns
    =======
    ras : numpy array
        The range of RAs corresponding to NAXIS1 of the FITS file (usually deg)
    decs : numpy array
        The range of Decs corresponding to NAXIS2 of the FITS file (usually deg)
    wavelengths : numpy array
        The range of wavelengths corresponding to NAXIS3 of the FITS file (usually :math:`\\mathring{\\mathrm{A}}`)
    """

    ##Grab the world coord system
    wcs = WCS(header)

    ##Make an array containing all the pixel values in the x axis, using
    ##the NAXIS1 header value
    ra_pix = np.arange(int(header['NAXIS1']))
    ##Only use y z = 0,0 to return just dec
    ##Anytime we don't care about an output, we set it equal to '_'
    ras, _, _ = wcs.all_pix2world(ra_pix,0,0,0)

    dec_pix = np.arange(int(header['NAXIS2']))
    ##Only use x, z = 0,0 to return just dec
    _, decs, _ = wcs.all_pix2world(0,dec_pix,0,0)

    wave_pix = np.arange(int(header['NAXIS3']))
    ##Only use x, y = 0,0 to return just wavelengths
    _, _, wavelengths = wcs.all_pix2world(0,0,wave_pix,0)

    return ras, decs, wavelengths

def plot_spectra(wavelengths, spectra, smoothed=False):
    """Plot a spectra as a function of wavelength on a 1D plot
    (`wavelengths` x-axis, `spectra` y-axis).
    If `smoothed` is given, plot the smoothed spectra result over
    the plot also. Plot is saved to 'input_spectra.png'.
    
    Parameters
    ----------
    wavelengths : numpy array
        Array containing wavelengths to plot (:math:`\\mathring{\\mathrm{A}}`)
    spectra : numpy array
        Array containing flux densities :math:`(\\mathrm{erg}\\,\\mathrm{s}^{-1}\\,\\mathrm{cm}^{-2}\\,\\mathring{\\mathrm{A}}^{-1})`
    smoothed : numpy array
        Optional array containing smoothed flux densities :math:`(\\mathrm{erg}\\,\\mathrm{s}^{-1}\\,\\mathrm{cm}^{-2}\\,\\mathring{\
        \mathrm{A}}^{-1})`
    
    Returns
    -------
    """

    ##Plotting the spectra
    fig, ax = plt.subplots(1,1,figsize=(12,5))

    ax.plot(wavelengths, spectra, 'gray', label='FITS data', alpha=0.7)

    if type(smoothed) == np.ndarray:
        ax.plot(wavelengths, smoothed, 'k', label='Smoothed')

    for colour,key in enumerate(emiss_line_dict.keys()):
        ax.axvline(emiss_line_dict[key], color="C{:d}".format(colour), linestyle='--', label=key)

    ax.legend()
    ax.set_ylabel('$F_{\lambda}$ ($\mathrm{erg}\,\mathrm{s}^{-1} \mathrm{cm}^{-2}$ Å$^{-1}$)')
    ax.set_xlabel('Wavelength (Å)')

    fig.savefig('input_spectra.png',bbox_inches='tight')
    plt.close()
    
    return
    
def spectra_subset(wavelengths, spectra, width=250, line_cent=False,
                   lower_wave=False, upper_wave=False):
    """Takes the given spectra and wavelengths (in
    angstrom) and crops. If line_cent is provided, crop about line_cent
    to the given width (defaults to 250 Angstrom). Alternatively,
    can manually specifiy a lower wavelength bounary (lower_wave)
    and upper wavelength bounary (upper_wave). lower_wave and upper_wave
    will overwrite any limits set via line_cent and width.
    
    Parameters
    ----------
    wavelengths : numpy array
        Array containing wavelengths to plot (:math:`\\mathring{\\mathrm{A}}`)
    spectra : numpy array
        Array containing flux densities :math:`(\\mathrm{erg}\\,\\mathrm{s}^{-1}\\,\\mathrm{cm}^{-2}\\,\\mathring{\\mathrm{A}}^{-1})`
    width : float
        Width about `line_cent` to crop the spectra to (:math:`\\mathring{\\mathrm{A}}`)
    line_cent : float
        Central wavelength (:math:`\\mathring{\\mathrm{A}}`) to crop about
    lower_wave: float
        Alternatively, specify a lower wavelength bound to crop to (:math:`\\mathring{\\mathrm{A}}`)
    upper_wave: float
        Alternatively, specify an upper wavelength bound to crop to (:math:`\\mathring{\\mathrm{A}}`)
    
    Returns
    -------
    spectra[indexes] : numpy array
        A cropped subset of `spectra`
    wavelengths[indexes] : numpy array
        A cropped subset of `wavelengths`
    """

    ##Exit if reasonable combo of arguements not given
    if not line_cent and not lower_wave and not upper_wave:
        sys.exit("ERROR in spectra_subset - if line_cent is not set, both lower_wave and upper_wave must be set. Exiting.")

    if line_cent:
        lower_crop = line_cent - width / 2
        upper_crop = line_cent + width / 2

    if lower_wave:
        lower_crop = lower_wave

    if upper_wave:
        upper_crop = upper_wave

    ##The where function returns an array of indexes where the boolean logic is true
    indexes = np.where((wavelengths >= lower_crop) & (wavelengths <= upper_crop))
    ##Use our result from where to crop the spectra and wavelengths
    return spectra[indexes], wavelengths[indexes]

def do_lmfit(lm_model, trim_wavelengths, trim_spectra, line_cent):
    """Takes an lmfit.models instance `lm_model` and fits an emission line
    centred at `line_cent`, as described by `trim_spectra`, `trim_wavelenghts`
    
    Parameters
    ----------
    lm_model : `lmfit.model`
        The lmfit model to be used in fitting
    trim_wavelengths : numpy array
        Array containing wavelengths to fit (:math:`\\mathring{\\mathrm{A}}`)
    trim_spectra : numpy array
        Array containing flux densities to fit :math:`(\\mathrm{erg}\\,\\mathrm{s}^{-1}\\,
        \\mathrm{cm}^{-2}\\,\\mathring{\\mathrm{A}^{-1})`
    line_cent : float
        Central wavelength (:math:`\\mathring{\\mathrm{A}}`) to fit about
    lower_wave: float
        Alternatively, specify a lower wavelength bound to crop to (:math:`\\mathring{\\
        mathrm{A}}`)
    upper_wave: float
        Alternatively, specify an upper wavelength bound to crop to (:math:`\\mathring{\\
        mathrm{A}}`)
    
    Returns
    -------
    fit : `lmfit.model.ModelResult`
        The results of the `lmfit` fitting
    """
    
    ##You can combine lmfit models just by adding them like this:
    emiss_and_power_model = lm_model(prefix='emission_') + PowerLawModel(prefix='power_')

    ##You create the params in the same way as you would a single model
    emiss_and_power_params = emiss_and_power_model.make_params()

    ##Start the fitting at the line centre that we want
    emiss_and_power_params['emission_center'].set(value=line_cent)

    ##Do the fit
    fit = emiss_and_power_model.fit(trim_spectra, emiss_and_power_params, x=trim_wavelengths)

    ##Print the report
    print(fit.fit_report())
    
    return fit

def do_fit_plot(rest_wavelengths, spectra, trim_wavelengths, trim_spectra, fit):
    """Plots the given input spectra (`rest_wavelengths`, `spectra`) the subset
    of data that was used for fitting (`trim_wavelenghts, `trim_spectra`), and
    the fit result out of lmfit using `fit`. Saves the figure as
    'fit_results.png'.
    
    Parameters
    ==========
    rest_wavelengths : numpy array
        Array containing all avaible wavelengths (:math:`\\mathring{\\mathrm{A}}`)
    spectra : numpy array
        Array containing all avaiable flux densities :math:`(\\mathrm{erg}\\,\\mathrm{s}^{-1}\\,\\mathrm{cm}^{-2}\\,\\mathring{\\mathrm{A}}
        ^{-1})`
    trim_wavelengths : numpy array
        Array containing wavelengths used for fit (:math:`\\mathring{\\mathrm{A}}`)
    trim_spectra : numpy array
        Array containing flux densities used for fit :math:`(\\mathrm{erg}\\,\\mathrm{s}^{-1}\\,\\mathrm{cm}^{-2}\\,\\mathring{\\mathrm{A}}
        ^{-1})`
    fit : `lmfit.model.ModelResult`
        The results of the `lmfit` fitting
    
    Returns
    =======
    """

    ##Plotting the fit
    fig, axs = plt.subplots(1,2,figsize=(12,4))

    axs[0].plot(rest_wavelengths, spectra, 'gray', label='Input spectra',lw=2.0, alpha=0.5)

    for ax in axs:
        ax.plot(trim_wavelengths, trim_spectra, 'k', label='Fitted spectra',lw=2.0)
        ax.plot(trim_wavelengths,fit.best_fit,lw=2.0, label='Fit result')

        ax.set_xlabel('Rest wavelength (Å)')
        ax.legend()

    axs[0].set_ylabel('$F_{\lambda}$ ($\mathrm{erg}\,\mathrm{s}^{-1} \mathrm{cm}^{-2}$ Å$^{-1}$)')

    fig.savefig('fit_results.png',bbox_inches='tight')
    plt.close()
    
    return
