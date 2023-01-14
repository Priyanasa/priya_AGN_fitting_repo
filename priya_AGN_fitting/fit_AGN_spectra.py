from astropy.io import fits
from scipy.signal import savgol_filter
from lmfit.models import PowerLawModel, GaussianModel, VoigtModel, SkewedGaussianModel, SkewedVoigtModel
import argparse
import fitting_functions as ff

def get_parser():
    """Wrapper function for the argparse commands"""
    
    ##This sets up a parser, which can read inputs given at the command line
    parser = argparse.ArgumentParser(description='A script to read in a spectral \
             AGN data and fit a given emission line at a given x,y pixel \
             location ')

    ##Adding arguments
    parser.add_argument('fitsfile', type=str,
                        help='3D FITS file with AGN data')
    parser.add_argument('x_pix', type=int,
                        help='The x coord of the data to fit (pixel coords)')
    parser.add_argument('y_pix', type=int,
                        help='The x coord of the data to fit (pixel coords)')
    parser.add_argument('--redshift', type=float, default=2.73763,
                        help='Redshift of the AGN - defaults to 2.73763')
    parser.add_argument('--smooth_data', action='store_true',
                        help='If added, apply a SavGol filter to the spectra')
    parser.add_argument('--savgol_wl', type=int, default=21,
                        help='Window length in pixels to use in SavGol smooting. \
                        Note this parameter must be an odd integer. Default=21')
    parser.add_argument('--savgol_poly', type=int, default=1,
                        help='Order of polynomial used in SavGol smooting.')
    parser.add_argument('--emission_line', default='CIV', type=str,
                        help="Try fitting this known line. Options are: (Angstroms).: \
                        'Ly_alpha' (1216), \
                        'NV' (1240), \
                        'SiIV' (1397), \
                        'CIV' (1549), \
                        'HeII' (1640). \
                        Defaults to 'CIV'")
    parser.add_argument('--fitting_model', default='SkewedVoigtModel', type=str,
                        help="Which lmfit model to use when fitting emission. Options \
                        are: GaussianModel, VoigtModel, SkewedGaussianModel, \
                        SkewedVoigtModel. Defaults to SkewedVoigtModel")
    parser.add_argument('--fit_width', default=250, type=int,
                        help="Width of data to fit in Angstrom")
    parser.add_argument('--lower_wave', default=False, type=int,
                        help="Instead of fitting a specific width of spectra with \
                        --fit_width, specify a lower boundary (in Angstrom) to fit to")
    parser.add_argument('--upper_wave', default=False, type=int,
                        help="Instead of fitting a specific width of spectra with \
                        --fit_width, specify an upper boundary (in Angstrom) to fit to")
    parser.add_argument('--no_spectra_plot', action='store_true',
                        help="If added, do not create input spectra plot")
    parser.add_argument('--no_fit_plot', action='store_true',
                        help="If added, do not create fiting plot, just report fit results")
    
    return parser

if __name__ == '__main__':
    import sys
   
    #Getting argparse
    parser = get_parser()
    
    ##This grabs all of the arguments out of the parser. Now all the arguments
    ##are attributes of args, which we can access and use
    args = parser.parse_args()

    ##This is only going to work on FITS data that has the same structure as
    ##1009_629_2006A_individual_cubes_3D.fits
    with fits.open(args.fitsfile) as hdu:
    ##Get the data from the 2nd hdu entry
        data_cube = hdu[1].data
        fits_header = hdu[1].header

    ##Check whether our input arguments make sense, and exit with error message
    ##if not
    if args.x_pix < 0 or args.x_pix >= fits_header['NAXIS1']:
        sys.exit("ERROR: --x_pix={:d} is outside the range of the FITS NAXIS1={:d}".format(args.x_pix,fits_header['NAXIS1']))

    if args.y_pix < 0 or args.y_pix >= fits_header['NAXIS2']:
        sys.exit("ERROR: --y_pix={:d} is outside the range of the FITS NAXIS2={:d}".format(args.y_pix,fits_header['NAXIS2']))

    if args.fitting_model not in ["GaussianModel", "VoigtModel", "SkewedGaussianModel", "SkewedVoigtModel"]:
        print("Model specified --fitting_model={:s} is not supported. Defaulting to GaussianModel".format(args.fitting_model))
        args.fitting_model = "GaussianModel"

    if args.emission_line not in ['Ly_alpha', 'NV', 'SiIV', 'CIV', 'HeII']:
        print("Line specified in --emission_line={:s} is not supported. Defaulting to CIV".format(args.emission_line))
        args.emission_line = "CIV"

    ##Grab some relevant coords from the header
    ras, decs, wavelengths = ff.get_coords_from_header(fits_header)

    ##Remove redshifting on the spectra
    rest_wavelengths = wavelengths / (args.redshift + 1)

    ##Cut to the specific pixel
    spectra = data_cube[:,args.y_pix,args.x_pix]

    ##Do the smoothing if asked for it
    if args.smooth_data:
        smoothed_spectra = savgol_filter(spectra, args.savgol_wl, args.savgol_poly)
    else:
        ##If no smoothing, set this to False for use in future functions
        smoothed_spectra = False

    ##Plot input data if requested
    if args.no_spectra_plot:
        pass
    else:
        ff.plot_spectra(rest_wavelengths, spectra, smoothed=smoothed_spectra)

    ##Set lm_model to the requested function form for the emission peak
    if args.fitting_model == "GaussianModel":
        lm_model = GaussianModel
    elif args.fitting_model == "VoigtModel":
        lm_model = VoigtModel
    elif args.fitting_model == "SkewedGaussianModel":
        lm_model = SkewedGaussianModel
    elif args.fitting_model == "SkewedVoigtModel":
        lm_model = SkewedVoigtModel

    ##Set the line_cent using the emiss_line_dict
    line_cent = ff.emiss_line_dict[args.emission_line]

    ##Specify which data to fit, based on whether we are smoothing or not
    if args.smooth_data:
        spectra_to_fit = smoothed_spectra
    else:
        spectra_to_fit = spectra

    ##Trim the spectra as required
    trim_spectra, trim_wavelengths = ff.spectra_subset(rest_wavelengths, spectra_to_fit,
                                     width=args.fit_width, line_cent=line_cent,
                                     lower_wave=args.lower_wave,
                                     upper_wave=args.upper_wave)

    ##Print out a summary of fitting settings
    print('Fitting using the following inputs:')
    print('\tfitting_model: {:s}'.format(args.fitting_model))
    print('\temission_line: {:s} ({:.1f} Angstrom)'.format(args.emission_line, line_cent))
    print('\twavelength bounds: {:.1f} to {:.1f} (Angstrom)'.format(trim_wavelengths[0], trim_wavelengths[-1]))

    ##Do the fit
    fit = ff.do_lmfit(lm_model, trim_wavelengths, trim_spectra, line_cent)

    ##Plot fit if required
    if args.no_fit_plot:
        pass
    else:
        ff.do_fit_plot(rest_wavelengths, spectra, trim_wavelengths, trim_spectra, fit)
