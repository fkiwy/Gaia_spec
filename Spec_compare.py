import os
from os import listdir
from os.path import isfile, join
import tempfile
import warnings
import numpy as np
from astropy.io import ascii
from astropy import units as u
from astropy.nddata import StdDevUncertainty
from specutils import Spectrum1D
from specutils.analysis import template_comparison
import matplotlib.pyplot as plt
from gaiaxpy import calibrate


def compare_spectra(sources, plot_spectrum=True, plot_template=True, spectrum_uncertainty=False, template_uncertainty=False,
                    output_dir=tempfile.gettempdir(), file_format='pdf'):
    """
    Compares Gaia spectra to templates

    Parameters
    ----------
    sources : list
        A list of Gaia source identifiers to download spectra for.
        At least one is required e.g. [1657463068195202432]
    plot_spectrum : bool, optional
        Whether to plot the downloaded Gaia spectra. The default is True.
    plot_template : bool, optional
        Whether to overplot the matched templates. The default is True.
    spectrum_uncertainty : bool, optional
        Whether to plot the spectra uncertainty. The default is False.
    template_uncertainty : bool, optional
        Whether to plot the templates uncertainty. The default is False.
    output_dir : str, optional
        The directory where the plots will be saved to. The default is tempfile.gettempdir().
    file_format : str, optional
        Output file format: pdf, png, eps, etc.. The default is 'pdf'.

    Returns
    -------
    A list containing one row per source including the Gaia source id, the Gaia spectrum, the matched spectral type and the matched template.
    All spectra and tempaltes are instances of specutils.Spectrum1D.
    Wavelength, flux and error arrays can be extracted as follows:
    - Spectrum1D.spectral_axis.value
    - Spectrum1D.flux.value
    - Spectrum1D.uncertainty.array

    """

    def normalize(data, mini, maxi):
        dmin = np.nanmin(data)
        dmax = np.nanmax(data)
        return (maxi - mini) * (data - dmin) / (dmax - dmin) + mini

    def create_spectrum_figure(spectrum, spectrum_label, template, template_label, template_error):
        fontsize = 18
        plt.rcParams.update({'font.size': fontsize})
        fig, ax = plt.subplots(figsize=(11, 8))

        ax.plot(spectrum.spectral_axis, spectrum.flux, c='navy', zorder=0, label=spectrum_label)
        if spectrum_uncertainty:
            ax.plot(spectrum.spectral_axis, spectrum.flux.value + spectrum.uncertainty.array, lw=0.5, c='blue', zorder=1)
            ax.plot(spectrum.spectral_axis, spectrum.flux.value - spectrum.uncertainty.array, lw=0.5, c='blue', zorder=1)

        if plot_template:
            ax.plot(template.spectral_axis, template.flux, c='red', label=f'Template {template_label}', zorder=2)
            if template_uncertainty:
                ax.plot(template.spectral_axis, template.flux.value + template_error, lw=0.5, c='pink', zorder=3)
                ax.plot(template.spectral_axis, template.flux.value - template_error, lw=0.5, c='pink', zorder=3)

        ax.set_xlabel('Wavelength [nm]')
        ax.set_ylabel('Flux [W$\\cdot$nm$^{-1}$$\\cdot$m$^{-2}$]')
        ax.grid(color='grey', alpha=0.5, linestyle='-.', linewidth=0.2, axis='both')
        plt.legend(loc='best')
        plt.title('Gaia low-resolution BP-RP spectrum', fontdict={'fontsize': fontsize})
        file_name = f'{template_label}_{spectrum_label}.{file_format}'
        plt.savefig(os.path.join(output_dir, file_name), dpi=600, bbox_inches='tight')
        plt.close()

    def compare_spectrum(spectrum, spectrum_label, templates):
        template_spectra = [t[0] for t in templates]
        template_labels = [t[1] for t in templates]

        template, _, smallest_chi_index, _, _ = template_comparison.template_match(
            observed_spectrum=spectrum, spectral_templates=template_spectra)

        template_error = template_spectra[smallest_chi_index].uncertainty.array
        spectrum_error = spectrum.uncertainty.array
        template_error = normalize(template_error, np.nanmin(spectrum_error), np.nanmax(spectrum_error))

        template_label = template_labels[smallest_chi_index]
        print('Matched template:', template_label)

        if plot_spectrum:
            create_spectrum_figure(spectrum, spectrum_label, template, template_label, template_error)

        return template_label, template

    warnings.simplefilter('ignore', category=Warning)

    # Spectral units
    wavelength_unit = u.nm
    flux_unit = u.W / u.nm / u.m**2

    # Load template spectra
    templates = []

    template_dir = 'templates/'
    template_paths = [template_dir + f for f in listdir(template_dir) if isfile(join(template_dir, f))]

    for template_path in template_paths:
        template = ascii.read(template_path, format='ipac')

        wavelength = template['Wavelength']
        flux = template['Flux']
        flux_error = template['Flux_error']

        wavelength = wavelength.data * wavelength_unit
        flux = flux * flux_unit
        flux_error = flux_error * flux_unit

        template_spectrum = Spectrum1D(spectral_axis=wavelength, flux=flux, uncertainty=StdDevUncertainty(flux_error))
        template_label = template_path[template_path.rfind('_')+1:template_path.rfind('.')]
        templates.append((template_spectrum, template_label))

    # Download Gaia spectra
    spectra = []

    output_data, output_sampling = calibrate(sources, save_file=False)
    wavelength = output_sampling * wavelength_unit

    for _, row in output_data.iterrows():
        source_id = row['source_id']
        flux = row['flux']
        flux_error = row['flux_error']

        flux = flux * flux_unit
        flux_error = flux_error * flux_unit

        spectrum = Spectrum1D(spectral_axis=wavelength, flux=flux, uncertainty=StdDevUncertainty(flux_error))
        spectra.append((spectrum, str(source_id)))

    # Compare Gaia to template spectra
    total_spectra = len(spectra)
    spectrum_number = 1

    results = []

    print('\nCompare Gaia to template spectra:')
    for spectrum in spectra:
        spectrum_data = spectrum[0]
        spectrum_label = spectrum[1]
        print('\n' + str(spectrum_number) + '/' + str(total_spectra) + ' ' + spectrum_label)
        template_label, template = compare_spectrum(spectrum_data, spectrum_label, templates)
        results.append([spectrum_label, spectrum, template_label, template])
        spectrum_number += 1

    return results
