import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table
# from astropy.stats import sigma_clip
from astroquery.simbad import Simbad
from astroquery.gaia import Gaia
from gaiaxpy import calibrate


spts = ['B', 'A', 'F', 'G', 'K', 'M', 'L']

subtypes = []

for spt in spts:
    for i in range(10):
        subtypes.append(f'{spt}{i}V')

for spt in subtypes:
    print(f'\nCreating {spt} template ...')

    simbad = Simbad()
    simbad.TIMEOUT = 300
    simbad.ROW_LIMIT = 5000
    simbad.add_votable_fields('ids')
    results = simbad.query_criteria(f"sptype='{spt}' & spqual in ('A', 'B', 'C')")

    if results is None:
        print('  Number of sources downloaded from SIMBAD: 0')
        continue

    print('  Number of sources downloaded from SIMBAD:', len(results))

    sources = []
    i = 0
    for row in results:
        ids = row['IDS']

        gaia_id = None
        ids = ids.split('|')
        for id_ in ids:
            if 'Gaia DR3' in id_:
                gaia_id = id_[9:]
                break

        if gaia_id:
            sources.append(gaia_id)
            i += 1

    print('  Number of sources with a Gaia identifier:', len(sources))

    if len(sources) == 0:
        continue

    ids = ','.join(sources)
    job = Gaia.launch_job(f'select top 1000 source_id from gaiadr3.gaia_source where source_id in ({ids}) and has_xp_continuous = \'true\'')
    sources_with_spectrum = job.get_results()['source_id'].data.tolist()

    number_of_spectra = len(sources_with_spectrum)
    print('  Number of sources with a Gaia spectrum  :', number_of_spectra)

    if number_of_spectra == 0:
        continue

    output_data, output_sampling = calibrate(sources_with_spectrum, output_path='data', save_file=False)

    fluxes = []

    for index, row in output_data.iterrows():
        fluxes.append(row['flux'])

    # fluxes = np.array(fluxes)
    # fluxes = sigma_clip(fluxes, axis=0, sigma=3, maxiters=None)
    # flux = np.median(fluxes, axis=0)

    flux = np.mean(fluxes, axis=0)
    flux_error = np.std(fluxes, axis=0) / np.sqrt(np.size(fluxes, axis=0))
    wavelength = output_sampling

    flux_unit = 'W/(nm*m^2)'
    template = Table([wavelength, flux, flux_error], names=['Wavelength', 'Flux', 'Flux_error'], units=['nm', flux_unit, flux_unit])
    template.write(f'templates/Gaia_template_{spt}.dat', format='ipac', overwrite=True)

    fontsize = 18
    plt.rcParams.update({'font.size': fontsize})
    fig, ax = plt.subplots(figsize=(11, 8))
    ax.plot(wavelength, flux, label=f'Template {spt}')
    ax.plot(wavelength, flux + flux_error, lw=0.5, c='gray', zorder=1, label='Standard error')
    ax.plot(wavelength, flux - flux_error, lw=0.5, c='gray', zorder=1)
    ax.set_xlabel('Wavelength [nm]')
    ax.set_ylabel('Flux [W$\\cdot$nm$^{-1}$$\\cdot$m$^{-2}$]')
    ax.grid(color='grey', alpha=0.8, linestyle='-.', linewidth=0.2, axis='both')
    plt.legend(loc='best')
    plt.title(f'Number of Gaia spectra for {spt}: {number_of_spectra}', fontdict={'fontsize': fontsize})
    plt.savefig(f'templates/figures/Gaia_template_{spt}.pdf', dpi=600, bbox_inches='tight')
    plt.close()
