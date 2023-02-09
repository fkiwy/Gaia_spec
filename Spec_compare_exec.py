from Spec_compare import compare_spectra

sources = [58925512487888896, 378026395579692672, 200296663143599104, 1204754033223498624]  # K5, M0, M5, M9

results = compare_spectra(sources, plot_spectrum=True, plot_template=True, spectrum_uncertainty=True, template_uncertainty=True, file_format='png', output_dir='output')

for result in results:
    print('Gaia source id:', result[0], 'Matched tempalte:', result[2])

"""
Console output:
Gaia source id: 58925512487888896 Matched tempalte: K5V
Gaia source id: 378026395579692672 Matched tempalte: M0V
Gaia source id: 200296663143599104 Matched tempalte: M5V
Gaia source id: 1204754033223498624 Matched tempalte: M9V
"""
