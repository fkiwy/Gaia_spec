from Spec_compare import compare_spectra

def test_compare_spectra():
    sources = [58925512487888896, 378026395579692672, 200296663143599104, 1204754033223498624]  # K5, M0, M5, M9

    results = compare_spectra(sources, plot_spectrum=False, plot_template=True, spectrum_uncertainty=True, template_uncertainty=True, file_format='png')

    assert results[0][2] == 'K5V'
    assert results[1][2] == 'M0V'
    assert results[2][2] == 'M5V'
    assert results[3][2] == 'M9V'
