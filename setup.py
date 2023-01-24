from setuptools import setup

setup(name='Gaia_spec',
      version='0.0.1',
      description='Compare Gaia spectra to templates',
      url='https://github.com/fkiwy/Gaia_spec',
      author='Frank Kiwy',
      author_email='frank.kiwy@outlook.com',
      license='MIT',
      packages=['Gaia_spec'],
      install_requires=['numpy', 'astropy', 'matplotlib', 'specutils', 'GaiaXPy'],
      zip_safe=False,
      include_package_data=True)
