from setuptools import setup

setup(name='Gaia_spec',
      version='1.0.0',
      description='Compare Gaia spectra to templates',
      url='https://github.com/fkiwy/Gaia_spec',
      author='Frank Kiwy',
      author_email='frank.kiwy@outlook.com',
      license='MIT',
      py_modules=['Spec_compare'],
      install_requires=['numpy', 'astropy', 'matplotlib', 'specutils', 'GaiaXPy'],
      include_package_data=True)
