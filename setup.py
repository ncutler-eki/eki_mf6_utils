from setuptools import setup, find_packages

setup(
    name='eki_mf6_utils',
    version='0.0.1',
    install_requires=['numpy',
                      'matplotlib',
                      'fiona',
                      'rasterio',
                      'flopy @ git+https://github.com/mmaneta/flopy.git@nonzero_lgr',
                      'jupyterlab'],
    packages=find_packages(),
    url='',
    license='',
    author='Marco Maneta',
    author_email='mmaneta@ekiconsult.com',
    description='Utilities to manipulate MODFLOW 6 models'
)
