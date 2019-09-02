from setuptools import setup

setup(
    name='vsomm_modeler',
    version='2.0.0',
    packages=['vsomm_modeler'],
    scripts=['bin/vsomm-modeler'],
    url='http://somm.boku.ac.at/',
    license='GPLv3.0',
    author='Yerko Escalona',
    author_email='yerko.escalona@boku.ac.at',
    description='Simulation and Modelling Art',
    include_package_data=True,
)
