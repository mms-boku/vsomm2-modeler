VSOMM2-modeler
==============

The Vienna Soil Organic Matter **Modeler** is a PYTHON program able to create condensed phase models of humic substances.

Dependencies
------------

Download vsomm2-building_blocks:
```bash
$ git clone https://gromos.boku.ac.at/yescalona/vsomm2-building_blocks.git
```

Install python dependencies:
```bash
$ pip install -r requirements
```



Optionally, install SMArt for GROMACS compatibility
```bash
$ git clone https://bitbucket.org/durozlikovski/smart.git
$ cd smart
$ pip install .
```

Install
-------

Set environmental variables in your *.bashrc*.
```bash
export GROMOSXX_BIN="path_to_gromosXX/bin/"
export GROMOSPP_BIN="path_to_gromos++/bin/"
export VSOMM_BUILDING_BLOCKS="path_to_vsomm2-building_blocks/"
```

Install using pip
```bash
$ cd vsomm2-modeler
$ pip install .
```

For GROMACS files, add the next line to your *.bashrc*:
```bash
export GROMACS_BIN="path_to_gromacs/bin/"
```
