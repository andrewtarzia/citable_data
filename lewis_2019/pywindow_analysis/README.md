structures:
-----------

* from DFT calculations:
    * 1_D
    * 4/5/6_C
* from XRD:
    * 1492902.cif (DOI: 10.5517/ccdc.csd.cc1m3h4c)
    * 768969.cif (DOI: 10.5517/cctt5g7)

method:
-------

* note that the majority of this code was taken from my github repo: https://github.com/andrewtarzia/atools
* the following was done using Python 3.7, pywindow 0.0.3, the atomic simulation environment (ASE) 3.17.0
* script 1 (single_molecule.py):
    * takes an XYZ file of a single molecule and performs full pywindow analysis.
    * performed on DFT structures
    * produces {name}.json and {name}.pdb (with atoms placed at pore COM and window COM) for each {name}.xyz

* script 2 (XRD.py)
    * takes a PDB of the unit cell in a CIF and extracts all individual molecules in that unit cell
    * for each unit cell, it performs a full analysis from pywindow
    * produces {name}_{molecule_no}.json {name}_{molecule_no}.pdb (with atoms placed at pore COM and window COM) for each molecule in {name}.cif
