In the following subdirectories are the input and outputs of cage and face analysis for:

chemrxiv: 

Published: 10.1039/d2sc03856k

Data DOI:

.. image:: https://zenodo.org/badge/DOI/10.5281/zenodo.6795512.svg
   :target: https://doi.org/10.5281/zenodo.6795512

Associated code DOI:

.. image:: https://zenodo.org/badge/DOI/10.5281/zenodo.6795472.svg
   :target: https://doi.org/10.5281/zenodo.6795472

cage_library directory:
    * contains the structures from xTB optimisation as "_opt.mol"
    
complex_library directory:
    * contains the optimised structures of both complexes


ligand_library directory:
    * contains `_opt.mol` input ligand structures for cage construction
    * for cap, the input was provided manually in `manual/` directory
    * in `face_analysis` directory:
        * contains manual_complex directory, with necessary input for face construction
        * _long_properties.json files contains the measurements for the named face (in file name)
        * _long_lopt.mol files contain the optimised structure of the named face, on which analysis was performed
        * `long` corresponds to the longer restricted optimization discussed in the SI.


xray_structures directory:
    * analysis directory:
        * contains input .pdb files for xray structure (as single molecules) used in analysis
        * contains `all_xray_csv_data.csv`, which has all data needed on xray structures.
