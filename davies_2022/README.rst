In the following subdirectories are the input and outputs of cage and face analysis for:

chemrxiv: 

Published:

Data DOI:

.. image:: https://zenodo.org/badge/DOI/10.5281/zenodo.7328029.svg
   :target: https://doi.org/10.5281/zenodo.7328029

Associated code DOI:

.. image:: https://zenodo.org/badge/DOI/10.5281/zenodo.7328036.svg
   :target: https://doi.org/10.5281/zenodo.7328036

NOTES:
    * the naming convention differs from manuscript:


.. list-table:: naming convention
   :widths: 25 25 50
   :header-rows: 1

   * - manuscript tetra-aniline
     - computational label
     - xtal-labe;
   * - A 
     - 5
     - 370
   * - B
     - 16
     - 326
   * - C
     - 12
     - 235
   * - D
     - 3
     - 301
   * - E
     - 8
     - 257
   * - F
     - 2
     - 354

* computational labels are often preceded by `quad2_` or `cl1_quad2_`
* much of the analysis was not used in the manuscript but remains part of the accumulated data

cage_library directory:
    * _CS.json: information on all cages in the set of diastereomers - properties and whether they optimized successfully.
    * _ligand_measures.json: information on the ligand associated with a set of cage diastereomers.
    * _measures.json: represenets a cleaned up collation of all measures the diastereomers made from a given ligand
    * C_NAME_optc.mol: optimized (at xTB level) structure of each cage.
    * set_dft_run directory contains the input and output of the CP2K optimisations of one set of diastereomers

        
    
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
