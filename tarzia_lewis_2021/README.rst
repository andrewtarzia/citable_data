In the following subdirectories are the input and output of Gaussian calculations + structures from screening for this publication:

chemrxiv: 10.26434/chemrxiv.14604294

Published: 10.1002/anie.202106721

Data DOI:

.. image:: https://zenodo.org/badge/204879639.svg
   :target: https://zenodo.org/badge/latestdoi/204879639

Associated code DOI:

.. image:: https://zenodo.org/badge/190806715.svg
   :target: https://zenodo.org/badge/latestdoi/190806715


screening_structures directory:
    * contains the structures from xTB optimisation that were used in ranking in structures.tar.gz as `.mol` files
    * all_cage_results.txt contains their properties for ranking

single_point_dft directory:
    * contains the structures and output of single point DFT calculations
    * during the revision process, we confirmed (based on reviewer suggestions) our DFT validation results using ORCA 4.2.1 with PBE0 and B97-3c in the gas phase. These results were consistent with our previous ones, so were not added to the manuscript. But are useful for future work!
        * These results are in the s_orca directory.
        
free_energy_calculations directory:
    * during the revision process, it was suggested to calculate the free energies using the xTB method (low-cost) and compare that to the total energies we use.
    * the script `run_gfn2_free_energy.py` in the unsymm_match code repository does this for top candidate ligands using the stko.XTB class.
        * for each structure, the free energy is output to a .fey file.
