In the following subdirectories are the input and output of DFT and conformer calculations for this publication:

XX

sub-directories:

* intermediate_calculations/
    * Input intermediate structures, incl. amine, aldehyde and water precursors, as `.mol` files.
    * Output structures from GFN2-xTB optimisation as `_opt.mol` and `_opt.xyz` files.
    * GFN2-xTB energies for different solvents as `_SOLV.ey` files.
    * `.csv` files containing:
        * all_total_energies.csv : Total energy of all intermediates from all methods.
        * all_formation_energies.csv : Formation energy of all intermediates from all methods.
        * all_formation_energies_per_imine.csv : Formation energy per imine bond formed of all intermediates from all methods.

* intermediate_calculations/gaussian_calculations/
	* Gaussian16 input and output for DFT calculations on the intermediates.
	* Structures are the GFN2-xTB optimised structures.
	* File names match the directory above as such `NAME_opt_SOLV_METHOD.gau/log`
	    * SOLV include `gas`, `dcm` (dichloromethane) and `cfm` (chloroform).
	    * METHOD inclde `pbe` (PBE1PBE/Def2TZVP) and `mp2` (MP2/aug-cc-pVDZ).

* intermediate_calculations/orca_calculations/
	* Orca input and output for DFT calculations on the intermediates.
	* Structures for B97-3c optimisation are the GFN2-xTB optimised structures.
	    * Result for each intermediate are saved as `_o_NAME_opt_B97-3c.xyz`.
	* Structures for the MP2 single-point are the B97-3c output.
	* File names match the directory above as such `o_NAME_opt_METHOD.in/out`
	    * All calculations performed in gas phase.
	    * METHOD inclde `B97-3c` (composite method) and `MP2` (RI-MP2 cc-pVTZ).
