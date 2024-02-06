In the following subdirectories are the input and output of GFN2-xTB and DFT calculations for this publication:

chemrxiv: XXX

Published: XXXX

Data DOI:

ZZZ

Associated code DOI:

ZZZ

Code repository: https://github.com/andrewtarzia/simple_het_construction


data directory:
    * a spreadsheet with all final energy values and exchange energy calculations
    * CSD Survey data, NPd_survey_data_261119.csv, for Pd centres




Naming convention for file conversions:

- l1: 1DBF
- l2: 1Ph
- l3: 1Th
- la: 2DBF
- lb: 2Py
- lc: 2Ph
- ld: 2Th

Structure naming convention: 
    `mX` implies homoleptic cage with X Pd atoms, `cis`/`trans` are the cis/trans heteroleptic cages, respectively

structures/xtb directory:
    * contains the structures from GFN2-xTB/ALPB(DMSO) optimisations of stk-generated structures   

structures/opt_*METHOD*_SP_*METHOD*_06-02-2024 directories:
    * All DFT was run by Victor Posligua
    * contains the input files (.com), output files (.log) and structure files (.xyz/.mol) of DFT optimisations and single point energy calculations with each method
    * you’ll find 8 different folders:
        opt_PBE0_SP_PBE0_06-02-2024
        opt_PBE0_SP_B3LYP_06-02-2024
        opt_PBE0_SP_B97D3_06-02-2024
        opt_PBE0_SP_HSE_06-02-2024
        opt_B3LYP_SP_B3LYP_06-02-2024
        opt_B97D3_SP_B97D3_06-02-2024
        opt_HSE_SP_HSE_06-02-2024
        opt_GFN2-xTB_SP_PBE0_06-02-2024
    * When the opt method and SP method are the same, the final structure is included in .mol and .xyz formats
    * However, if opt method is different from the SP method, the final structure is not included because only a single-point energy calculation was run. 
    * For example, there are no .mol or .xyz files for 'opt_PBE0_SP_B3LYP_06-02-2024’ since the structure is already in 'opt_PBE0_SP_PBE0_06-02-2024’.

