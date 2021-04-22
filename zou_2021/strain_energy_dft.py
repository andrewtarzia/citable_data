from os.path import exists
import json
import sys
import pandas as pd

import stk

from atools import get_organic_linkers, calculate_energy


def calculate_ligand_SE(
    org_ligs,
    sdict,
    extracted_energies,
    free_energies,
    output_json,
):

    strain_energies = {}
    binol_f = sdict['binol']
    long_f = sdict['long']
    binol_natoms = sdict['batoms']
    long_natoms = sdict['latoms']
    # Iterate over ligands.
    for lig in org_ligs:
        stk_lig = org_ligs[lig]
        # kJ/mol.
        E_extracted = extracted_energies[lig] * 2625.5

        # Is binol or is long?
        sgt = stk_lig.get_num_atoms()
        if sgt == binol_natoms:
            opt_name = binol_f
        elif sgt == long_natoms:
            opt_name = long_f

        # kJ/mol.
        E_free = free_energies[opt_name] * 2625.5
        # Add to list the strain energy:
        # (E(extracted) - E(optimised/free))
        lse = E_extracted - E_free
        # kJ/mol.
        strain_energies[lig] = lse

    # Write data.
    with open(output_json, 'w') as f:
        json.dump(strain_energies, f)

    return strain_energies



def main():

    csv_file = 'low_e_dft_spe.csv'

    data = pd.read_csv(csv_file)
    free_ligs = ['longC', 'longD', 'binolA', 'binolB']
    free_ey_dict = {
        i: float(list(data[data['structure'] == i]['energy'])[0])
        for i in free_ligs
    }
    print(data)

    structures = {
        'fe4dd': {
            'binol': 'binolB', 'long': 'longC',
            'batoms': 86, 'latoms': 46,
        },
        'fe4ll': {
            'binol': 'binolB', 'long': 'longC',
            'batoms': 86, 'latoms': 46,
        },
        'fe6dd': {
            'binol': 'binolA', 'long': 'longD',
            'batoms': 80, 'latoms': 56,
        },
        'fe6ll': {
            'binol': 'binolA', 'long': 'longD',
            'batoms': 80, 'latoms': 56,
        },
        'fe7dd': {
            'binol': 'binolB', 'long': 'longD',
            'batoms': 86, 'latoms': 56,
        },
        'fe7ll': {
            'binol': 'binolB', 'long': 'longD',
            'batoms': 86, 'latoms': 56,
        },
    }
    cage_dir = 'orca_opt'
    low_c_dir = 'low_e_bb_confs'
    type_suffix = '_xtb_dft'
    metal_atom_no = 26

    all_data = {}
    for struct in structures:
        filename = f'{cage_dir}/{struct}{type_suffix}.mol'
        sdict = structures[struct]
        stk_mol = stk.BuildingBlock.init_from_file(filename)

        org_ligs, _ = get_organic_linkers(
            cage=stk_mol,
            metal_atom_nos=(metal_atom_no, ),
            file_prefix=f'{struct}{type_suffix}_sg'
        )

        org_ligs_energies = {}
        for i in org_ligs:
            dd = data[data['structure'] == i.replace('.mol', '')]
            org_ligs_energies[i] = float(list(dd['energy'])[0])

        print(org_ligs)
        print(org_ligs_energies)

        lse_dict = calculate_ligand_SE(
            org_ligs=org_ligs,
            sdict=sdict,
            extracted_energies=org_ligs_energies,
            free_energies=free_ey_dict,
            output_json=f'{struct}_dftstrain.json',
        )
        all_data[struct] = lse_dict

    print(all_data)
    for i in all_data:
        print(i, sum(all_data[i].values()))


if __name__ == '__main__':
    main()