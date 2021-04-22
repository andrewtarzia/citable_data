from os.path import exists
import json

import stk

from atools import (
    get_organic_linkers, calculate_energy, read_gfnx2xtb_eyfile
)


def calculate_ligand_SE(
    org_ligs,
    sdict,
    opt_mol_dir,
    output_json,
):
    """
    Calculate the strain energy of each ligand in the cage.

    Parameters
    ----------
    org_lig : :class:`dict` of :class:`stk.BuildingBlock`
        Dictionary of building blocks where the key is the file name,
        and the value is the stk building block.

    output_json : :class:`str`
        File name to save output to to avoid reruns.

    Returns
    -------
    strain_energies : :class:`dict`
        Strain energies for each ligand.

    """

    strain_energies = {}
    binol_f = sdict['binol']
    long_f = sdict['long']
    binol_natoms = sdict['batoms']
    long_natoms = sdict['latoms']
    # Iterate over ligands.
    for lig in org_ligs:
        stk_lig = org_ligs[lig]
        ey_file = lig.replace('mol', 'ey')

        # Calculate energy of extracted ligand.
        if not exists(ey_file):
            calculate_energy(
                name=lig.replace('.mol', ''),
                mol=stk_lig,
                ey_file=ey_file
            )
        # Read energy.
        # kJ/mol.
        E_extracted = read_gfnx2xtb_eyfile(ey_file)

        # Is binol or is long?
        sgt = stk_lig.get_num_atoms()
        if sgt == binol_natoms:
            opt_lig_n = binol_f+'_dopt'
            opt_filename = binol_f+'_dopt'
            opt_lig_ey = binol_f+'_dopt.ey'
        elif sgt == long_natoms:
            opt_lig_n = long_f+'_dopt'
            opt_filename = long_f+'_dopt'
            opt_lig_ey = long_f+'_dopt.ey'

        # Calculate energy of optimised ligand.
        # Load in lowest energy conformer.
        opt_mol = stk.BuildingBlock.init_from_file(
            f'{opt_mol_dir}/{opt_filename}.mol'
        )
        if not exists(opt_lig_ey):
            calculate_energy(
                name=opt_lig_n,
                mol=opt_mol,
                ey_file=opt_lig_ey
            )
        # Read energy.
        # kJ/mol.
        E_free = read_gfnx2xtb_eyfile(opt_lig_ey)
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
        for lig in org_ligs:
            stk_lig = org_ligs[lig]
            stk_lig.write(lig.replace('.mol', '.xyz'))
        lse_dict = calculate_ligand_SE(
            org_ligs=org_ligs,
            sdict=sdict,
            opt_mol_dir=low_c_dir,
            output_json=f'{struct}.json',
        )
        all_data[struct] = lse_dict

    print(all_data)


if __name__ == '__main__':
    main()