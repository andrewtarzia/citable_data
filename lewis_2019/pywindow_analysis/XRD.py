#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Script to run pywindow on a series of CIFs.

Author: Andrew Tarzia

Date Created: 31 Oct 2019

"""

import glob
import pywindow as pw
from os.path import isfile
from ase.io import read


def convert_CIF_2_PDB(file, wstruct=True):
    """
    Convert CIF to PDB file, save and return structure.

    """
    pdb_file = file.replace('.cif', '.pdb')
    print('converting:', file, 'to', pdb_file)
    if isfile(pdb_file) is False:
        try:
            structure = read(file)
        except IndexError:
            print('ASE load failed with IndexError. Skipping...')
            if wstruct:
                return None, None
            else:
                return None
        except ValueError:
            print('ASE load failed with IndexError. Skipping...')
            if wstruct:
                return None, None
            else:
                return None
        structure.write(pdb_file)
        print('conversion done.')
    structure = read(pdb_file)
    if wstruct:
        return pdb_file, structure
    else:
        return pdb_file


def rebuild(file, overwrite=False):
    """
    Rebuild the PDB system, output and reread.

    """
    out_file = file.replace('.pdb', '_rebuild.pdb')
    if isfile(out_file) is False or overwrite is True:
        print('rebuilding:', file)
        molsys = pw.MolecularSystem.load_file(file)
        rebuild_molsys = molsys.rebuild_system()
        # output
        rebuild_molsys.dump_system(out_file,
                                   include_coms=False,
                                   override=True)
        print('rebuild done.')
    else:
        rebuild_molsys = pw.MolecularSystem.load_file(out_file)
    return rebuild_molsys


def modularize(file):
    """
    Rebuild pyWindow MolecularSystem from file and modularize.

    """
    rebuilt_structure = rebuild(file=file)
    if len(rebuilt_structure.system['coordinates']) == 0:
        print(
            f'WARN: {file} failed rebuild using pyWindow, return None'
        )
        return None
    rebuilt_structure.make_modular()
    return rebuilt_structure


def analyze_rebuilt(rebuilt_structure, file_prefix, atom_limit,
                    include_coms=True, verbose=False):
    """
    Run all desired analysis on each molecule in rebuilt structure.
        (modified version of Example6 of pywindow examples.)

    Keyword Arguments:
        rebuilt_structure (Molecule) - Pywindow rebuilt unitcell with
            cage molecules
        file_prefix (str) - file naming convention
        atom_limit (int) - number of atoms used as cutoff for analysis
        include_coms (bool) - whether output PDB includes window COMs
        verbose (bool) - [Default = False]

    Returns:
        result_dict (dictionary) - dictionary of window information
            for all cages

    """
    result_dict = {}
    for molecule in rebuilt_structure.molecules:
        mol = rebuilt_structure.molecules[molecule]
        if mol.no_of_atoms < atom_limit:
            continue
        print(
            'Analysing molecule {0} out of {1}'.format(
                molecule + 1, len(rebuilt_structure.molecules)
            )
        )
        try:
            analysis = mol.full_analysis()
        except ValueError:
            print(
                f'WARN: {file_prefix}_{molecule} failed pywindow '
                'full_analysis.'
            )
            continue
        if verbose:
            print(analysis, '\n')
        # Each molecule can be saved separately
        mol.dump_molecule(
            file_prefix + "_{0}.pdb".format(molecule),
            include_coms=include_coms,
            override=True)
        mol.dump_properties_json(
            file_prefix + "_{0}.json".format(molecule),
            override=True)
        # output COM, window size and COM, and COM of pore optimized
        result_dict[molecule] = (
            analysis['pore_diameter_opt']
        )
    # print(result_dict)
    print('analysis done.')
    return result_dict


def main():
    files = glob.glob(f'*.cif')
    for file in files:
        print(f'doing {file}')
        # Convert to PDB.
        pdb_file = convert_CIF_2_PDB(file, wstruct=False)
        # Modularize with pywindow into independant molecules.
        rebuilt = modularize(pdb_file)
        r_dict = analyze_rebuilt(
            rebuilt_structure=rebuilt,
            file_prefix=file.replace('.cif', '_out'),
            atom_limit=50,
            include_coms=True,
            verbose=False
        )
        print(r_dict)


if __name__ == "__main__":
    main()
