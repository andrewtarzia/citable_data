import stk
import stko
import os
import pathlib
import json
import numpy as np


def gulp_path():
    return "/home/atarzia/software/gulp-6.1/Src/gulp"


def xtb_path():
    return "/home/atarzia/miniconda3/envs/simple_het/bin/xtb"


def optimisation_sequence(mol, name, charge, calc_dir):
    gulp1_output = os.path.join(calc_dir, f"{name}_gulp1.mol")
    gulp2_output = os.path.join(calc_dir, f"{name}_gulp2.mol")
    gulpmd_output = os.path.join(calc_dir, f"{name}_gulpmd.mol")
    xtbopt_output = os.path.join(calc_dir, f"{name}_xtb.mol")

    if not os.path.exists(gulp1_output):
        output_dir = os.path.join(calc_dir, f"{name}_gulp1")
        CG = True
        print(f"UFF4MOF optimisation 1 of {name} CG: {CG}")
        gulp_opt = stko.GulpUFFOptimizer(
            gulp_path=gulp_path(),
            maxcyc=1000,
            metal_FF={46: "Pd4+2"},
            metal_ligand_bond_order="",
            output_dir=output_dir,
            conjugate_gradient=CG,
        )
        gulp_opt.assign_FF(mol)
        gulp1_mol = gulp_opt.optimize(mol=mol)
        gulp1_mol.write(gulp1_output)
    else:
        print(f"loading {gulp1_output}")
        gulp1_mol = mol.with_structure_from_file(gulp1_output)

    if not os.path.exists(gulp2_output):
        output_dir = os.path.join(calc_dir, f"{name}_gulp2")
        CG = False
        print(f"UFF4MOF optimisation 2 of {name} CG: {CG}")
        gulp_opt = stko.GulpUFFOptimizer(
            gulp_path=gulp_path(),
            maxcyc=1000,
            metal_FF={46: "Pd4+2"},
            metal_ligand_bond_order="",
            output_dir=output_dir,
            conjugate_gradient=CG,
        )
        gulp_opt.assign_FF(gulp1_mol)
        gulp2_mol = gulp_opt.optimize(mol=gulp1_mol)
        gulp2_mol.write(gulp2_output)
    else:
        print(f"loading {gulp2_output}")
        gulp2_mol = mol.with_structure_from_file(gulp2_output)

    if not os.path.exists(gulpmd_output):
        output_dir = os.path.join(calc_dir, f"{name}_gulpmd")

        print(f"UFF4MOF equilib MD of {name}")
        gulp_MD = stko.GulpUFFMDOptimizer(
            gulp_path=gulp_path(),
            metal_FF={46: "Pd4+2"},
            metal_ligand_bond_order="",
            output_dir=os.path.join(calc_dir, f"{name}_gulpmd"),
            integrator="leapfrog verlet",
            ensemble="nvt",
            temperature=1000,
            timestep=0.25,
            equilbration=0.5,
            production=0.5,
            N_conformers=2,
            opt_conformers=False,
            save_conformers=False,
        )
        gulp_MD.assign_FF(gulp2_mol)
        gulpmd_mol = gulp_MD.optimize(mol=gulp2_mol)

        print(f"UFF4MOF production MD of {name}")
        gulp_MD = stko.GulpUFFMDOptimizer(
            gulp_path=gulp_path(),
            metal_FF={46: "Pd4+2"},
            metal_ligand_bond_order="",
            output_dir=os.path.join(calc_dir, f"{name}_gulpmd"),
            integrator="leapfrog verlet",
            ensemble="nvt",
            temperature=1000,
            timestep=0.75,
            equilbration=0.5,
            production=200.0,
            N_conformers=40,
            opt_conformers=True,
            save_conformers=False,
        )
        gulp_MD.assign_FF(gulpmd_mol)
        gulpmd_mol = gulp_MD.optimize(mol=gulpmd_mol)
        gulpmd_mol.write(gulpmd_output)
    else:
        print(f"loading {gulpmd_output}")
        gulpmd_mol = mol.with_structure_from_file(gulpmd_output)

    if not os.path.exists(xtbopt_output):
        output_dir = os.path.join(calc_dir, f"{name}_xtbopt")
        print(f"final xtb optimisation of {name}")
        xtb_opt = stko.XTB(
            xtb_path=xtb_path(),
            output_dir=output_dir,
            gfn_version=2,
            num_cores=6,
            charge=charge,
            opt_level="extreme",
            num_unpaired_electrons=0,
            max_runs=1,
            calculate_hessian=False,
            unlimited_memory=True,
            solvent="DMSO",
            solvent_grid="tight",
            solvent_model="alpb",
        )
        xtbopt_mol = xtb_opt.optimize(mol=gulpmd_mol)
        xtbopt_mol.write(xtbopt_output)
    else:
        print(f"loading {xtbopt_output}")
        xtbopt_mol = mol.with_structure_from_file(xtbopt_output)

    final_mol = mol.with_structure_from_file(xtbopt_output)
    return final_mol


def host_optimisation_sequence(mol, name, charge, calc_dir, opt_level=None):
    if opt_level is None:
        opt_level = "extreme"

    xtbopt_output = os.path.join(calc_dir, f"{name}_xtb.mol")

    if not os.path.exists(xtbopt_output):
        output_dir = os.path.join(calc_dir, f"{name}_xtbopt")
        print(f"xtb optimisation of {name}")
        xtb_opt = stko.XTB(
            xtb_path=xtb_path(),
            output_dir=output_dir,
            gfn_version=2,
            num_cores=6,
            charge=charge,
            opt_level=opt_level,
            num_unpaired_electrons=0,
            max_runs=1,
            calculate_hessian=False,
            unlimited_memory=True,
            solvent="DMSO",
            solvent_grid="tight",
            solvent_model="alpb",
        )
        xtbopt_mol = xtb_opt.optimize(mol=mol)
        xtbopt_mol.write(xtbopt_output)
    else:
        print(f"loading {xtbopt_output}")
        xtbopt_mol = mol.with_structure_from_file(xtbopt_output)

    return mol.with_structure_from_file(xtbopt_output)


def get_xtb_energies(mol, name, charge, calc_dir):
    xtb_output = os.path.join(calc_dir, f"{name}_xtb.ey")
    output_dir = os.path.join(calc_dir, f"{name}_xtbhessey")
    if os.path.exists(xtb_output):
        # Read .ey file.
        with open(xtb_output, "r") as f:
            res_dict = json.load(f)

    else:
        print(f"running free energy of {name}")
        xtb = stko.XTBEnergy(
            xtb_path=xtb_path(),
            output_dir=output_dir,
            gfn_version=2,
            charge=charge,
            unlimited_memory=False,
            num_cores=1,
            calculate_free_energy=True,
            solvent="DMSO",
            solvent_grid="tight",
            solvent_model="alpb",
        )
        xtb_results = xtb.get_results(mol)
        total_energy = xtb_results.get_total_energy()
        total_free_energy = xtb_results.get_total_free_energy()
        res_dict = {
            "totale": total_energy[0],
            "totale_unit": total_energy[1],
            "fe": total_free_energy[0],
            "fe_unit": total_free_energy[1],
        }
        # Save to .ey file.
        with open(xtb_output, "w") as f:
            json.dump(res_dict, f)

    return res_dict


class AromaticCNCFactory(stk.FunctionalGroupFactory):
    """
    A subclass of stk.SmartsFunctionalGroupFactory.

    """

    def __init__(self, bonders=(1,), deleters=()):
        """
        Initialise :class:`.AromaticCNCFactory`.

        """

        self._bonders = bonders
        self._deleters = deleters

    def get_functional_groups(self, molecule):
        generic_functional_groups = stk.SmartsFunctionalGroupFactory(
            smarts="[#6]~[#7X2]~[#6]",
            bonders=self._bonders,
            deleters=self._deleters,
        ).get_functional_groups(molecule)
        for fg in generic_functional_groups:
            atom_ids = (i.get_id() for i in fg.get_atoms())
            atoms = tuple(molecule.get_atoms(atom_ids))
            yield AromaticCNC(
                carbon1=atoms[0],
                nitrogen=atoms[1],
                carbon2=atoms[2],
                bonders=tuple(atoms[i] for i in self._bonders),
                deleters=tuple(atoms[i] for i in self._deleters),
            )


class AromaticCNC(stk.GenericFunctionalGroup):
    """
    Represents an N atom in pyridine functional group.

    The structure of the functional group is given by the pseudo-SMILES
    ``[carbon][nitrogen][carbon]``.

    """

    def __init__(self, carbon1, nitrogen, carbon2, bonders, deleters):
        """
        Initialize a :class:`.Alcohol` instance.

        Parameters
        ----------
        carbon1 : :class:`.C`
            The first carbon atom.

        nitrogen : :class:`.N`
            The nitrogen atom.

        carbon2 : :class:`.C`
            The second carbon atom.

        bonders : :class:`tuple` of :class:`.Atom`
            The bonder atoms.

        deleters : :class:`tuple` of :class:`.Atom`
            The deleter atoms.

        """

        self._carbon1 = carbon1
        self._nitrogen = nitrogen
        self._carbon2 = carbon2
        atoms = (carbon1, nitrogen, carbon2)
        super().__init__(atoms, bonders, deleters)

    def get_carbon1(self):
        """
        Get the first carbon atom.

        Returns
        -------
        :class:`.C`
            The first carbon atom.

        """

        return self._carbon1

    def get_carbon2(self):
        """
        Get the second carbon atom.

        Returns
        -------
        :class:`.C`
            The second carbon atom.

        """

        return self._carbon2

    def get_nitrogen(self):
        """
        Get the nitrogen atom.

        Returns
        -------
        :class:`.N`
            The nitrogen atom.

        """

        return self._nitrogen

    def clone(self):
        clone = super().clone()
        clone._carbon1 = self._carbon1
        clone._nitrogen = self._nitrogen
        clone._carbon2 = self._carbon2
        return clone

    def with_atoms(self, atom_map):
        clone = super().with_atoms(atom_map)
        clone._carbon1 = atom_map.get(
            self._carbon1.get_id(),
            self._carbon1,
        )
        clone._nitrogen = atom_map.get(
            self._nitrogen.get_id(),
            self._nitrogen,
        )
        clone._carbon2 = atom_map.get(
            self._carbon2.get_id(),
            self._carbon2,
        )
        return clone

    def __repr__(self):
        return (
            f"{self.__class__.__name__}("
            f"{self._carbon1}, {self._nitrogen}, {self._carbon2}, "
            f"bonders={self._bonders})"
        )


def add_cl_opt(
    cage_name,
    isomer_name,
    cage,
    structure_dir,
    calculation_dir,
    energies,
):
    host_name = f"{cage_name}_{isomer_name}_cl"
    guest = stk.host_guest.Guest(
        building_block=stk.BuildingBlock("[Cl]"),
    )
    host_guest = stk.ConstructedMolecule(
        topology_graph=stk.host_guest.Complex(
            host=stk.BuildingBlock.init_from_molecule(cage),
            guests=guest,
        ),
    )
    host_guest.write(structure_dir / f"{host_name}_unopt.mol")

    # Optimise structures.
    host_guest = host_optimisation_sequence(
        mol=host_guest,
        name=host_name,
        charge=3,
        calc_dir=calculation_dir,
    )
    energies[cage_name][f"{isomer_name}_cl"] = get_xtb_energies(
        mol=host_guest,
        name=host_name,
        charge=3,
        calc_dir=calculation_dir,
    )
    host_guest.write(structure_dir / f"{host_name}_optc.mol")
    return host_guest, energies


def add_bf4_to_cl_opt(
    cage_name,
    isomer_name,
    host_guest,
    structure_dir,
    calculation_dir,
    energies,
):
    host_name = f"{cage_name}_{isomer_name}_cl_bf4"

    # Get two guest positions.
    centroid = host_guest.get_centroid()
    shift = 5
    ion_positions = []
    pd_vectors = []
    for atom in host_guest.get_atoms():
        x, y, z = next(host_guest.get_atomic_positions(atom.get_id()))
        if atom.get_atomic_number() == 46:
            # Shift xyz further out.
            xyz = np.array((x, y, z))
            vector = xyz - centroid
            vector = vector / np.linalg.norm(vector)
            xyz = xyz + vector * shift
            ion_positions.append(xyz)
            pd_vector = xyz - np.array((x, y, z))
            pd_vectors.append(pd_vector)

    assert len(ion_positions) == 2
    bf4_guest1 = stk.host_guest.Guest(
        building_block=stk.BuildingBlock("[B-](F)(F)(F)F"),
        displacement=ion_positions[0],
        end_vector=-pd_vectors[0],
    )
    bf4_guest2 = stk.host_guest.Guest(
        building_block=stk.BuildingBlock("[B-](F)(F)(F)F"),
        displacement=ion_positions[1],
        end_vector=-pd_vectors[1],
    )
    host_guest = stk.ConstructedMolecule(
        topology_graph=stk.host_guest.Complex(
            host=stk.BuildingBlock.init_from_molecule(host_guest),
            guests=(bf4_guest1, bf4_guest2),
        ),
    )
    host_guest.write(structure_dir / f"{host_name}_unopt.mol")

    # Optimise structures.
    host_guest = host_optimisation_sequence(
        mol=host_guest,
        name=host_name,
        charge=1,
        calc_dir=calculation_dir,
        opt_level="extreme",
    )
    energies[cage_name][f"{isomer_name}_cl_bf4"] = get_xtb_energies(
        mol=host_guest,
        name=host_name,
        charge=1,
        calc_dir=calculation_dir,
    )
    host_guest.write(structure_dir / f"{host_name}_optc.mol")
    return host_guest, energies


def get_gulp1_cage(mol, name, calc_dir):
    gulp1_output = os.path.join(calc_dir, f"{name}_gulp1.mol")

    if not os.path.exists(gulp1_output):
        raise FileNotFoundError(f"{gulp1_output} not found!")
    else:
        print(f"loading {gulp1_output}")
        gulp1_mol = mol.with_structure_from_file(gulp1_output)
    return gulp1_mol.with_centroid(np.array((0, 0, 0)))


def add_cl_and_bf4_opt(
    cage_name,
    isomer_name,
    cage,
    structure_dir,
    calculation_dir,
    energies,
):
    gulp_cage = get_gulp1_cage(
        mol=cage,
        name=f"{cage_name}_{isomer_name}",
        calc_dir=calculation_dir,
    )

    host_name = f"{cage_name}_{isomer_name}_g1_cl_bf4"

    # Get two guest positions.
    centroid = gulp_cage.get_centroid()
    shift = 5
    ion_positions = []
    pd_vectors = []
    for atom in gulp_cage.get_atoms():
        x, y, z = next(gulp_cage.get_atomic_positions(atom.get_id()))
        if atom.get_atomic_number() == 46:
            # Shift xyz further out.
            xyz = np.array((x, y, z))
            vector = xyz - centroid
            vector = vector / np.linalg.norm(vector)
            xyz = xyz + vector * shift
            ion_positions.append(xyz)
            pd_vector = xyz - np.array((x, y, z))
            pd_vectors.append(pd_vector)

    assert len(ion_positions) == 2
    bf4_guest1 = stk.host_guest.Guest(
        building_block=stk.BuildingBlock("[B-](F)(F)(F)F"),
        displacement=ion_positions[0],
        end_vector=-pd_vectors[0],
    )
    bf4_guest2 = stk.host_guest.Guest(
        building_block=stk.BuildingBlock("[B-](F)(F)(F)F"),
        displacement=ion_positions[1],
        end_vector=-pd_vectors[1],
    )
    cl_guest = stk.host_guest.Guest(
        building_block=stk.BuildingBlock("[Cl]"),
        displacement=(ion_positions[0] + ion_positions[1]) / 2,
    )
    host_guest = stk.ConstructedMolecule(
        topology_graph=stk.host_guest.Complex(
            host=stk.BuildingBlock.init_from_molecule(gulp_cage),
            guests=(bf4_guest1, bf4_guest2, cl_guest),
        ),
    )
    host_guest.write(structure_dir / f"{host_name}_unopt.mol")
    # Optimise structures.
    host_guest = host_optimisation_sequence(
        mol=host_guest,
        name=host_name,
        charge=1,
        calc_dir=calculation_dir,
        opt_level="extreme",
    )
    energies[cage_name][f"{isomer_name}_g1_cl_bf4"] = get_xtb_energies(
        mol=host_guest,
        name=host_name,
        charge=1,
        calc_dir=calculation_dir,
    )
    host_guest.write(structure_dir / f"{host_name}_optc.mol")
    return host_guest, energies


def main():

    structure_dir = pathlib.Path("./structures_alpb")
    if not structure_dir.exists():
        os.mkdir(structure_dir)

    calculation_dir = pathlib.Path("./calculations_alpb")
    if not calculation_dir.exists():
        os.mkdir(calculation_dir)

    ligands = {
        "l1q": stk.BuildingBlock(
            smiles="C1C=C(COC(=O)C2C=C3C(C=CC=C3)=NC=2)C=NC=1",
            functional_groups=(AromaticCNCFactory(),),
        ),
        "l1p": stk.BuildingBlock(
            smiles="C1C=C(COC(=O)C2C=CC(C)=NC=2)C=NC=1",
            functional_groups=(AromaticCNCFactory(),),
        ),
        "l2q": stk.BuildingBlock(
            smiles=("C1C=C(C#CC2C=C(C#CC3C=NC4=C(C=CC=C4)C=3)C=CC=2)C=NC=1"),
            functional_groups=(AromaticCNCFactory(),),
        ),
        "tq": stk.BuildingBlock(
            smiles="N1C=C(C2C=CC=C(C3C=CC=NC=3)C=2)C=C2C=1C=CC=C2",
            functional_groups=(AromaticCNCFactory(),),
        ),
        "tp": stk.BuildingBlock(
            smiles="N1C=C(C2=CC(C3C=CC=NC=3)=CC=C2)C=CC=1C",
            functional_groups=(AromaticCNCFactory(),),
        ),
    }

    pd = stk.BuildingBlock(
        smiles="[Pd+2]",
        functional_groups=(
            stk.SingleAtom(stk.Pd(0, charge=2)) for i in range(4)
        ),
        position_matrix=np.array([[0, 0, 0]]),
    )

    energies = {}
    for i in ligands:
        energies[i] = {}

        opt_file = structure_dir / f"{i}_opt.mol"
        # Use manual orienation.
        ligands[i] = stko.UFF().optimize(ligands[i])
        ligands[i].write(structure_dir / f"{i}.mol")
        if not opt_file.exists():
            raise FileNotFoundError(
                "You need to manually optimise ligand directions."
            )
        ligands[i] = ligands[i].with_structure_from_file(str(opt_file))
        ligands[i].write(opt_file)

        # Construct Isomers.
        c_orientations = {
            "A": {2: 0, 3: 0, 4: 0, 5: 0},
            "B": {2: 1, 3: 0, 4: 0, 5: 0},
            "C": {2: 1, 3: 1, 4: 0, 5: 0},
            "D": {2: 1, 3: 0, 4: 1, 5: 0},
        }

        for c_a in c_orientations:
            c_orientation = c_orientations[c_a]
            name_ = f"{i}_{c_a}"

            # Build cage.
            cage = stk.ConstructedMolecule(
                topology_graph=stk.cage.M2L4Lantern(
                    building_blocks=(pd, ligands[i]),
                    optimizer=stk.MCHammer(),
                    vertex_alignments=c_orientation,
                    reaction_factory=stk.DativeReactionFactory(
                        stk.GenericReactionFactory(
                            bond_orders={
                                frozenset({AromaticCNC, stk.SingleAtom}): 9,
                            }
                        )
                    ),
                ),
            )
            cage.write(structure_dir / f"{name_}_unopt.mol")

            # Optimise structures.
            cage = optimisation_sequence(
                mol=cage,
                name=name_,
                charge=4,
                calc_dir=calculation_dir,
            )
            energies[i][c_a] = get_xtb_energies(
                mol=cage,
                name=name_,
                charge=4,
                calc_dir=calculation_dir,
            )
            cage.write(structure_dir / f"{name_}_optc.mol")

            if "l1" in i:
                host_guest, energies = add_cl_opt(
                    cage_name=i,
                    isomer_name=c_a,
                    cage=cage,
                    structure_dir=structure_dir,
                    calculation_dir=calculation_dir,
                    energies=energies,
                )

                # Now add two BF4 ions on the outside
                host_guest, energies = add_bf4_to_cl_opt(
                    cage_name=i,
                    isomer_name=c_a,
                    host_guest=host_guest,
                    structure_dir=structure_dir,
                    calculation_dir=calculation_dir,
                    energies=energies,
                )

                # Want to do ion-based xtb opt from first Gulp structure too,
                # before twist comes into structure.
                host_guest, energies = add_cl_and_bf4_opt(
                    cage_name=i,
                    isomer_name=c_a,
                    cage=cage,
                    structure_dir=structure_dir,
                    calculation_dir=calculation_dir,
                    energies=energies,
                )

    with open("xtb_energies.json", "w") as f:
        json.dump(energies, f, indent=4)
    print(energies)


if __name__ == "__main__":
    main()
