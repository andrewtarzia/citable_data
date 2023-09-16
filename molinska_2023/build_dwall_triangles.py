import stk
import stko
import os
import pathlib
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


def mod_optimisation_sequence(mol, name, charge, calc_dir):
    gulp1_output = os.path.join(calc_dir, f"{name}_gulp1.mol")
    gulp2_output = os.path.join(calc_dir, f"{name}_gulp2.mol")
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
        xtbopt_mol = xtb_opt.optimize(mol=gulp2_mol)
        xtbopt_mol.write(xtbopt_output)
    else:
        print(f"loading {xtbopt_output}")
        xtbopt_mol = mol.with_structure_from_file(xtbopt_output)

    final_mol = mol.with_structure_from_file(xtbopt_output)
    return final_mol


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


def main():

    structure_dir = pathlib.Path("./structures_alpb")
    if not structure_dir.exists():
        os.mkdir(structure_dir)

    calculation_dir = pathlib.Path("./calculations_alpb")
    if not calculation_dir.exists():
        os.mkdir(calculation_dir)

    ligands = {
        "c3p": stk.BuildingBlock(
            smiles="C1C=CC(C(OC2C=C(OC(C3C=NC(C)=CC=3)=O)C=CC=2)=O)=CN=1",
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

    for i in ligands:
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
            "cis": {3: 1, 4: 1, 5: 0, 6: 0, 7: 0, 8: 0},
            "trans": {3: 0, 4: 1, 5: 1, 6: 0, 7: 1, 8: 0},
        }

        for c_a in c_orientations:
            c_orientation = c_orientations[c_a]
            name_ = f"{i}_{c_a}"

            # Build cage.
            cage = stk.ConstructedMolecule(
                topology_graph=stk.cage.M3L6(
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
                charge=6,
                calc_dir=calculation_dir,
            )
            cage.write(structure_dir / f"{name_}_optc.mol")

            # Load manual structure, and optimise again.
            mod_name_ = f"{i}_{c_a}_mod"
            cage = cage.with_structure_from_file(
                str(structure_dir / f"{name_}_mod.mol")
            )
            cage = mod_optimisation_sequence(
                mol=cage,
                name=mod_name_,
                charge=6,
                calc_dir=calculation_dir,
            )
            cage = cage.write(structure_dir / f"{name_}_mod_optc.mol")


if __name__ == "__main__":
    main()
