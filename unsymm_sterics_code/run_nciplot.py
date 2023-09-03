import stk
import os
import pathlib
import subprocess as sp
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl


def nci_path():
    return (
        "/home/atarzia/workingspace/unsymm-sterics/nciplot/nciplot/"
        "src_nciplot_4.2/nciplot"
    )


def plot_g(oname):
    dat_file = f"{oname}.dat"
    out_file = f"{oname}.png"
    data = np.loadtxt(dat_file).T

    print(len(data[0]))

    vmin = -0.06
    vmax = 0.06

    fig, ax = plt.subplots(figsize=(8, 5))

    ax.scatter(
        x=data[0],
        y=data[1],
        c=data[0],
        edgecolors="none",
        s=20,
        vmin=vmin,
        vmax=vmax,
        cmap="jet",
    )

    ax.tick_params(axis="both", which="major", labelsize=16)
    ax.set_xlabel(r"sign($\lambda_2$)$\rho$ (a.u.)", fontsize=16)
    ax.set_ylabel("RDG", fontsize=16)

    ax.set_xlim(-0.07, 0.07)
    ax.set_ylim(0, 1)

    cbar_ax = fig.add_axes([1.01, 0.2, 0.02, 0.7])
    cmap = mpl.cm.jet
    norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
    cbar = fig.colorbar(
        mpl.cm.ScalarMappable(norm=norm, cmap=cmap),
        cax=cbar_ax,
        orientation="vertical",
    )
    cbar.ax.tick_params(labelsize=16)
    cbar.set_label("", fontsize=16)

    fig.tight_layout()
    fig.savefig(
        out_file,
        dpi=120,
        bbox_inches="tight",
    )
    plt.close()


def run_nci(structure, name, calc_dir):

    if not os.path.exists(calc_dir):
        os.mkdir(calc_dir)
    init_dir = os.getcwd()
    os.chdir(calc_dir)

    range_rule = "RANGE 3\n-0.1 -0.015\n-0.015 0.015\n0.015 -0.1\n"

    runs = {
        "full": ("", range_rule),
        # "full_range": ("", range_rule),
    }

    radius = 2

    # Find Pd atoms for radii.
    centroid = structure.get_centroid()
    shift = 2
    xyz_string = []
    for atom in structure.get_atoms():
        x, y, z = next(structure.get_atomic_positions(atom.get_id()))
        if atom.get_atomic_number() == 46:
            # Shift xyz further out.
            xyz = np.array((x, y, z))
            vector = xyz - centroid
            vector = vector / np.linalg.norm(vector)
            xyz = xyz + vector * shift
            runs[str(atom.get_id())] = (
                (
                    f"RADIUS {round(xyz[0], 2)} {round(xyz[1], 2)} "
                    f"{round(xyz[2], 2)} {radius}\n"
                ),
                range_rule,
            )
            # runs[f"{str(atom.get_id())}_range"] = (
            #     (
            #         f"RADIUS {round(xyz[0], 2)} {round(xyz[1], 2)} "
            #         f"{round(xyz[2], 2)} {radius}\n"
            #     ),
            #     range_rule,
            # )
            xyz_string.append(
                f"N {round(xyz[0], 2)} {round(xyz[1], 2)} "
                f"{round(xyz[2], 2)}"
            )

    with open("with_points.xyz", "w") as f:
        f.write(f"{len(xyz_string)}\n\n")
        f.write("\n".join(xyz_string))

    for run_name in runs:
        print(f"running {name}, {run_name}")
        oname = f"{name}_{run_name}"
        if not os.path.exists(f"{oname}.COMPLETE"):
            radius_string = runs[run_name][0]
            range_string = runs[run_name][1]
            input_xyz = "input_xyz.xyz"
            structure.write(input_xyz)
            input_string = (
                "1\n"
                "input_xyz.xyz\n"
                "FINE\n"
                f"{radius_string}"
                "INTEGRATE\n"
                f"{range_string}"
                f"ONAME {oname}\n"
            )
            input_file = f"{oname}_inp.in"
            with open(input_file, "w") as f:
                f.write(input_string)
            output_file = f"{oname}_out.out"
            cmd = f"{nci_path()} {input_file}"

            with open(output_file, "w") as f:
                # Note that sp.call will hold the program until completion
                # of the calculation.
                sp.call(
                    cmd,
                    stdin=sp.PIPE,
                    stdout=f,
                    stderr=sp.PIPE,
                    # Shell is required to run complex arguments.
                    shell=True,
                )
            # Write a done file, to avoid reruns.
            with open(f"{oname}.COMPLETE", "w") as f:
                f.write("yeah it is done!")
            print("done")

        print("plotting")
        plot_g(oname)
        print("done")

    os.chdir(init_dir)


def main():

    structure_dir = pathlib.Path("./structures_alpb")
    if not structure_dir.exists():
        os.mkdir(structure_dir)
    calculation_dir = pathlib.Path("./calculations_alpb")
    if not calculation_dir.exists():
        os.mkdir(calculation_dir)
    xtal_dir = pathlib.Path("./xtals")
    if not xtal_dir.exists():
        os.mkdir(xtal_dir)

    xtal_structures = (
        "faux_PM1_solo",
        "PM1_41Cl_repeat_extraction_cl",
        "PM1_41Cl_repeat_extraction_cl_bf4",
        "PM1_41Cl_repeat_extraction_solo",
        "PM1_41Cl_repeat_extraction_surrounds",
        "PMI_43Cl_repeat_extraction_cl",
        "PMI_43Cl_repeat_extraction_solo",
        "PMI_43Cl_repeat_extraction_surrounds",
    )
    for xtal in xtal_structures:
        mol_file = xtal_dir / f"{xtal}.mol"
        loaded_mol = stk.BuildingBlock.init_from_file(str(mol_file))
        run_nci(
            structure=loaded_mol,
            name=xtal,
            calc_dir=calculation_dir / f"{xtal}_nci",
        )

    ligand_names = ("l1q", "l1p", "l2q", "tq", "tp")
    c_orientations = ("A", "B", "C", "D")

    for i in ligand_names:
        for c_a in c_orientations:
            name_ = f"{i}_{c_a}"
            host_name = f"{i}_{c_a}_cl"
            cage_file = structure_dir / f"{name_}_optc.mol"
            hg_file = structure_dir / f"{host_name}_optc.mol"

            cage = stk.BuildingBlock.init_from_file(str(cage_file))
            run_nci(
                structure=cage,
                name=name_,
                calc_dir=calculation_dir / f"{name_}_nci",
            )
            if "l1" in i:
                hg = stk.BuildingBlock.init_from_file(str(hg_file))
                run_nci(
                    structure=hg,
                    name=host_name,
                    calc_dir=calculation_dir / f"{host_name}_nci",
                )


if __name__ == "__main__":
    main()
