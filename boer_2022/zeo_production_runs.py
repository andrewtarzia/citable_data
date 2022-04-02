import os
import matplotlib.pyplot as plt


def get_values(lines, file):
    values = {}
    line_of_int = lines[0]
    nl = line_of_int.replace(file, '').replace('@  ', '').rstrip()
    nl = [i for i in nl.split(' ') if i]
    values['volume'] = float(nl[1])
    values['density'] = float(nl[3])
    values['asa_m2g-1'] = float(nl[9])
    values['nasa_m2g-1'] = float(nl[-1])

    return values


def get_res_values(lines, file):
    values = {}
    line_of_int = lines[0]
    nl = line_of_int.replace(file, '').replace('@  ', '').rstrip()
    nl = [i for i in nl.split(' ') if i]
    print(nl)
    values['di'] = float(nl[0])
    values['df'] = float(nl[1])
    values['dif'] = float(nl[2])

    return values


def main():

    crystals = [
        'Alpha_solvent_removed.cif', 
        'Beta_solvent_removed_mod_p1_avo.cif',
        'alpha_OLD_waters_removed_mod.cif',
        'C_phenyl.cif', 
        'C_phenyl_desolvated.cif',
    ]
    zeo_path = '/home/atarzia/software/zeo++-0.3/network'
    probes = [1.0, 1.1, 1.2, 1.3, 1.325, 1.4, 1.5, 1.55, 1.6, 1.7, 1.82, 1.9]
    sampl = 20000

    results = {}
    for cryst in crystals:
        resoutput = cryst.replace('.cif', f'_p.res')
        cmd = (
            f'{zeo_path} -ha -res '
            f'{resoutput} {cryst}'
        )
        if not os.path.exists(resoutput):
            os.system(cmd)

        with open(resoutput, 'r') as f:
            lines = f.readlines()
        values = get_res_values(lines, resoutput)
        print(f'{cryst}: {values}')
    
        results[cryst] = {}
        for probe in probes:
            results[cryst][probe] = {}
            output = cryst.replace('.cif', f'_{probe}_{sampl}_p.out')
            cmd = (
                f'{zeo_path} -ha -sa {probe} {probe} {sampl} '
                f'{output} {cryst}'
            )
            if not os.path.exists(output):
                os.system(cmd)

            with open(output, 'r') as f:
                lines = f.readlines()

            values = get_values(lines, output)
            results[cryst][probe][sampl] = values

    for cryst in results:
        suffix = cryst.replace('.cif', '')
        fig, ax = plt.subplots(figsize=(8, 5))
        xs = []
        acc = []
        nacc = []
        totals = []
        for probe in results[cryst]:
            da = results[cryst][probe][sampl]
            asa = da['asa_m2g-1']
            nasa = da['nasa_m2g-1']
            print(
                f'{cryst}, {probe}, {sampl}: '
                f'ASA={asa}, NASA={nasa}'
            )
            xs.append(probe)
            acc.append(asa)
            nacc.append(nasa)
            totals.append(asa+nasa)

        ax.plot(
            xs,
            acc,
            c='gold',
            marker='o',
            # edgecolor='k',
            # s=120,
            lw=3,
            markersize=12,
            label='access.'
        )
        ax.plot(
            xs,
            nacc,
            c='skyblue',
            marker='X',
            # edgecolor='k',
            # s=120,
            lw=3,
            markersize=12,
            label='nonaccess.'
        )
        ax.plot(
            xs,
            totals,
            c='k',
            marker='P',
            # edgecolor='k',
            # s=120,
            lw=3,
            markersize=12,
            label='total'
        )
        ax.legend(fontsize=16)

        ax.tick_params(axis='both', which='major', labelsize=16)
        ax.set_xlabel('probe radius', fontsize=16)
        ax.set_ylabel('surface area [m$^2$ g$^{-1}$]', fontsize=16)
        ax.set_xlim(0.9, 2.0)
        ax.set_ylim(-1, max(totals)+50)

        fig.tight_layout()
        fig.savefig(f'SA_vs_probe_{suffix}.pdf', dpi=720, bbox_inches='tight')
        plt.close()


if __name__ == '__main__':
    main()
