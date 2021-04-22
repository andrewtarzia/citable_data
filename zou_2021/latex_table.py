import sys
import pandas as pd


def pairs():
    return {
        'fe4ll': 'fe4dd',
        'fe4dd': 'fe4ll',
        'fe6ll': 'fe6dd',
        'fe6dd': 'fe6ll',
        'fe7ll': 'fe7dd',
        'fe7dd': 'fe7ll',
        'zn1ll': 'zn1ld',
        'zn1ld': 'zn1ll',
        'zn2ld': 'zn2ll',
        'zn2ll': 'zn2ld',
    }


def main():
    csv_file = sys.argv[1]

    data = pd.read_csv(csv_file)
    print(data)

    rel_energies = []
    for idx, row in data.iterrows():
        s = str(row['structure'])
        pair = pairs()[s]
        print(s, pair)
        e = float(row['energy'])
        pair_e = float(data[data['structure'] == pair]['energy'])
        print(e, pair_e)
        r = e-min([e, pair_e])
        rel_energies.append(r)

    data['relenergy'] = rel_energies
    data['relenergy_kjmol'] = [i*2625.5 for i in rel_energies]

    print(data.to_latex(
        index=False,
        columns=['structure', 'relenergy_kjmol'],
    ))


if __name__ == '__main__':
    main()