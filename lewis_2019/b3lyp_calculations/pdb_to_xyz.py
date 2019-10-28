from ase.io import read
from ase.io.xyz import write_xyz
import glob

for i in glob.glob('*.pdb'):
    print(i)
    a = read(i)
    cmt = 'This is an XYZ structure.'
    write_xyz(i.replace('.pdb', '.xyz'), images=a, comment=cmt,
              columns=['symbols', 'positions'])

