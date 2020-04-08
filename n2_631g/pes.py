import os
import sys

def getEnergy(fname):
  fh = open(fname, 'r')
  for line in fh:
      pass
  return float(line.split()[1])

vmc = 'mpirun -np 24 /projects/anma2640/VMC/rVMC/bin/VMC'

#for d in [1.8, 2.0, 2.2, 2.5, 3.0, 3.6, 4.2, 5.0]:
for d in [1.6]:
  print(f'd = {d}', flush = True)
  os.system(f'rm {d} -rf; mkdir {d}; cd {d}; python ../prepVMC.py {d} > pyscf.out;')
  for wave in ["rhf", "uhf", "ghf", "agp", "pfaffian"]:
    os.system(f'''cd {d}; mkdir {wave}; cd {wave};
                  ln -s ../FCIDUMP; ln -s ../bestDet;
                  ln -s ../{wave}.txt ./{"hf" if wave in ["rhf","uhf","ghf"] else "pairMat"}.txt;
                  ln -s ../../{wave}.json;'''
                  + vmc + f' {wave}.json > {wave}.out ')
    ene = getEnergy(f'{d}/{wave}/{wave}.out')
    print(f'{wave}:  {ene}', flush = True)
