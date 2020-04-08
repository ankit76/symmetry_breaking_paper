import os
import sys

def getEnergy(fname):
  fh = open(fname, 'r')
  for line in fh:
      pass
  return float(line.split()[1])

vmc = 'mpirun -np 24 /projects/anma2640/VMC/rVMC/bin/VMC'

#for d in [1.0, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.6, 3.0, 4.0, 5.0, 6.0]:
for d in [4.0]:
  print(f'd = {d}', flush = True)
  os.system(f'rm {d} -rf; mkdir {d}; cd {d}; python ../prepVMC.py {d} > pyscf.out; cp ghf.txt kghf.txt')
  for wave in ["agp", "ghf", "kghf", "pfaffian"]:
    os.system(f'''cd {d}; mkdir {wave}; cd {wave};
                  ln -s ../FCIDUMP;
                  ln -s ../{wave}.txt ./{"hf" if wave in ["ghf", "kghf"] else "pairMat"}.txt;
                  ln -s ../../{wave}.json;'''
                  + vmc + f' {wave}.json > {wave}.out ')
    ene = getEnergy(f'{d}/{wave}/{wave}.out')
    print(f'{wave}:  {ene}', flush = True)
