import os
import sys

def getEnergy(fname):
  fh = open(fname, 'r')
  for line in fh:
      pass
  return float(line.split()[1])

vmc = 'mpirun -np 24 /projects/anma2640/VMC/rVMC/bin/VMC'
gfmc = 'mpirun -np 24 /projects/anma2640/VMC/rVMC/bin/GFMC'

#for U in [4, 8, 10, 20]:
for U in [4]:
  print(f'U = {U}', flush = True)
  os.system(f'rm {U} -rf; mkdir {U}; cd {U}; python ../prepVMC.py {U} > pyscf.out; cp ghf.txt kghf.txt; cp pfaffian.txt kpfaffian.txt')
  for wave in ["ghf", "kghf", "pfaffian", "kpfaffian"]:
    os.system(f'''cd {U}; mkdir {wave}; cd {wave};
                  ln -s ../FCIDUMP;
                  ln -s ../{wave}.txt ./{"hf" if wave in ["ghf", "kghf"] else "pairMat"}.txt;
                  ln -s ../../{wave}.json;
                  ln -s ../../{wave}GFMC.json;'''
                  + vmc + f' {wave}.json > {wave}.out;'
                  + gfmc + f' {wave}GFMC.json > {wave}GFMC.out;')
    eneVMC = getEnergy(f'{U}/{wave}/{wave}.out')
    eneGFMC = getEnergy(f'{U}/{wave}/{wave}GFMC.out')
    print(f'{wave}:  {eneVMC}, {eneGFMC}', flush = True)
