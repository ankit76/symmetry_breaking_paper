import os
import sys

def getEnergy(fname):
  fh = open(fname, 'r')
  for line in fh:
      pass
  return float(line.split()[1])

vmc = 'mpirun -np 24 /projects/anma2640/VMC/rVMC/bin/VMC'
gfmc = 'mpirun -np 24 /projects/anma2640/VMC/rVMC/bin/GFMC'

#for d in [1.6, 1.8, 2.5]:
for d in [1.6]:
  print(f'd = {d}', flush = True)
  print('wave     eneVMC       eneGFMC', flush = True)
  os.system(f'rm {d} -rf; mkdir {d}; cd {d}; python ../prepVMC.py {d} > pyscf.out;')
  for wave in ["rhf", "uhf", "ghf"]:
    os.system(f'''cd {d}; mkdir {wave}; cd {wave};
                  ln -s ../FCIDUMP; ln -s ../bestDet;
                  ln -s ../{wave}.txt ./{"hf" if wave in ["rhf","uhf","ghf"] else "pairMat"}.txt;
                  ln -s ../../{wave}.json;
                  ln -s ../../{wave}GFMC.json;'''
                  + vmc + f' {wave}.json > {wave}.out;'
                  + gfmc + f' {wave}GFMC.json > {wave}GFMC.out ')
    eneVMC = getEnergy(f'{d}/{wave}/{wave}.out')
    eneGFMC = getEnergy(f'{d}/{wave}/{wave}GFMC.out')
    print(f'{wave}   {eneVMC}, {eneGFMC}', flush = True)
  for wave in ["agp", "pfaffian"]:
    os.system(f'''cd {d}; mkdir {wave}; cp ghf/BestDeterminant.txt {wave}/; cd {wave};
                  ln -s ../FCIDUMP;
                  ln -s ../{wave}.txt ./{"hf" if wave in ["rhf","uhf","ghf"] else "pairMat"}.txt;
                  ln -s ../../{wave}.json;
                  ln -s ../../{wave}GFMC.json;'''
                  + vmc + f' {wave}.json > {wave}.out;'
                  + gfmc + f' {wave}GFMC.json > {wave}GFMC.out ')
    eneVMC = getEnergy(f'{d}/{wave}/{wave}.out')
    eneGFMC = getEnergy(f'{d}/{wave}/{wave}GFMC.out')
    print(f'{wave}   {eneVMC}, {eneGFMC}', flush = True)
