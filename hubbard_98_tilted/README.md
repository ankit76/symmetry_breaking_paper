# Tilted 98 site Hubbard model with periodic boundary conditions

To get values for different U values

```bash
python sweepU.py > sweepU.out
```

The path to the VMC and GFMC executables needs to be updated in the sweepU.py script.
As an example, results for U = 2 are in the folder 2. It may be necessary to increase the number of stochastic samples to get converged results in some cases.
