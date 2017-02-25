
# Nalu Wind Energy Utilities

This repository contains various pre- and post-processing utilities to be used
with Nalu CFD solver. 

## Pre-processing utilities

The preprocessor `nalu_preprocess` can perform various tasks on an exodus input
mesh database. See `examples/nalu_preprocess.yaml` for an example input deck.
Currently, the following preprocessor options are available.

**3-D Utilities**

- `init_abl_fields` : Initialize ABL velocity and temperature profiles based on
  user-defined tables.
- `generate_planes` : Generate horizontal sampling planes at desired heights for
  averaging statistics or source terms to drive ABL profiles.

**2-D Utilities**

- `calc_ndtw2d` : Calculate the nearest distance to wall for a 2-D airfoil-like
  geometry for use with RANS wall models.

# Build instructions

  ```
  git clone <repo_url>
  cd nalu_utils
  mkdir build
  cd build
  cp ../doconfig.sh .
  # Edit doconfig.sh to set Nalu/Trilinos paths appropriately
  ./doconfig.sh
  make
  ```

Once compiled, the executables are available in `build/preprocessing/` directory. For
example, `build/preprocessing/nalu_preprocess`. Example execution

```
mpiexec -np 1 preprocessing/nalu_preprocess -i nalu_preprocess.yaml
```
