
# Misc. Nalu Wind Energy Utilities

- `gen_zplanes` : Generate sampling planes at given heights for ABL forcing algorithm
- `ndtw2d` : Brute-force wall distance calculation algorithm for 2-D airfoils


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
