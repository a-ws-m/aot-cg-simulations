Bilayer simulations
===================

After using cookiecutter to make a copy of this template, copy in the topology
(.itp) and structure (.gro) files for the model you want to use. Next, run
`setup-bilayer.sh`. This script creates an initial configuration for a bilayer
simulation in a new subdirectory. The script has one positional argument, which
specifies how many molecules of AOT to include (e.g. `setup-bilayer.sh 400`).

Within the new subdirectory, the `initial-eq.sh` script will perform an initial
equilibration, which will hopefully self-assemble into a bilayer in the
*xy*-plane. If the bilayer is not in this plane, re-run the simulation until it
is. Afterwards, the `bilayer-eq.sh` script will perform longer simulations, from
which the bilayer properties can be determined.