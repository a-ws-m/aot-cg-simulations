Bilayer simulations
===================

After using cookiecutter to make a copy of this template, copy in the topology
(.itp) and structure (.gro) files for the model you want to use. Next, run
`setup.sh`. This script creates an initial configuration for a bilayer
simulation in a new subdirectory. The script has two optional arguments: `-nsol`
to set the number of AOT molecules, and `-box_sf` to scale the box size (default
has a side of 12 nm). A `run.sh` script will also be created in the new
subdirectory, which performs the simulations.