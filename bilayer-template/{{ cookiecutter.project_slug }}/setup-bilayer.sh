#!/bin/zsh

solute="{{ cookiecutter.topology_name }}.gro"
nsol=330
solvent_box="water.gro"
solute_name="AOT"
counterion="na.gro"
counterion_name="NA"
solvent_name="W"
solvent_atoms="1"

if [ $# -eq 1 ]; then
    nsol="$1" # Set $nsol to the first command line argument
fi

atoms_per_solute=$(sed '2q;d' $solute)

many_mol_gro="many-${solute}"

module load GROMACS/2023.1-foss-2022a-CUDA-11.7.0

# Generate box with several AOT molecules
gmx insert-molecules -ci ${solute} -nmol $nsol -radius 0.21 -try 100 -box 10.0 10.0 10.0 -o ${many_mol_gro}
total_solute_atoms=$(sed '2q;d' ${many_mol_gro})
solute_molecules=$(( $total_solute_atoms / $atoms_per_solute))

# Solvate the molecules in water
gmx solvate -cp ${many_mol_gro} -cs ${solvent_box} -radius 0.21 -box 10.0 10.0 15.0 -o initial.gro

# Add counterions
gmx insert-molecules -ci ${counterion} -f initial.gro -nmol $solute_molecules -replace "resname W" -try 100 -radius 0.21 -o neutral.gro

# Set up topology file
cp system_empty.top system.top
solvent_lines=$(grep $solvent_name neutral.gro | wc -l)
solvent_molecules=$(expr $solvent_lines / $solvent_atoms )
echo "$solute_name               $solute_molecules" >> system.top
echo "$solvent_name               $solvent_molecules" >> system.top
echo "$counterion_name               $solute_molecules" >> system.top

# Calculate weight percent
aot_mass=$((444.6 * $solute_molecules))
water_mass=$((4 * 18 * $solvent_molecules))
wt_percent=$(( 100 * $aot_mass / ($water_mass + $aot_mass) ))
printf "AOT: %d\nWater: %d\nWeight percent: %.2f%%" $solute_molecules $solvent_molecules $wt_percent
