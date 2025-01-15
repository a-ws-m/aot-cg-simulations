#!/bin/zsh

set -e

solute="{{ cookiecutter.topology_name }}.gro"
solute_top="{{ cookiecutter.topology_name }}.itp"
nsol=330
solvent_box="water.gro"
solute_name="AOT"
counterion="na.gro"
counterion_name="NA"
solvent_name="W"
solvent_atoms="1"
box_dim=12.0
box_sf=1

while [[ $# -gt 0 ]]; do
    case "$1" in
        -nsol)
            nsol="$2"
            shift 2
            ;;
        -box_sf)
            box_sf="$2"
            shift 2
            ;;
        *)
            echo "Invalid argument: $1"
            exit 1
            ;;
    esac
done

box_dim=$(($box_dim * $box_sf))

atoms_per_solute=$(sed '2q;d' $solute)

many_mol_gro="many-${solute}"

module load GROMACS/2023.1-foss-2022a-CUDA-11.7.0

# Generate box with several AOT molecules
gmx insert-molecules -ci ${solute} -nmol $nsol -radius 0.21 -try 100 -box $box_dim $box_dim $box_dim -o ${many_mol_gro}
total_solute_atoms=$(sed '2q;d' ${many_mol_gro})
solute_molecules=$(( $total_solute_atoms / $atoms_per_solute))

# Solvate the molecules in water
gmx solvate -cp ${many_mol_gro} -cs ${solvent_box} -radius 0.21 -o initial.gro

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
printf "AOT: %d\nWater: %d\nWeight percent: %.2f%%\n" $solute_molecules $solvent_molecules $wt_percent

generated_files=(initial.gro neutral.gro ${many_mol_gro} system.top)
files_to_link=(martini_em.mdp martini_eq.mdp martini_run.mdp martini.ff $solute $solute_top)

# Get user input
vared -p "Are you happy with the percentage? [y/n] " -c response

if [[ "$response" == "y" ]]; then
    # Copy all generated files to a new directory
    new_directory=$(printf "%.2f" $wt_percent | tr '.' '_')-percent
    mkdir $new_directory
    for file in $generated_files; do
        cp $file $new_directory
    done
    cp run-template.sh $new_directory/run.sh
    sed -i "s/%/$new_directory/g" $new_directory/run.sh
    cd $new_directory
    for file in $files_to_link; do
        ln -s ../$file
    done
else
    # Delete all generated files
    for file in $generated_files; do
        rm $file
    done
fi
