#!/bin/zsh
# Remove NA from aot-na.itp and aot-na.gro files in directory

# Remove the line with NA from aot-na.itp
sed -i '/NA/d' aot-na.itp
# And the same for aot-na.gro
sed -i '/NA/d' aot-na.gro
# And reduce the atom count in the second line by 1
original_num_atoms=$(sed '2q;d' aot-na.gro)
new_num_atoms=$(( $original_num_atoms - 1 ))
sed -i "2s/$original_num_atoms/$new_num_atoms/" aot-na.gro

# Finally, rename both files
mv aot-na.itp {{cookiecutter.topology_name}}.itp
mv aot-na.gro {{cookiecutter.topology_name}}.gro