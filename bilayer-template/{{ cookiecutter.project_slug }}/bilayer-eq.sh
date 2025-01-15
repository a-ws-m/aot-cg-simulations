#!/bin/bash
#SBATCH --partition=molscieng
#SBATCH --nodes=1
#SBATCH --ntasks=126
#SBATCH --output=bilayer_eq_%J_stdout.txt
#SBATCH --error=bilayer_eq_%J_stderr.txt
#SBATCH --time=0
#SBATCH --job-name={{ cookiecutter.friendly_top_name }}BilayerEq
#SBATCH --mail-type=ALL
#SBATCH --chdir=~/aot/cg/bilayer/{{ cookiecutter.project_slug }}

module load GROMACS/2023.1-foss-2022a-CUDA-11.7.0

gmx grompp -p system.top -c 3-run.gro -f martini_eq_semiiso.mdp -o 4-bilayer-eq.tpr -po 4-bilayer-eq.mdp -maxwarn 1
mpirun gmx_mpi mdrun -v -cpi -deffnm 4-bilayer-eq >> mdrun.log 2>&1

gmx grompp -p system.top -c 4-bilayer-eq.gro -f martini_run_semiiso.mdp -o 5-bilayer-run.tpr -po 5-bilayer-run.mdp -maxwarn 1
mpirun gmx_mpi mdrun -cpi -v -deffnm 5-bilayer-run >> mdrun.log 2>&1
