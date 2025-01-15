#!/bin/bash
#SBATCH --partition=molscieng
#SBATCH --nodes=1
#SBATCH --ntasks=126
#SBATCH --output=run_%J_stdout.txt
#SBATCH --error=run_%J_stderr.txt
#SBATCH --time=0
#SBATCH --job-name={{ cookiecutter.friendly_top_name }}-%
#SBATCH --chdir=~/aot/cg/{{ cookiecutter.project_slug }}/%

module load GROMACS/2023.1-foss-2022a-CUDA-11.7.0

gmx grompp -p system.top -c neutral.gro -f martini_em.mdp  -o 1-min.tpr -po 1-min.mdp -maxwarn 2 
mpirun gmx_mpi mdrun -cpi -v -deffnm 1-min >> mdrun.log 2>&1

gmx grompp -p system.top -c 1-min.gro   -f martini_eq.mdp  -o 2-eq.tpr  -po 2-eq.mdp -maxwarn 1
mpirun gmx_mpi mdrun -cpi -v -deffnm 2-eq  >> mdrun.log 2>&1

gmx grompp -p system.top -c 2-eq.gro    -f martini_run.mdp -o 3-run.tpr -po 3-run.mdp -maxwarn 1
mpirun gmx_mpi mdrun -cpi -v -deffnm 3-run >> mdrun.log 2>&1
