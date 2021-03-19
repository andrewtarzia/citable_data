#PBS -N _s_3d1_a
#PBS -l select=1:ncpus=32:mem=124gb
#PBS -l walltime=72:00:00
module load gaussian/g16-c01-avx

cd $PBS_O_WORKDIR

g16 < s_3d1_a.gau > s_3d1_a.log
