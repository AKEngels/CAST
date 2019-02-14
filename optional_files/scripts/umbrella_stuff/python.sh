#PBS -l nodes=1:ppn=4                                                           
#PBS -l mem=4020MB
#PBS -l walltime=999:00:00
#PBS -N TITEL

cd $PBS_O_WORKDIR
cp -pr * $TMPDIR
cd $TMPDIR

python umbrella2d_helper.py

cp -pr * $PBS_O_WORKDIR
