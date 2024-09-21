#!/usr/bin/env bash
#
#SBATCH -J fst.ag # A single job name for the array
#SBATCH -c 1
#SBATCH -N 1 # on one node
#SBATCH -t 12:00:00 #<= this may depend on your resources
#SBATCH --mem 80G #<= this may depend on your resources
#SBATCH -p bluemoon
#SBATCH -o ./slurmOutputSims/sims.%A_%a.out # Standard output
#SBATCH --array=1-608%50

module load slim/4.2.2 

met=reco_wins.master.txt

chr=$( cat $met  | sed "${SLURM_ARRAY_TASK_ID}q;d" | awk '{ print $1 }' )
reco_bin=$( cat $met  | sed "${SLURM_ARRAY_TASK_ID}q;d" | awk '{ print $2 }' )
startw=$( cat $met  | sed "${SLURM_ARRAY_TASK_ID}q;d" | awk '{ print $3 }' )
endw=$( cat $met  | sed "${SLURM_ARRAY_TASK_ID}q;d" | awk '{ print $4 }' )
Lenght=$( cat $met  | sed "${SLURM_ARRAY_TASK_ID}q;d" | awk '{ print $5 }' )
rho=$( cat $met  | sed "${SLURM_ARRAY_TASK_ID}q;d" | awk '{ print $6 }' )

###
echo $chr  $reco_bin  $startw $endw $Lenght $rho

mkdir output/${chr}_${startw}_${endw}_drift

for i in {1..100}
do
echo $i
slim \
-d "META='/gpfs2/scratch/jcnunez/fst_brent/simulations/win_meta/${chr}_${startw}_${endw}.meta.txt'" \
-d "DATA='/gpfs2/scratch/jcnunez/fst_brent/simulations/win_data/${chr}_${startw}_${endw}.data.txt'" \
-d "ith=${i}" \
introgression_drift.6_2_2024.slim
done
