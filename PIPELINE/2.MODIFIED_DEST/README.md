# Modified DEST pipeline

Because of VACC dynamics. The flag "--home $1" was added to singularity. Also argument $1 must be a hard path

Other modifications have been made to "RunDocker" for optimization into the VACC cluster

## Run
wd=.
sbatch --array=1-1 \
${wd}/DEST_freeze1/mappingPipeline/scripts/runDocker.sh \
/users/j/n/jnunez2/scratch/test_DEST \
/netfiles02/lockwood_lab/IntrogressionProject/IntrogressionRawData/MERGED_READS_fixed \
${wd}/pipeline_output \
/netfiles02/lockwood_lab/IntrogressionProject/IntrogressionRawData/sampleId_meta_data.txt

