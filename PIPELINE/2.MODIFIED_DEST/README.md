# Modified DEST pipeline

Because of VACC dynamics. The flag "--home $1" was added to singularity. Also argument $1 must be a hard path

Other modifications have been made to "RunDocker" for optimization into the VACC cluster

## Run
```
wd=.
sbatch --array=1-20 \
--mail-type=ALL \
${wd}/runDocker_RUN.sh \
/users/j/n/jnunez2/scratch/test_DEST \
/netfiles02/lockwood_lab/IntrogressionProject/IntrogressionRawData/RawReads_Curated_2023 \
${wd}/pipeline_output \
/netfiles02/lockwood_lab/IntrogressionProject/IntrogressionRawData/mapping_meta_data_Jan5.2023.txt
```
