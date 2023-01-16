#!/usr/bin/env bash
#
#SBATCH -J dockerMap # A single job name for the array
#SBATCH -c 55
#SBATCH -N 1 # on one node
#SBATCH -t 30:00:00 #<= this may depend on your resources
#SBATCH --mem 55G #<= this may depend on your resources
#SBATCH -o ./slurmOutput/RunDest.%A_%a.out # Standard output
#SBATCH -e ./slurmOutput/RunDest.%A_%a.err # Standard error
#SBATCH -p bluemoon

pwd
echo $SLURM_CPUS_PER_TASK 

#sbatch --array=1-$( cat jobs.pt1.csv | wc -l  ) \
#1.scatter.sh  \
#jobs.pt1.csv 

#sbatch --array=1-$( cat jobs.pt2.csv | wc -l  ) \
#1.scatter.sh  \
#jobs.pt2.csv 

job_file=$1
#SLURM_ARRAY_TASK_ID=1

jobid=$( cat $job_file | sed "${SLURM_ARRAY_TASK_ID}q;d" )
echo $jobid
job=$( echo $jobid | sed 's/,/_/g')
echo $job

### slice out details
    chr=$( echo $jobid | cut -f1 -d',' )
    start=$( echo $jobid | cut -f2 -d',' )
    stop=$( echo $jobid | cut -f3 -d',' )
	echo $chr $start $stop


### load packages
module load spack/spack-0.18.1
spack load bcftools@1.14
spack load r@4.2.1 r-sf
spack load python

module load parallel-20190222-gcc-7.3.0-fhj23tz
module load vcftools-0.1.14-gcc-7.3.0-g4ntkwp

##### params
  popSet="PoolSeq"
  method="PoolSNP"
  maf="001"
  mac=5
  version=test
  script_dir="/users/j/n/jnunez2/scratch/test_DEST/poolSNP/snpCalling/innerscripts"
  wd="/users/j/n/jnunez2/scratch/test_DEST/poolSNP/output"
  syncPath1orig="/users/j/n/jnunez2/scratch/test_DEST/pipeline_output/*/*masked.sync.gz"
  syncPath2orig="/users/j/n/jnunez2/scratch/test_DEST/pipeline_output/*/*masked.sync.gz"
#####

## Run params
#  popSet=${1}
#  method=${2}
#  maf=${3}
#  mac=${4}
#  version=${5}
#  jobid=${6}
#  job=$( echo $jobid | sed 's/_/,/g')
#  script_dir=${7}

## working & temp directory
  outdir="${wd}/sub_vcfs"
    if [ ! -d $outdir ]; then
        mkdir $outdir
    fi

## get list of SNYC files based on popSet & method
### full list
  #syncPath1orig="${9}/*/*masked.sync.gz"
  #syncPath2orig="${10}/*/*masked.sync.gz"

### target pops
  if [[ "${popSet}" == "PoolSeq" ]]; then
    syncPath1=""
    syncPath2=${syncPath2orig}
  elif [[ "${popSet}" == "all" ]]; then
    syncPath1=${syncPath1orig}
    syncPath2=${syncPath2orig}
  fi

echo $( ls ${syncPath1} ${syncPath2})

## get job
#   job=$( cat ${wd}/${jobs} | sed "${SLURM_ARRAY_TASK_ID}q;d" )
#   jobid=$( echo ${job} | sed 's/,/_/g' )
  echo $job

## set up RAM disk
  #[ ! -d /dev/shm/$USER/ ] && mkdir /dev/shm/$USER/
  #[ ! -d /dev/shm/$USER/${SLURM_JOB_ID} ] && mkdir /dev/shm/$USER/${SLURM_JOB_ID}
  
  tmpdir=${wd}/tmp_dir_${job}
    if [ ! -d $tmpdir ]; then
        mkdir $tmpdir
    fi

echo "Temp dir is $tmpdir"

## get sub section
  subsection () {
  
  	#######debug
  	#syncFile="/users/j/n/jnunez2/scratch/test_DEST/pipeline_output/ME_bo_09_fall.r2/ME_bo_09_fall.r2.masked.sync.gz"
  	#chrslice=${job}
  	#echo $syncFile
  	#echo $chrslice
  	#tmpdir=${tmpdir}
  	#echo $tmpdir
  	#####
   
    syncFile=${1}
    chrslice=${2}
    tmpdir=${3}

	echo $chrslice
	
    pop=$( echo ${syncFile} | rev | cut -f1 -d'/' | rev | sed 's/.masked.sync.gz//g' )
	echo $pop
    chr=$( echo $chrslice | cut -f1 -d'_' )
    start=$( echo $chrslice | cut -f2 -d'_' )
    stop=$( echo $chrslice | cut -f3 -d'_' )
	echo $chr $start $stop

    echo ${pop}_${chrslice}
	
	echo "tabix step now"
    tabix -b 2 -s 1 -e 2 \
    ${syncFile} \
    ${chr}:${start}-${stop} > ${tmpdir}/${pop}_${chrslice}

  }
  
  export -f subsection

  echo "subset"

  if [[ "${method}" == "SNAPE" ]]; then
    echo "SNAPE" ${method}
    parallel -j 1 subsection ::: $( ls ${syncPath1} ${syncPath2} | grep "SNAPE" | grep "monomorphic" ) ::: ${job} ::: ${tmpdir}
  elif [[ "${method}" == "PoolSNP" ]]; then
    echo "PoolSNP" ${method}
    parallel -j 1 subsection ::: $( ls ${syncPath1} ${syncPath2} | grep -v "SNAPE" ) ::: ${job} ::: ${tmpdir}
  fi

#ls 
### paste function
  echo "paste"
  Rscript --no-save --no-restore ${script_dir}/paste.R ${job} ${tmpdir} ${method}

### run through SNP calling
  echo "SNP calling"

  if [[ "${method}" == "SNAPE" ]]; then
    echo $method
    cat ${tmpdir}/allpops.${method}.sites | python ${script_dir}/PoolSNP/PoolSnp.py \
    --sync - \
    --min-cov 4 \
    --max-cov 0.95 \
    --miss-frac 0.5 \
    --min-count 0 \
    --min-freq 0 \
    --posterior-prob 0.9 \
    --SNAPE \
    --names $( cat ${tmpdir}/allpops.${method}.names |  tr '\n' ',' | sed 's/,$//g' )  > ${tmpdir}/${job}.${popSet}.${method}.${maf}.${mac}.${version}.vcf

  elif [[ "${method}"=="PoolSNP" ]]; then
    echo $method

    cat ${tmpdir}/allpops.${method}.sites | python ${script_dir}/PoolSNP/PoolSnp.py \
    --sync - \
    --min-cov 4 \
    --max-cov 0.95 \
    --min-count ${mac} \
    --min-freq 0.${maf} \
    --miss-frac 0.5 \
    --names $( cat ${tmpdir}/allpops.${method}.names |  tr '\n' ',' | sed 's/,$//g' )  > ${tmpdir}/${job}.${popSet}.${method}.${maf}.${mac}.${version}.vcf
  fi

### compress and clean up
  echo "compress and clean"

  cat ${tmpdir}/${job}.${popSet}.${method}.${maf}.${mac}.${version}.vcf | vcf-sort | bgzip -c > ${outdir}/${job}.${popSet}.${method}.${maf}.${mac}.${version}.vcf.gz
  tabix -p vcf ${outdir}/${job}.${popSet}.${method}.${maf}.${mac}.${version}.vcf.gz

  rm -fr ${tmpdir}

### done
  echo "done"
