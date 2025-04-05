#! /bin/bash


module load vcftools

# Print help text and exit.
if [ "$1" = "-h" ] || [ "$1" = "--help" ]; then
  echo "script input_vcf samps window step pop"
  exit 1
fi

#### calculate Taj d and Pi
input_vcf=${1}
samps=${2}
window=${3}
step=${4}
pop=${5}

### Evaluate 
if [ $samps = "all" ]; then
echo "samps set as 'all'"

vcftools \
--gzvcf $input_vcf \
--window-pi $window \
--window-pi-step $step \
--out Pi_${pop}_W_${window}_S_${step}

vcftools \
--gzvcf $input_vcf \
--TajimaD $window \
--out D_${pop}_W_${window}_S_${step}

else
echo "samps set as 'file'"

vcftools \
--gzvcf $input_vcf \
--keep $samps \
--window-pi $window \
--window-pi-step $step \
--out Pi_${pop}_W_${window}_S_${step}

vcftools \
--gzvcf $input_vcf \
--keep $samps \
--TajimaD $window \
--out D_${pop}_W_${window}_S_${step}

fi

### collect files
mkdir results

### clan up steps
mv D_${pop}_W_${window}_S_${step}.Tajima.D results
rm D_${pop}_W_${window}_S_${step}.log

mv Pi_${pop}_W_${window}_S_${step}.windowed.pi results
rm Pi_${pop}_W_${window}_S_${step}.log
