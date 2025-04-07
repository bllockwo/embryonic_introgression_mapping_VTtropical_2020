#! /bin/bash

module load vcftools

# Print help text and exit.
if [ "$1" = "-h" ] || [ "$1" = "--help" ]; then
  echo "script test(pi_D, fst) input_vcf window step pop samps1 samps2(optional)"
  exit 1
fi

#### calculate Taj d and Pi
test1=${1}
input_vcf=${2}
window=${3}
step=${4}
pop=${5}
samps1=${6}
samps2=${7}

### Evaluate 
if [[ "$test1" = "pi_D" ]]; then
echo "evaluate pi/D"

vcftools \
--gzvcf $input_vcf \
--keep $samps1 \
--window-pi $window \
--window-pi-step $step \
--from-bp 20546633 \
--to-bp 20556633 \
--chr 2R \
--out Pi_${pop}_W_${window}_S_${step}

vcftools \
--gzvcf $input_vcf \
--keep $samps1 \
--TajimaD $window \
--from-bp 20546633 \
--to-bp 20556633 \
--chr 2R \
--out D_${pop}_W_${window}_S_${step}

### collect files
mkdir results_pi_D
### clwan up steps
mv D_${pop}_W_${window}_S_${step}.Tajima.D results
rm D_${pop}_W_${window}_S_${step}.log
mv Pi_${pop}_W_${window}_S_${step}.windowed.pi results
rm Pi_${pop}_W_${window}_S_${step}.log

elif [[ "$test1" = "fst" ]]; then
echo "evaluate fst"
vcftools \
--gzvcf $input_vcf \
--weir-fst-pop $samps1 \
--weir-fst-pop $samps2 \
--fst-window-size $window \
--fst-window-step  $step \
--from-bp 20546633 \
--to-bp 20556633 \
--chr 2R \
--out fst_${pop}_W_${window}_S_${step}

### collect files
mkdir results_fst
### clwan up steps
mv fst_${pop}_W_${window}_S_${step}.windowed.weir.fst results_fst
rm fst_${pop}_W_${window}_S_${step}.log
else 
echo "not considered"
fi


echo "done"