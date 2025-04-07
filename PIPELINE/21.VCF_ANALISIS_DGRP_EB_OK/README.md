
### Create sample guides
#cat Sample_guides/CC_samps_VA.txt Sample_guides/AA_samps_VA.txt > Sample_guides/All_samps_VA.txt 
#cat Sample_guides/CC_samps_DGRP.txt Sample_guides/AA_samps_DGRP.txt > Sample_guides/All_samps_DGRP.txt 
#cat Sample_guides/CC_samps_ME.txt Sample_guides/AA_samps_ME.txt > Sample_guides/All_samps_ME.txt 
#cat Sample_guides/CC_samps_Netherlands.txt Sample_guides/AA_samps_Netherlands.txt > Sample_guides/All_samps_Netherlands.txt 
#cat Sample_guides/CC_samps_PA.txt Sample_guides/AA_samps_PA.txt > Sample_guides/All_samps_PA.txt 
#cat Sample_guides/CC_samps_Zimb.txt Sample_guides/AA_samps_Zimb.txt > Sample_guides/All_samps_Zimb.txt 

### calculate pi/D for all
for pop in VA ME Netherlands PA Zimb;
do
echo $pop

if [ "$pop" = "VA" ];
then
  VCF=top2R_1MB_VA.recode.vcf.gz
elif [ "$pop" = "DGRP" ]; 
then
  VCF=top2R_1MB_DGRP.recode.vcf.gz
elif [ "$pop" = "ME" ] || [ "$pop" = "Netherlands" ] || [ "$pop" = "PA" ] || [ "$pop" = "Zimb" ]; 
then
  VCF=top2R_1MB_World.recode.vcf.gz 
else "not considered"
fi

#run the code

for subset in  All CC AA; 
do 
 ./3.pi_D_fst_calculator.sh \
 pi_D \
 $VCF \
 30000 \
 15000 \
 ${pop}-${subset} \
 Sample_guides/${subset}_samps_${pop}.txt \


done
done

### calculate fst
for pop in VA ME Netherlands PA Zimb;
do
echo $pop

if [ "$pop" = "VA" ];
then
  VCF=top2R_1MB_VA.recode.vcf.gz
elif [ "$pop" = "DGRP" ]; 
then
  VCF=top2R_1MB_DGRP.recode.vcf.gz
elif [ "$pop" = "ME" ] || [ "$pop" = "Netherlands" ] || [ "$pop" = "PA" ] || [ "$pop" = "Zimb" ]; 
then
  VCF=top2R_1MB_World.recode.vcf.gz 
else "not considered"
fi

#run the code

 ./3.pi_D_fst_calculator.sh \
 fst \
 $VCF \
 30000 \
 15000 \
 ${pop}-${subset} \
 Sample_guides/AA_samps_${pop}.txt \
 Sample_guides/CC_samps_${pop}.txt 

done


### Zoom analysis
VCF=VA_World_merged.vcf.gz
for pop in VA ME Netherlands PA Zimb;
do
echo $pop

#run the code
for subset in  All CC AA; 
do 
 ./3.1.pi_D_fst_calculator.zoom.sh \
 pi_D \
 $VCF \
 1000 \
 500 \
 ${pop}-${subset} \
 Sample_guides/${subset}_samps_${pop}.txt \

done
done

### Calculate fst
VCF=VA_World_merged.vcf.gz

for pop2 in VA ME Netherlands PA Zimb;
do
for pop1 in VA ME Netherlands PA Zimb;
do

echo $pop1 $pop2

if [[ "$pop1" = "$pop2" ]]; then
echo "same"
else
#run the code

 ./3.1.pi_D_fst_calculator.zoom.sh \
 fst \
 $VCF \
 1000 \
 500 \
 ${pop1}v${pop2}-vs \
 Sample_guides/All_samps_${pop1}.txt \
 Sample_guides/All_samps_${pop2}.txt 

fi

done
done
