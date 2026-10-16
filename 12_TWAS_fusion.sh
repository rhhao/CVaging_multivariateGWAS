tissue=$1
ls ${PWAS}/TWAS/weights/*.pos|grep -v "nofilter"|while read weightfile
do
   for chr in `seq 1 22`
   do
      Rscript FUSION.assoc_test.R \
         --sumstats CVaging_commonfactor_multivariatesGWAS_twas_input.gz \
         --weights $weightfile \
         --weights_dir ${PWAS}/TWAS/weights/ \
         --ref_ld_chr ${PWD}/LDREF/1000G.EUR. \
         --chr $chr \
         --out ${PWD}/TWAS_output/CVaging_commonfactorGWAS_TWAS_${tissue}.${chr}.dat \
         --PANELN ${PWD}/TWAS/panel_sample_size.txt \
         --GWASN 1783422
   done
done
