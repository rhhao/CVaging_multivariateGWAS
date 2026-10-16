weightfile="${PWD}/PWAS_weight/INTERVAL.pos"
cohort="INTERVAL"
for chr in `seq 1 22`
do
   Rscript FUSION.assoc_test.R \
      --sumstats CVaging_commonfactor_multivariatesGWAS_pwas_input.gz \
      --weights $weightfile \
      --weights_dir ${PWD}/PWAS_weight/ \
      --ref_ld_chr ${PWD}/LDREF/1000G.EUR. \
      --chr $chr \
      --out ${PWD}/PWAS_output/CVaging_commonfactorGWAS_PWAS_${cohort}.${chr}.dat \
      --GWASN 1783422
done

