find ./ -name "*.besd" | while read file
do
   tissue=`basename $file`
   tissue=${tissue%.*}
   smr \
      --bfile ${PWD}/1000G_EUR_Phase3_plink/1000G.EUR.QC.1-22 \
      --peqtl-smr 0.05 \
      --gwas-summary CVaging_commonfactor_multivariatesGWAS_res.ma \
      --beqtl-summary ${tissue}/${tissue} \
      --out ${PWD}/SMR_output/CVaging_commonfactorGWAS_SMR_${tissue} \
      --diff-freq 0.2 \
      --thread-num 10
done
