##conda activate ldsc

###run munge_sumstats.py
python ~/software/ldsc/munge_sumstats.py \
  --sumstats ../genomicSEM/CVaging_commonfactor_multivariatesGWAS_res.txt.gz \
  --merge-alleles ~/data/ldsc/w_hm3.snplist \
  --N 1783422 \
  --keep-pval \
  --p P \
  --out CVaging_commonfactorGWAS_munge


### ldsc
python ~/software/ldsc/ldsc.py \
  --h2 CVaging_commonfactorGWAS_munge.sumstats.gz \
  --ref-ld-chr ~/data/ldsc/eur_w_ld_chr/ \
  --w-ld-chr ~/data/ldsc/eur_w_ld_chr/ \
  --out CVaging_commonfactorGWAS_h2 





