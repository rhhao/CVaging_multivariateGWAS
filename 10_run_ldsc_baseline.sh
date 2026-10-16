#conda activaate ldsc
python ldsc.py \
  --h2 CVaging_commonfactorGWAS_munge.sumstats.gz \
  --ref-ld-chr baseline_v1.2/baseline. \
  --out CVaging_commonfactorGWAS_baseline_enrich \
  --overlap-annot  \
  --frqfile-chr 1000G_Phase3_frq/1000G.EUR.QC. \
  --w-ld-chr weights_hm3_no_hla/weights. \
  --print-coefficients
