#!/bin/bash

magma \
   --annotate \
   --snp-loc CVaging_commonfactorGWAS_SNPloc.txt \
   --gene-loc NCBI37.3.gene.loc \
   --out CVaging_commonfactorGWAS_NCBI37.3
wait
magma \
   --bfile magma/g1000_eur/g1000_eur \
   --pval CVaging_commonfactorGWAS_pval.txt \
   N=1783422 \
   --gene-annot CVaging_commonfactorGWAS_NCBI37.3.genes.annot \
   --out CVaging_commonfactorGWAS_magma_gene_annot
