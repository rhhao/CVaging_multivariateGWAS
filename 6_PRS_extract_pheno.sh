#!/bin/bash

INPUT_FILE = "ICD10_Date_2026_0307_long.txt"
OUTPUT_PHENO = "prsice2_phenotype_HF.txt"
FID_COL = 1
IID_COL = 1
ICD10_COL = 2
SEPARATOR = "\t"

ICD10_CODES = "I50.0|I50.1|I50.9|I11.0|I13.0|I13.2"


echo "Processing $INPUT_FILE ..."

awk -v FID="$FID_COL" \
    -v IID="$IID_COL" \
    -v ICD="$ICD10_COL" \
    -v codes="$ICD10_CODES" \
    -v sep="$SEPARATOR" \
    -v out="$OUTPUT_PHENO" \
'
BEGIN {
    FS = sep;
    OFS = "\t";
    split(codes, code_arr, "|")
}

{
    if (NF == 0 || $1 ~ /^#/ || NR == 1) next;

    fid_val = $FID;
    iid_val = $IID;
    icd_val = $ICD;

    if (FID == IID) {
        id = iid_val;
    } else {
        id = fid_val " " iid_val;
    }

    if (!(id in pheno)) {
        pheno[id] = 1;
        first_icd[id] = icd_val;
    }

    for (i in code_arr) {
        pattern = code_arr[i];
        if (icd_val ~ "^" pattern) {
            pheno[id] = 2;
            matched_icd[id] = icd_val;
            break;
        }
    }
}

END {
    print "FID\tIID\tPhenotype" > out;

    for (id in pheno) {
        split(id, id_parts, " ");
        if (length(id_parts) == 1) {
            fid_out = id;
            iid_out = id;
        } else {
            fid_out = id_parts[1];
            iid_out = id_parts[2];
        }
        
        # output format：FID IID Phenotype
        print fid_out, iid_out, pheno[id] > out;
    }
    

    print "Done！" > "/dev/stderr";
    printf "Phenotype file save to: %s\n", out > "/dev/stderr";
    
}
' "$INPUT_FILE"

