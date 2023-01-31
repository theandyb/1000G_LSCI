#!/bin/sh

# Generate per-subtype singleton files

awk -F, '{if($9 == "AT_CG")print(substr($8,1,21))}' chr*_annotated.csv > AT_CG.txt
awk -F, '{if($9 == "AT_GC")print(substr($8,1,21))}' chr*_annotated.csv > AT_GC.txt
awk -F, '{if($9 == "AT_TA")print(substr($8,1,21))}' chr*_annotated.csv > AT_TA.txt
awk -F, '{if($9 == "GC_AT")print(substr($8,1,21))}' chr*_annotated.csv > GC_AT.txt
awk -F, '{if($9 == "GC_TA")print(substr($8,1,21))}' chr*_annotated.csv > GC_TA.txt
awk -F, '{if($9 == "GC_CG")print(substr($8,1,21))}' chr*_annotated.csv > GC_CG.txt
awk -F, '{if($9 == "cpg_GC_AT")print(substr($8,1,21))}' chr*_annotated.csv > cpg_GC_AT.txt
awk -F, '{if($9 == "cpg_GC_TA")print(substr($8,1,21))}' chr*_annotated.csv > cpg_GC_TA.txt
awk -F, '{if($9 == "cpg_GC_CG")print(substr($8,1,21))}' chr*_annotated.csv> cpg_GC_CG.txt

# Generate per-subtype control files

awk -F, '{if($4 == "AT_CG")print(substr($8,1,21))}' chr*_at.csv | sed 's/"//g' > AT_CG.txt
awk -F, '{if($4 == "AT_GC")print(substr($8,1,21))}' chr*_at.csv | sed 's/"//g' > AT_GC.txt
awk -F, '{if($4 == "AT_TA")print(substr($8,1,21))}' chr*_at.csv | sed 's/"//g' > AT_TA.txt
awk -F, '{if($4 == "GC_AT")print(substr($8,1,21))}' chr*_gc.csv | sed 's/"//g' > GC_AT.txt
awk -F, '{if($4 == "GC_TA")print(substr($8,1,21))}' chr*_gc.csv | sed 's/"//g' > GC_TA.txt
awk -F, '{if($4 == "GC_CG")print(substr($8,1,21))}' chr*_gc.csv | sed 's/"//g' > GC_CG.txt
awk -F, '{if($4 == "cpg_GC_AT")print(substr($8,1,21))}' chr*_gc.csv | sed 's/"//g' > cpg_GC_AT.txt
awk -F, '{if($4 == "cpg_GC_TA")print(substr($8,1,21))}' chr*_gc.csv | sed 's/"//g' > cpg_GC_TA.txt
awk -F, '{if($4 == "cpg_GC_CG")print(substr($8,1,21))}' chr*_gc.csv | sed 's/"//g' > cpg_GC_CG.txt

## All gc versions
awk -F, '{if($4 == "GC_AT")print(substr($8,1,21))}' chr*_gc_all.csv | sed 's/"//g' > all_GC_AT.txt
awk -F, '{if($4 == "GC_TA")print(substr($8,1,21))}' chr*_gc_all.csv | sed 's/"//g' > all_GC_TA.txt
awk -F, '{if($4 == "GC_CG")print(substr($8,1,21))}' chr*_gc_all.csv | sed 's/"//g' > all_GC_CG.txt
awk -F, '{if($4 == "cpg_GC_AT")print(substr($8,1,21))}' chr*_gc_all.csv | sed 's/"//g' >> all_GC_AT.txt
awk -F, '{if($4 == "cpg_GC_TA")print(substr($8,1,21))}' chr*_gc_all.csv | sed 's/"//g' >> all_GC_TA.txt
awk -F, '{if($4 == "cpg_GC_CG")print(substr($8,1,21))}' chr*_gc_all.csv | sed 's/"//g' >> all_GC_CG.txt

# closest and furthest controls
awk -F, '{if($4 == "AT_CG")print(substr($8,1,21))}' chr*_at.csv.min > AT_CG_min.txt
awk -F, '{if($4 == "AT_GC")print(substr($8,1,21))}' chr*_at.csv.min > AT_GC_min.txt
awk -F, '{if($4 == "AT_TA")print(substr($8,1,21))}' chr*_at.csv.min > AT_TA_min.txt
awk -F, '{if($4 == "GC_AT")print(substr($8,1,21))}' chr*_gc.csv.min > GC_AT_min.txt
awk -F, '{if($4 == "GC_TA")print(substr($8,1,21))}' chr*_gc.csv.min > GC_TA_min.txt
awk -F, '{if($4 == "GC_CG")print(substr($8,1,21))}' chr*_gc.csv.min > GC_CG_min.txt
awk -F, '{if($4 == "cpg_GC_AT")print(substr($8,1,21))}' chr*_gc.csv.min > cpg_GC_AT_min.txt
awk -F, '{if($4 == "cpg_GC_TA")print(substr($8,1,21))}' chr*_gc.csv.min > cpg_GC_TA_min.txt
awk -F, '{if($4 == "cpg_GC_CG")print(substr($8,1,21))}' chr*_gc.csv.min > cpg_GC_CG_min.txt

awk -F, '{if($4 == "AT_CG")print(substr($8,1,21))}' chr*_at.csv.max > AT_CG_max.txt
awk -F, '{if($4 == "AT_GC")print(substr($8,1,21))}' chr*_at.csv.max > AT_GC_max.txt
awk -F, '{if($4 == "AT_TA")print(substr($8,1,21))}' chr*_at.csv.max > AT_TA_max.txt
awk -F, '{if($4 == "GC_AT")print(substr($8,1,21))}' chr*_gc.csv.max > GC_AT_max.txt
awk -F, '{if($4 == "GC_TA")print(substr($8,1,21))}' chr*_gc.csv.max > GC_TA_max.txt
awk -F, '{if($4 == "GC_CG")print(substr($8,1,21))}' chr*_gc.csv.max > GC_CG_max.txt
awk -F, '{if($4 == "cpg_GC_AT")print(substr($8,1,21))}' chr*_gc.csv.max > cpg_GC_AT_max.txt
awk -F, '{if($4 == "cpg_GC_TA")print(substr($8,1,21))}' chr*_gc.csv.max > cpg_GC_TA_max.txt
awk -F, '{if($4 == "cpg_GC_CG")print(substr($8,1,21))}' chr*_gc.csv.max > cpg_GC_CG_max.txt

# 6 Category Versions
# Singletons
awk -F, '{if($4 == "A>C")print($3)}' chr*_annotated.csv > A_C.txt
awk -F, '{if($4 == "A>G")print($3)}' chr*_annotated.csv > A_G.txt
awk -F, '{if($4 == "A>T")print($3)}' chr*_annotated.csv > A_T.txt
awk -F, '{if($4 == "T>A")print($3)}' chr*_annotated.csv > T_A.txt
awk -F, '{if($4 == "T>C")print($3)}' chr*_annotated.csv > T_C.txt
awk -F, '{if($4 == "T>G")print($3)}' chr*_annotated.csv > T_G.txt
awk -F, '{if($4 == "C>A" && substr($9, 1, 3) == "cpg")print($3)}' chr*_annotated.csv > cpg_C_A.txt
awk -F, '{if($4 == "C>G" && substr($9, 1, 3) == "cpg")print($3)}' chr*_annotated.csv > cpg_C_G.txt
awk -F, '{if($4 == "C>T" && substr($9, 1, 3) == "cpg")print($3)}' chr*_annotated.csv > cpg_C_T.txt
awk -F, '{if($4 == "C>A" && substr($9, 1, 3) != "cpg")print($3)}' chr*_annotated.csv > C_A.txt
awk -F, '{if($4 == "C>G" && substr($9, 1, 3) != "cpg")print($3)}' chr*_annotated.csv > C_G.txt
awk -F, '{if($4 == "C>T" && substr($9, 1, 3) != "cpg")print($3)}' chr*_annotated.csv > C_T.txt
awk -F, '{if($4 == "G>A" && substr($9, 1, 3) == "cpg")print($3)}' chr*_annotated.csv > cpg_G_A.txt
awk -F, '{if($4 == "G>C" && substr($9, 1, 3) == "cpg")print($3)}' chr*_annotated.csv > cpg_G_C.txt
awk -F, '{if($4 == "G>T" && substr($9, 1, 3) == "cpg")print($3)}' chr*_annotated.csv > cpg_G_T.txt
awk -F, '{if($4 == "G>A" && substr($9, 1, 3) != "cpg")print($3)}' chr*_annotated.csv > G_A.txt
awk -F, '{if($4 == "G>C" && substr($9, 1, 3) != "cpg")print($3)}' chr*_annotated.csv > G_C.txt
awk -F, '{if($4 == "G>T" && substr($9, 1, 3) != "cpg")print($3)}' chr*_annotated.csv > G_T.txt

# Controls
awk -F, '{if($4 == "AT_CG" && $5 == "A")print($3)}' chr*_at.csv > A_C.txt
awk -F, '{if($4 == "AT_GC" && $5 == "A")print($3)}' chr*_at.csv > A_G.txt
awk -F, '{if($4 == "AT_TA" && $5 == "A")print($3)}' chr*_at.csv > A_T.txt
awk -F, '{if($4 == "AT_CG" && $5 == "T")print($3)}' chr*_at.csv > T_G.txt
awk -F, '{if($4 == "AT_GC" && $5 == "T")print($3)}' chr*_at.csv > T_C.txt
awk -F, '{if($4 == "AT_TA" && $5 == "T")print($3)}' chr*_at.csv > T_A.txt
awk -F, '{if($4 == "GC_TA" && $5 == "C")print($3)}' chr*_gc.csv > C_A.txt
awk -F, '{if($4 == "GC_AT" && $5 == "C")print($3)}' chr*_gc.csv > C_T.txt
awk -F, '{if($4 == "GC_CG" && $5 == "C")print($3)}' chr*_gc.csv > C_G.txt
awk -F, '{if($4 == "GC_TA" && $5 == "G")print($3)}' chr*_gc.csv > G_T.txt
awk -F, '{if($4 == "GC_AT" && $5 == "G")print($3)}' chr*_gc.csv > G_A.txt
awk -F, '{if($4 == "GC_CG" && $5 == "G")print($3)}' chr*_gc.csv > G_C.txt
awk -F, '{if($4 == "cpg_GC_TA" && $5 == "C")print($3)}' chr*_gc.csv > cpg_C_A.txt
awk -F, '{if($4 == "cpg_GC_AT" && $5 == "C")print($3)}' chr*_gc.csv > cpg_C_T.txt
awk -F, '{if($4 == "cpg_GC_CG" && $5 == "C")print($3)}' chr*_gc.csv > cpg_C_G.txt
awk -F, '{if($4 == "cpg_GC_TA" && $5 == "G")print($3)}' chr*_gc.csv > cpg_G_T.txt
awk -F, '{if($4 == "cpg_GC_AT" && $5 == "G")print($3)}' chr*_gc.csv > cpg_G_A.txt
awk -F, '{if($4 == "cpg_GC_CG" && $5 == "G")print($3)}' chr*_gc.csv > cpg_G_C.txt
