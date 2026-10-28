#!/usr/bin/env bash
set -euo pipefail
cd ~/adesola_1/reanalysis

# 1. Add headers to the raw BLAST file (no header present in NCBI export)
{ echo "qseqid,sseqid,pident,length,mismatch,gapopen,qstart,qend,sstart,send,evalue,bitscore"; \
  cat DHJ5C5HR014-Alignment-HitTable.csv; } > blast/lacto_blast_v2.csv

# 2. Keep top 10 hits per ASV (drop the tail of ~80 hits/ASV)
awk -F',' 'NR==1{print; next} {if (count[$1] < 10) {print; count[$1]++}}' \
  blast/lacto_blast_v2.csv > blast/lacto_top10.csv

echo "=== summary ==="
echo "Total hits in raw file: $(wc -l < blast/lacto_blast_v2.csv)"
echo "Kept top-10 hits:       $(wc -l < blast/lacto_top10.csv)"
echo "Unique ASVs:            $(tail -n +2 blast/lacto_top10.csv | cut -d',' -f1 | sort -u | wc -l)"
echo "Unique subject acc.:    $(tail -n +2 blast/lacto_top10.csv | cut -d',' -f2 | sort -u | wc -l)"
