#!/usr/bin/env bash
# 18_build_supplementary.sh
# Builds all supplementary files in ~/adesola_1/submission/Supplementary/
# with correct names, formatted for submission.
set -euo pipefail

ROOT="$HOME/adesola_1"
CSV="$ROOT/reanalysis/outputs/csv"
TAB="$ROOT/reanalysis/outputs/tables"
SUPP="$ROOT/submission/Supplementary"

echo "=============================================="
echo "Building supplementary files"
echo "  Source CSVs : $CSV"
echo "  Source tables: $TAB"
echo "  Output      : $SUPP"
echo "=============================================="

mkdir -p "$SUPP"

# -----------------------------------------------------------------------------
# Check sources exist
# -----------------------------------------------------------------------------
echo
echo "[0] Checking source files..."
MISSING=0
for f in \
  "$CSV/read_tracking.csv" \
  "$CSV/read_tracking_postchim.csv" \
  "$CSV/audit_genus_per_sample.csv" \
  "$CSV/table_lactobacillus_species_calls.csv" \
  "$TAB/SupplementaryTableS3_species_calls.csv" \
  "$TAB/SupplementaryTableS4_species_contribution.csv"
do
  if [[ -f "$f" ]]; then
    echo "  OK   $(basename "$f")"
  else
    echo "  MISS $f"
    MISSING=$((MISSING + 1))
  fi
done
if [[ $MISSING -gt 0 ]]; then
  echo
  echo "ERROR: $MISSING source file(s) missing. Run 14_audit_clean.R and 15_figures_tables.R first."
  exit 1
fi

# -----------------------------------------------------------------------------
# S1 - Read tracking (combined from two source CSVs)
# -----------------------------------------------------------------------------
echo
echo "[S1] Read tracking (combined)..."
Rscript -e '
suppressPackageStartupMessages({})
CSV  <- file.path(Sys.getenv("HOME"), "adesola_1/reanalysis/outputs/csv")
OUT  <- file.path(Sys.getenv("HOME"), "adesola_1/submission/Supplementary",
                  "Supplementary_Table_S1_read_tracking.csv")
tr   <- read.csv(file.path(CSV, "read_tracking.csv"), stringsAsFactors = FALSE)
tr2  <- read.csv(file.path(CSV, "read_tracking_postchim.csv"), stringsAsFactors = FALSE)

out  <- data.frame(
  Sample_ID                       = tr$sample,
  Raw_reads                       = tr$raw_reads,
  After_primer_removal            = tr$after_primers,
  After_quality_filter            = tr$after_filter,
  Denoised_reads                  = tr$denoised,
  Non_chimeric_reads              = tr$nonchim,
  Percent_retained_from_raw       = tr$pct_final,
  Percent_non_chimeric_of_denoised = tr2$pct_of_denoised,
  stringsAsFactors = FALSE
)
write.csv(out, OUT, row.names = FALSE)
cat("    wrote", nrow(out), "rows\n")
'

# -----------------------------------------------------------------------------
# S2 - Genus-level relative abundance per sample
# -----------------------------------------------------------------------------
echo
echo "[S2] Genus-level relative abundance..."
cp "$CSV/audit_genus_per_sample.csv" \
   "$SUPP/Supplementary_Table_S2_genus_abundance.csv"
echo "    wrote $(($(wc -l < "$SUPP/Supplementary_Table_S2_genus_abundance.csv") - 1)) rows"

# -----------------------------------------------------------------------------
# S3 - Lactobacillus species assignments per ASV
# -----------------------------------------------------------------------------
echo
echo "[S3] Species assignments per ASV..."
cp "$TAB/SupplementaryTableS3_species_calls.csv" \
   "$SUPP/Supplementary_Table_S3_species_calls.csv"
echo "    wrote $(($(wc -l < "$SUPP/Supplementary_Table_S3_species_calls.csv") - 1)) rows"

# -----------------------------------------------------------------------------
# S4 - Species contribution summary
# -----------------------------------------------------------------------------
echo
echo "[S4] Species contribution to total Lactobacillus..."
cp "$TAB/SupplementaryTableS4_species_contribution.csv" \
   "$SUPP/Supplementary_Table_S4_species_contribution.csv"
echo "    wrote $(($(wc -l < "$SUPP/Supplementary_Table_S4_species_contribution.csv") - 1)) rows"

# -----------------------------------------------------------------------------
# Supplementary Figure S1 - Heatmap of top 20 genera
# -----------------------------------------------------------------------------
echo
echo "[S1 Fig] Heatmap of top 20 genera..."
Rscript "$(dirname "$0")/07b_build_heatmap.R"

# -----------------------------------------------------------------------------
# README
# -----------------------------------------------------------------------------
echo
echo "[README] Writing..."
cat > "$SUPP/README.md" <<'READMEEOF'
# Supplementary Materials

**Manuscript:** Vaginal microbiome profiles in women with and without HIV in Lagos, Nigeria: a descriptive pilot study using full-length PacBio 16S rRNA sequencing

**Authors:** Olalekan A, Beckley B, Oluwabunmi O, Agomma F, Ademuyiwa O

**Corresponding author:** Adesola Olalekan, adesolaolalekan@unilag.edu.ng

**BioProject accession:** PRJEB124225

**Code repository:** https://github.com/bbeckley-hub/Vaginal-Microbiome-Study-Nigeria

---

## Contents

| File | Description |
|---|---|
| `Supplementary_Table_S1_read_tracking.csv` | Read counts at each pipeline step per sample (raw → primer removal → quality filter → denoising → chimera removal) |
| `Supplementary_Table_S2_genus_abundance.csv` | Genus-level relative abundance (%) per sample (full table, all genera) |
| `Supplementary_Table_S3_species_calls.csv` | Species-level assignments per *Lactobacillus* ASV, including BLAST top hit, % identity, query coverage, second-best hit, final call, and reason |
| `Supplementary_Table_S4_species_contribution.csv` | Per-species contribution to total *Lactobacillus* reads across all samples |
| `Supplementary_Figure_S1_heatmap_top20_genera.png` | Heatmap of top 20 genera across samples (hierarchical clustering) |
| `Supplementary_Figure_S1_heatmap_top20_genera.pdf` | Vector version of the same figure |

## Source files and traceability

Every supplementary file is generated from analysis scripts available at the code repository. The mapping between supplement and source CSV is:

| Supplement | Source file in repository |
|---|---|
| S1 | `reanalysis/outputs/csv/read_tracking.csv` + `read_tracking_postchim.csv` |
| S2 | `reanalysis/outputs/csv/audit_genus_per_sample.csv` |
| S3 | `reanalysis/outputs/tables/SupplementaryTableS3_species_calls.csv` |
| S4 | `reanalysis/outputs/tables/SupplementaryTableS4_species_contribution.csv` |
| S1 Figure | `reanalysis/outputs/csv/audit_genus_per_sample.csv` |

## Software versions

- R 4.3.3
- DADA2 1.30.0
- cutadapt 5.2
- SILVA v138.1 nr99
- pheatmap (for Supplementary Figure S1)

## Notes

- **S1:** Primer removal achieved ≥97% read retention across all samples. Chimera removal retained 30.8%–80.7% of denoised reads, with the lowest retention in 01PND (35.6%) and 02PND (30.8%).
- **S2:** 19 ASVs (6.4% of total, low abundance) could not be assigned to a genus at the SILVA v138.1 confidence threshold and are not included in this table.
- **S3:** Of 84 *Lactobacillus* ASVs, 75 (89%) had at least one significant BLASTn hit; 9 (11%) did not. Of those 75, 10 are reported as candidate novel lineages and require further characterisation (see Discussion).
- **S4:** *L. iners* accounts for 74.1% of *Lactobacillus* reads; *L. crispatus* for 10.6%; 7.8% have no significant named BLAST hit.
- **Figure S1:** Rows are genera, columns are samples, values are raw relative abundance (%). Colour scale runs from white (low) to dark red (high). Samples cluster into two groups corresponding to the community states described in the manuscript.

## File reproducibility

All supplementary files can be regenerated from scratch by running:

Provided the upstream outputs (`seqtab_final.rds`, `taxonomy_final.rds`, and the CSVs in `reanalysis/outputs/`) are present.
READMEEOF

# -----------------------------------------------------------------------------
# Manifest
# -----------------------------------------------------------------------------
echo
echo "[MANIFEST] Writing..."
{
  echo "# Supplementary files manifest"
  echo "# Generated: $(date -Iseconds)"
  echo
  printf '%-60s %10s  %s\n' "FILE" "BYTES" "SHA256"
  find "$SUPP" -type f -not -name "MANIFEST.txt" | sort | while read -r f; do
    printf '%-60s %10s  %s\n' \
      "$(basename "$f")" \
      "$(stat -c%s "$f")" \
      "$(sha256sum "$f" | cut -d' ' -f1)"
  done
} > "$SUPP/MANIFEST.txt"

# -----------------------------------------------------------------------------
# Summary
# -----------------------------------------------------------------------------
echo
echo "=============================================="
echo "DONE. Supplementary files in:"
echo "  $SUPP"
echo
echo "Contents:"
ls -la "$SUPP"
echo
echo "Total size: $(du -sh "$SUPP" | cut -f1)"
echo "=============================================="
