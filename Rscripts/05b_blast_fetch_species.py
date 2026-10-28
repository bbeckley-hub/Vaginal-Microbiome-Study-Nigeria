#!/usr/bin/env python3
import csv, json, time, urllib.request, re

input_csv  = "blast/lacto_nt_final.csv"
output_tsv = "blast/accession_species_master.tsv"

accs = set()
with open(input_csv) as f:
    for row in csv.DictReader(f):
        acc = (row.get("sseqid") or "").strip()
        if acc:
            accs.add(acc)
accs = sorted(accs)
print(f"Fetching species for {len(accs)} unique accessions...")

species = {}
BATCH = 200
for i in range(0, len(accs), BATCH):
    batch = accs[i:i+BATCH]
    ids   = ",".join(batch)
    url   = (f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
             f"?db=nuccore&id={ids}&retmode=json")
    try:
        with urllib.request.urlopen(url, timeout=60) as r:
            data = json.load(r)
        n_resolved = 0
        for uid, info in data.get("result", {}).items():
            if uid == "uids":
                continue
            acc   = info.get("accessionversion", "")
            acc   = re.sub(r'\.\d+$', '', acc)   # strip version to match input
            title = info.get("title", "")
            words = title.split()
            if len(words) >= 2 and acc:
                species[acc] = f"{words[0]} {words[1]}"
                n_resolved += 1
        print(f"  batch {i//BATCH + 1}: {len(batch)} accessions, {n_resolved} resolved")
    except Exception as e:
        print(f"  batch {i//BATCH + 1} FAILED: {e}")
    time.sleep(0.4)

with open(output_tsv, "w") as f:
    for acc in accs:
        f.write(f"{acc}\t{species.get(acc, 'NA')}\n")

print(f"\nWrote {len(accs)} rows to {output_tsv}")
print(f"Resolved {len(species)} species names")
