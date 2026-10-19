#!/usr/bin/env bash
# APOE genotype (e2/e3/e4) for every sample in the QC'd target, from rs429358 and rs7412.
# Usage: bsub < 01_extract_apoe_status.sh   (writes ${WORKDIR}/apoe/APOE_status.tsv)
#BSUB -P re_gecip_neurology
#BSUB -q short
#BSUB -J apoe01 -o logs/%J_apoe01.out -e logs/%J_apoe01.err
set -euo pipefail
source "${PIPELINE_CONFIG:?export PIPELINE_CONFIG=/path/to/config.sh}"

module load ${PLINK19_MODULE}

BFILE="${WORKDIR}/merged/case-HC-DEE-EPI.QC"
OUT_DIR="${WORKDIR}/apoe"
OUT_FILE="${OUT_DIR}/APOE_status.tsv"

mkdir -p "$OUT_DIR"
TMP_DIR=$(mktemp -d)
trap 'rm -rf "$TMP_DIR"' EXIT

# rs ID first, then hg38 position, then hg19 position
find_snp() {
    local rsid="$1" pos_hg38="$2" pos_hg19="$3"
    local name
    name=$(awk -v r="$rsid" '$2 == r {print $2; exit}' "${BFILE}.bim")
    [ -n "$name" ] && { echo "$name"; return 0; }
    name=$(awk -v p="$pos_hg38" '$1 == 19 && $4 == p {print $2; exit}' "${BFILE}.bim")
    [ -n "$name" ] && { echo "$name"; return 0; }
    name=$(awk -v p="$pos_hg19" '$1 == 19 && $4 == p {print $2; exit}' "${BFILE}.bim")
    [ -n "$name" ] && { echo "$name"; return 0; }
    return 1
}

RS429358_NAME=$(find_snp rs429358 44908684 45411941) || {
    echo "ERROR: rs429358 not found in ${BFILE}.bim (searched rs ID, hg38, hg19)"; exit 1; }
RS7412_NAME=$(find_snp rs7412 44908822 45412079) || {
    echo "ERROR: rs7412 not found in ${BFILE}.bim (searched rs ID, hg38, hg19)"; exit 1; }

echo "rs429358 -> $RS429358_NAME"
echo "rs7412   -> $RS7412_NAME"

# count rs429358-C (defines e4) and rs7412-T (defines e2)
printf "%s\n%s\n" "$RS429358_NAME" "$RS7412_NAME" > "${TMP_DIR}/apoe_snps.txt"
printf "%s\tC\n%s\tT\n" "$RS429358_NAME" "$RS7412_NAME" > "${TMP_DIR}/apoe_alleles.txt"

plink --noweb \
    --bfile "$BFILE" \
    --extract "${TMP_DIR}/apoe_snps.txt" \
    --recode A \
    --recode-allele "${TMP_DIR}/apoe_alleles.txt" \
    --out "${TMP_DIR}/apoe" 2>&1 | tail -20

if [ ! -f "${TMP_DIR}/apoe.raw" ]; then
    echo "ERROR: PLINK extraction failed, see ${OUT_DIR}/apoe_plink.log"
    cp "${TMP_DIR}/apoe.log" "${OUT_DIR}/apoe_plink.log" 2>/dev/null || true
    exit 1
fi
head -1 "${TMP_DIR}/apoe.raw" | tr ' ' '\n' | head -10

# The double heterozygote (C=1, T=1) is e2/e4 or e1/e3 and cannot be resolved
# from unphased calls, so it is flagged and excluded downstream.
python3 - <<PYEOF
import csv
from collections import Counter

raw = "${TMP_DIR}/apoe.raw"
out = "${OUT_FILE}"
rs4_id = "${RS429358_NAME}"
rs7_id = "${RS7412_NAME}"

with open(raw) as f:
    header = f.readline().rstrip("\n").split()
c_col = next(c for c in header if c.startswith(rs4_id + "_"))
t_col = next(c for c in header if c.startswith(rs7_id + "_"))
assert c_col.endswith("_C"), f"expected _C suffix on {c_col}"
assert t_col.endswith("_T"), f"expected _T suffix on {t_col}"

def geno(c, t):
    lookup = {
        (0,0): "e3/e3",
        (0,1): "e2/e3",
        (0,2): "e2/e2",
        (1,0): "e3/e4",
        (1,1): "excluded_double_het",
        (2,0): "e4/e4",
    }
    return lookup.get((c,t), f"rare(C={c},T={t})")

rows = []
with open(raw) as f:
    r = csv.DictReader(f, delimiter=" ")
    for row in r:
        fid, iid = row["FID"].strip(), row["IID"].strip()
        cv, tv = row[c_col].strip(), row[t_col].strip()
        if cv in ("NA","") or tv in ("NA",""):
            rows.append((fid, iid, "NA", "NA", "NA", "NA", "NA"))
            continue
        c, t = int(float(cv)), int(float(tv))
        g = geno(c, t)
        e4_car = "e4_carrier" if g in ("e3/e4","e4/e4") else ("NA" if g=="excluded_double_het" else "non_e4_carrier")
        e2_car = "e2_carrier" if g in ("e2/e3","e2/e2") else ("NA" if g=="excluded_double_het" else "non_e2_carrier")
        e4_dose = "NA" if g == "excluded_double_het" else str({"e3/e3":0,"e2/e3":0,"e2/e2":0,"e3/e4":1,"e4/e4":2}.get(g,"NA"))
        e2_dose = "NA" if g == "excluded_double_het" else str({"e3/e3":0,"e2/e3":1,"e2/e2":2,"e3/e4":0,"e4/e4":0}.get(g,"NA"))
        rows.append((fid, iid, g, e4_car, e2_car, e4_dose, e2_dose))

with open(out, "w") as f:
    f.write("FID\tIID\tAPOE_genotype\tAPOE_e4_carrier\tAPOE_e2_carrier\tAPOE_e4_dose\tAPOE_e2_dose\n")
    for r in sorted(rows):
        f.write("\t".join(r) + "\n")

print(f"\nTotal samples: {len(rows)}")
print("Genotype distribution:")
for g, n in sorted(Counter(r[2] for r in rows).items()):
    print(f"  {g}: {n} ({100*n/len(rows):.1f}%)")
PYEOF

echo "Wrote ${OUT_FILE}"
