#!/usr/bin/env bash
# Drops monomorphic sites (AC=0) from the post-QC merged VCF and, if
# OUTPUT_COPY_DIR is set, copies the result there and checks the copy.
# Submitted by `burden_submit.sh acgt0`.
set -euo pipefail

: "${ACGT0_IN:?}"

module load bcftools/1.16
module load htslib/1.21

IN="$ACGT0_IN"
FINAL="${IN%.vcf.gz}_ACgt0.vcf.gz"
TMP="${FINAL}.partial"

test -s "$IN" || { echo "Missing $IN" >&2; exit 1; }
if [[ -n "${OUTPUT_COPY_DIR:-}" ]]; then
    test -d "$OUTPUT_COPY_DIR" || { echo "Copy folder missing: $OUTPUT_COPY_DIR" >&2; exit 1; }
fi

rm -f "$FINAL" "${FINAL}.tbi" "$TMP"

# PLINK-exported VCFs carry no INFO/AC and bcftools merge does not add it, so
# AC/AN/AF are filled from the genotypes before --min-ac.
echo "[$(date)] recomputing AC/AN, then filtering --min-ac 1"
bcftools +fill-tags "$IN" -Ou -- -t AC,AN,AF |
  bcftools view --threads 2 --min-ac 1 -Oz -o "$TMP"

if ! bcftools view -h "$TMP" >/dev/null 2>&1; then
  echo "Output looks corrupt; aborting" >&2; exit 2
fi
N_REC=$(bcftools view -H "$TMP" | wc -l)
N_SMP=$(bcftools query -l "$TMP" | wc -l)
echo "[$(date)] $N_REC records, $N_SMP samples"

if [[ "$N_REC" -eq 0 ]]; then
  echo "Output has 0 records, something is wrong upstream; aborting" >&2
  exit 6
fi

mv "$TMP" "$FINAL"
bcftools index --tbi --threads 2 "$FINAL"
echo "[$(date)] wrote $FINAL"

if [[ -z "${OUTPUT_COPY_DIR:-}" ]]; then
  exit 0
fi

COPY="${OUTPUT_COPY_DIR}/$(basename "$FINAL")"
echo "[$(date)] copying to $OUTPUT_COPY_DIR"
cp -f "$FINAL"        "$COPY"
cp -f "${FINAL}.tbi"  "${COPY}.tbi"

SRC_MD5=$(md5sum "$FINAL"       | awk '{print $1}')
DST_MD5=$(md5sum "$COPY"        | awk '{print $1}')
SRC_TBI=$(md5sum "${FINAL}.tbi" | awk '{print $1}')
DST_TBI=$(md5sum "${COPY}.tbi"  | awk '{print $1}')

if [[ "$SRC_MD5" != "$DST_MD5" ]] || [[ "$SRC_TBI" != "$DST_TBI" ]]; then
  echo "MD5 mismatch between source and copy" >&2
  echo "  vcf src=$SRC_MD5 dst=$DST_MD5"      >&2
  echo "  tbi src=$SRC_TBI dst=$DST_TBI"      >&2
  exit 3
fi

N_REC2=$(bcftools view -H "$COPY" | wc -l)
N_SMP2=$(bcftools query -l "$COPY" | wc -l)
[[ "$N_REC2" -eq "$N_REC" ]] || { echo "Record count mismatch in copy" >&2; exit 4; }
[[ "$N_SMP2" -eq "$N_SMP" ]] || { echo "Sample count mismatch in copy" >&2; exit 5; }

echo "[$(date)] copy verified: md5 OK, $N_REC2 records, $N_SMP2 samples"
echo "  source: $FINAL"
echo "  copy:   $COPY"
