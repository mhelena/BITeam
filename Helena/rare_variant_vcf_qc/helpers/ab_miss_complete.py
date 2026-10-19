#!/usr/bin/env python3
"""Per-site allele balance, missingness and completeness from a bcftools query stream.

stdin:  CHROM  POS  ID  then one "GT|DP|AD" field per sample
        (bcftools query -f '%CHROM\\t%POS\\t%ID[\\t%GT|%DP|%AD]\\n')
usage:  ab_miss_complete.py <stats.tsv> <fail_ids.txt>

A site fails if more than 5% of samples have DP=0, fewer than 50% have a
called genotype, or fewer than 25% of hets are balanced (binomial p >= 0.01
around 0.5).
"""
import math
import sys

out_stats = sys.argv[1]
out_fail  = sys.argv[2]

AB_RATIO_MIN  = 0.25
AB_BINOM_P    = 0.01
MISSING_MAX   = 0.05
COMPLETE_MIN  = 0.50
MIN_HET       = 1

# Smallest minor-allele count that still passes the two-sided binomial test,
# for every depth up to N_MAX. Deeper sites use the normal approximation.
N_MAX = 1000
thresh = [0] * (N_MAX + 1)
for N in range(1, N_MAX + 1):
    pow2 = 1 << N
    cutoff = (AB_BINOM_P * pow2) / 2.0
    c = 1
    cum = 1
    if cum >= cutoff:
        thresh[N] = 0
        continue
    found = False
    for m in range(1, N // 2 + 1):
        c = c * (N - m + 1) // m
        cum += c
        if cum >= cutoff:
            thresh[N] = m
            found = True
            break
    if not found:
        thresh[N] = N // 2 + 1

_z = 2.5758293035489004  # two-sided alpha = 0.01
def is_balanced(alt, total):
    if total <= 0:
        return False
    if total > N_MAX:
        delta = abs(alt / total - 0.5)
        return delta <= _z * 0.5 / math.sqrt(total)
    return min(alt, total - alt) >= thresh[total]

first_line = sys.stdin.readline()
stats_fh = open(out_stats, 'w')
stats_fh.write("CHROM\tPOS\tID\tn_called\tn_dp0\tn_het\tn_balanced\t"
               "miss_rate\tcomplete_rate\tab_ratio\tfilter\n")
fail_fh = open(out_fail, 'w')

if not first_line:
    # A contig can have no sites left after a sample subset (chrM in EUR-only
    # runs). Leave header-only outputs so the per-chromosome job still succeeds.
    stats_fh.close()
    fail_fh.close()
    print("empty input stream, wrote header-only stats and no fails", file=sys.stderr)
    sys.exit(0)

parts = first_line.rstrip('\n').split('\t')
n_samples = len(parts) - 3
print(f"Samples per row: {n_samples}", file=sys.stderr)

n_total      = 0
n_fail_miss  = 0
n_fail_comp  = 0
n_fail_ab    = 0

def process(parts):
    global n_total, n_fail_miss, n_fail_comp, n_fail_ab
    n_total += 1
    chrom, pos, vid = parts[0], parts[1], parts[2]

    n_called   = 0
    n_dp0      = 0
    n_het      = 0
    n_balanced = 0

    for s in parts[3:]:
        try:
            gt, dp_s, ad_s = s.split('|', 2)
        except ValueError:
            continue

        if dp_s in ('.', ''):
            dp = 0
        else:
            try:
                dp = int(dp_s)
            except ValueError:
                dp = 0
        if dp == 0:
            n_dp0 += 1

        if gt in ('./.', '.', '.|.', ''):
            continue
        if '|' in gt:
            gt = gt.replace('|', '/')
        n_called += 1

        if '/' in gt:
            a, b = gt.split('/', 1)
            if a != b and a != '.' and b != '.':
                if ',' in ad_s and '.' not in ad_s:
                    ad_parts = ad_s.split(',')
                    try:
                        n_ref = int(ad_parts[0])
                        n_alt = int(ad_parts[1])
                    except (ValueError, IndexError):
                        continue
                    total = n_ref + n_alt
                    if total > 0:
                        n_het += 1
                        if is_balanced(n_alt, total):
                            n_balanced += 1

    miss_rate     = n_dp0    / n_samples
    complete_rate = n_called / n_samples
    ab_ratio      = (n_balanced / n_het) if n_het >= MIN_HET else 1.0

    fr = []
    if miss_rate > MISSING_MAX:
        fr.append("miss"); n_fail_miss += 1
    if complete_rate < COMPLETE_MIN:
        fr.append("comp"); n_fail_comp += 1
    if n_het >= MIN_HET and ab_ratio < AB_RATIO_MIN:
        fr.append("ab"); n_fail_ab += 1

    filt = "PASS" if not fr else ",".join(fr)
    stats_fh.write(
        f"{chrom}\t{pos}\t{vid}\t{n_called}\t{n_dp0}\t{n_het}\t{n_balanced}\t"
        f"{miss_rate:.4f}\t{complete_rate:.4f}\t{ab_ratio:.4f}\t{filt}\n"
    )
    if fr:
        site_id = vid if vid != '.' else f"{chrom}:{pos}"
        fail_fh.write(site_id + '\n')

process(parts)
for line in sys.stdin:
    parts = line.rstrip('\n').split('\t')
    if len(parts) < 3:
        continue
    process(parts)
    if n_total % 500000 == 0:
        print(f"  processed {n_total:,} sites...", file=sys.stderr, flush=True)

stats_fh.close()
fail_fh.close()

print(f"Total processed:                  {n_total:,}", file=sys.stderr)
print(f"Fail missingness (>5% DP=0):      {n_fail_miss:,}", file=sys.stderr)
print(f"Fail completeness (<50% called):  {n_fail_comp:,}", file=sys.stderr)
print(f"Fail AB ratio (<0.25):            {n_fail_ab:,}", file=sys.stderr)
