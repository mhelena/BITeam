#!/usr/bin/env python3
"""Two-sided HWE exact test per site (Wigginton, Cutler & Abecasis 2005).

stdin:  CHROM  POS  then one GT per sample  (bcftools query -f '%CHROM\\t%POS[\\t%GT]\\n')
stdout: CHROM  POS  for sites with p < cutoff
usage:  hwe_exact.py [cutoff]    (default 1e-5)

Same p-value as PLINK --hardy without midp. Only 0/1 genotypes are counted;
multi-allelic and non-diploid calls are ignored.
"""
import sys

THRESH = float(sys.argv[1]) if len(sys.argv) > 1 else 1e-5


def hwe_exact_p(obs_het, obs_hom1, obs_hom2):
    n = obs_het + obs_hom1 + obs_hom2
    if n == 0:
        return 1.0
    obs_homc = max(obs_hom1, obs_hom2)
    obs_homr = min(obs_hom1, obs_hom2)
    rare = 2 * obs_homr + obs_het
    if rare == 0 or rare == 2 * n:
        return 1.0

    het_probs = [0.0] * (rare + 1)
    mid = rare * (2 * n - rare) // (2 * n)
    if (mid % 2) != (rare % 2):
        mid += 1

    het_probs[mid] = 1.0
    sum_p = 1.0

    curr_hets = mid
    curr_homr = (rare - mid) // 2
    curr_homc = n - curr_hets - curr_homr
    while curr_hets >= 2:
        het_probs[curr_hets - 2] = (
            het_probs[curr_hets] * curr_hets * (curr_hets - 1.0)
            / (4.0 * (curr_homr + 1) * (curr_homc + 1))
        )
        sum_p += het_probs[curr_hets - 2]
        curr_homr += 1
        curr_homc += 1
        curr_hets -= 2

    curr_hets = mid
    curr_homr = (rare - mid) // 2
    curr_homc = n - curr_hets - curr_homr
    while curr_hets <= rare - 2:
        het_probs[curr_hets + 2] = (
            het_probs[curr_hets] * 4.0 * curr_homr * curr_homc
            / ((curr_hets + 2.0) * (curr_hets + 1.0))
        )
        sum_p += het_probs[curr_hets + 2]
        curr_homr -= 1
        curr_homc -= 1
        curr_hets += 2

    p_obs = het_probs[obs_het]
    p_value = 0.0
    for i in range(rare + 1):
        if het_probs[i] <= p_obs:
            p_value += het_probs[i]
    return min(p_value / sum_p, 1.0)


def count_gts(fields):
    nhet = nhom_ref = nhom_alt = 0
    for gt in fields:
        if not gt or gt[0] == '.':
            continue
        if '|' in gt:
            a, sep, b = gt.partition('|')
        elif '/' in gt:
            a, sep, b = gt.partition('/')
        else:
            continue
        if a == '.' or b == '.':
            continue
        if a == '0' and b == '0':
            nhom_ref += 1
        elif a == '1' and b == '1':
            nhom_alt += 1
        elif (a == '0' and b == '1') or (a == '1' and b == '0'):
            nhet += 1
    return nhet, nhom_ref, nhom_alt


n_in = 0
n_fail = 0
for line in sys.stdin:
    parts = line.rstrip('\n').split('\t')
    if len(parts) < 3:
        continue
    n_in += 1
    chrom, pos = parts[0], parts[1]
    nhet, nhom_ref, nhom_alt = count_gts(parts[2:])
    p = hwe_exact_p(nhet, nhom_ref, nhom_alt)
    if p < THRESH:
        sys.stdout.write(f"{chrom}\t{pos}\n")
        n_fail += 1
    if n_in % 500000 == 0:
        print(f"  HWE: processed {n_in:,}, fails {n_fail:,}", file=sys.stderr, flush=True)

print(f"HWE: total {n_in:,}, fails {n_fail:,} (p<{THRESH})", file=sys.stderr)
