#!/usr/bin/env python3
# Reformat each GWAS in the trait table to SNP A1 A2 BETA P [N] for munge_sumstats.
# A1 is the effect allele. Also writes the munge N arguments and the trait order.
#
# usage: 01_prep_sumstats.py traits.tsv gwas_dir out_dir ref_ld_prefix

import gzip
import math
import os
import sys


def open_any(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def is_number(x):
    try:
        float(x)
        return True
    except ValueError:
        return False


def read_traits(path):
    with open(path) as f:
        hdr = f.readline().rstrip("\n").split("\t")
        rows = []
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            rows.append(dict(zip(hdr, line.rstrip("\n").split("\t"))))
    return rows


def hm3_chrbp_map(ref_ld, cache):
    # LD score files are hg19, so chr:bp input has to be GRCh37 as well
    m = {}
    if os.path.exists(cache):
        with open(cache) as f:
            for line in f:
                k, v = line.rstrip().split("\t")
                m[k] = v
        return m
    for chrom in range(1, 23):
        with gzip.open(f"{ref_ld}.{chrom}.l2.ldscore.gz", "rt") as f:
            hdr = f.readline().split()
            ic, isnp, ibp = hdr.index("CHR"), hdr.index("SNP"), hdr.index("BP")
            for line in f:
                x = line.split()
                m[f"{x[ic]}:{x[ibp]}"] = x[isnp]
    with open(cache, "w") as f:
        for k, v in m.items():
            f.write(f"{k}\t{v}\n")
    return m


def chrbp_key(s):
    parts = s.split(":")
    if len(parts) < 2:
        return None
    return parts[0].replace("chr", "") + ":" + parts[1]


def prep_one(t, gwas_dir, out_dir, chrbp):
    src = os.path.join(gwas_dir, t["file"])
    dst = os.path.join(out_dir, t["name"] + ".txt")

    # n is a column name, a fixed N, or cases/controls (e.g. 29944/52538)
    n = t["n"]
    n_from_col = not is_number(n) and "/" not in n
    if n_from_col:
        args = "--N-col N"
    elif "/" in n:
        cas, con = n.split("/")
        args = f"--N-cas {cas} --N-con {con}"
    else:
        args = f"--N {n}"

    use_log = t["effect_type"].lower() == "or"
    total = kept = unmapped = 0
    with open_any(src) as fin, open(dst, "w") as fout:
        first = fin.readline()
        tab = "\t" in first  # METAL .tbl files are space padded

        def split(s):
            return s.rstrip("\n").split("\t") if tab else s.split()

        hdr = split(first)
        ix = {c: i for i, c in enumerate(hdr)}
        cols = [ix[t["snp_col"]], ix[t["a1_col"]], ix[t["a2_col"]], ix[t["effect_col"]], ix[t["p_col"]]]
        if n_from_col:
            cols.append(ix[n])
            fout.write("SNP\tA1\tA2\tBETA\tP\tN\n")
        else:
            fout.write("SNP\tA1\tA2\tBETA\tP\n")

        for line in fin:
            total += 1
            f = split(line)
            try:
                vals = [f[i] for i in cols]
            except IndexError:
                continue
            snp, a1, a2, b, p = vals[:5]
            if b in ("NA", "") or p in ("NA", ""):
                continue
            if t["snp_type"] == "chrbp":
                key = chrbp_key(snp)
                snp = chrbp.get(key) if key else None
                if snp is None:
                    unmapped += 1
                    continue
            elif not snp.startswith("rs"):
                continue
            if use_log:
                if float(b) <= 0:
                    continue
                b = repr(math.log(float(b)))
            out = [snp, a1.upper(), a2.upper(), b, p] + vals[5:]
            fout.write("\t".join(out) + "\n")
            kept += 1

    msg = f"{t['name']}: kept {kept} of {total}"
    if t["snp_type"] == "chrbp":
        msg += f" ({unmapped} not in HapMap3)"
    print(msg, file=sys.stderr)
    return args


def main():
    if len(sys.argv) != 5:
        sys.exit("usage: 01_prep_sumstats.py traits.tsv gwas_dir out_dir ref_ld_prefix")
    traits_file, gwas_dir, out_dir, ref_ld = sys.argv[1:]
    prep_dir = os.path.join(out_dir, "prep")
    os.makedirs(prep_dir, exist_ok=True)

    traits = read_traits(traits_file)
    names = [t["name"] for t in traits]
    if len(set(names)) != len(names):
        sys.exit("duplicate trait names in " + traits_file)

    chrbp = {}
    if any(t["snp_type"] == "chrbp" for t in traits):
        chrbp = hm3_chrbp_map(ref_ld, os.path.join(prep_dir, "hm3_chrbp_rsid.tsv"))

    with open(os.path.join(prep_dir, "traits.list"), "w") as lst:
        for t in traits:
            args = prep_one(t, gwas_dir, prep_dir, chrbp)
            lst.write(f"{t['name']}\t{args}\n")


if __name__ == "__main__":
    main()
