#!/usr/bin/env python3
import sys
from pathlib import Path

def read_counts(path):
    result = {}
    with open(path) as handle:
        for line in handle:
            gene, value = line.rstrip("\n").split("\t")
            if not gene.startswith("__"):
                result[gene] = int(value)
    return result

paths = sorted(sys.argv[1:], key=lambda value: Path(value).name)
if not paths:
    sys.exit("usage: merge_counts.py SAMPLE.counts.txt [...]")
tables = [read_counts(path) for path in paths]
genes = list(tables[0])
for path, table in zip(paths, tables):
    if list(table) != genes:
        sys.exit(f"Gene IDs in {path} differ from the first count file")
samples = [Path(path).name.removesuffix(".counts.txt") for path in paths]
print("gene_id\t" + "\t".join(samples))
for gene in genes:
    print(gene + "\t" + "\t".join(str(table[gene]) for table in tables))
