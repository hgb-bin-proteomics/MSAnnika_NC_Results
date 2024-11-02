#!/usr/bin/env python3

# 2024 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

# the main purpose of this script is to extract proteins with 3 or more high-confidence
# PSMs and write their accessions to a json file

import json
import pandas as pd

from typing import Set

def get_proteins(target: pd.DataFrame, decoy: pd.DataFrame) -> Set[str]:
    colname = "Protein Accessions"
    proteins = dict() # Dict[str, int]
    for df in [target, decoy]:
        for i, row in df.iterrows():
            protein_accessions = [str(protein).strip() for protein in str(row[colname]).split(";")]
            for pa in protein_accessions:
                if pa in proteins:
                    proteins[pa] += 1
                else:
                    proteins[pa] = 1
    proteins_filtered = set() # Set[str]
    for key in proteins.keys():
        if proteins[key] > 2:
            proteins_filtered.add(key)
    return proteins_filtered

def main(argv = None) -> None:
    s1_target = pd.read_excel("S1_all_ng_linear_Amanda_deisotope_decharge.xlsx")
    s1_decoy = pd.read_excel("S1_all_ng_linear_Amanda_deisotope_decharge_decoys.xlsx")
    s2_target = pd.read_excel("S2_all_ng_linear_Amanda_deisotope_decharge.xlsx")
    s2_decoy = pd.read_excel("S2_all_ng_linear_Amanda_deisotope_decharge_decoys.xlsx")
    s3_target = pd.read_excel("S3_all_ng_linear_Amanda_deisotope_decharge.xlsx")
    s3_decoy = pd.read_excel("S3_all_ng_linear_Amanda_deisotope_decharge_decoys.xlsx")
    s1_proteins = get_proteins(s1_target, s1_decoy)
    print(f"Found {len(s1_proteins)} in S1.")
    s2_proteins = get_proteins(s2_target, s2_decoy)
    print(f"Found {len(s2_proteins)} in S2.")
    s3_proteins = get_proteins(s3_target, s3_decoy)
    print(f"Found {len(s3_proteins)} in S3.")
    proteins = s1_proteins.union(s2_proteins, s3_proteins)
    print(f"Found {len(proteins)} in total.")
    with open("proteins.json", "w", encoding = "utf-8") as f:
        json.dump(list(proteins), f)
        f.close()

if __name__ == "__main__":

    m = main()
