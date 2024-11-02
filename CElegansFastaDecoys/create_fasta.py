#!/usr/bin/env python3

# 2024 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

import json
from Bio import SeqIO

def main(argv = None) -> None:
    proteins = []
    with open("proteins.json", "r", encoding = "utf-8") as f:
        proteins = json.load(f)
        f.close()
    print(f"Loaded {len(proteins)} proteins from json (should 3266).")
    print(len(proteins) == 3266)
    fasta_proteins = list(SeqIO.parse("CElegans_all_2024-03-15/uniprotkb_proteome_UP000001940_2024_03_15.fasta", "fasta"))
    print(f"Loaded {len(fasta_proteins)} proteins from fasta (should be 26695).")
    print(len(fasta_proteins) == 26695)
    proteins_to_keep = list()
    for protein in fasta_proteins:
        for keep in proteins:
            if keep in protein.id:
                proteins_to_keep.append(protein)
                break
    print(f"Found {len(proteins_to_keep)} proteins to keep (should be 3266).")
    print(len(proteins_to_keep) == 3266)
    SeqIO.write(proteins_to_keep, "filtered_w_decoys_uniprotkb_proteome_UP000001940_2024_03_15.fasta", "fasta")
    print("Wrote fasta file: filtered_w_decoys_uniprotkb_proteome_UP000001940_2024_03_15.fasta")

if __name__ == "__main__":

    m = main()
