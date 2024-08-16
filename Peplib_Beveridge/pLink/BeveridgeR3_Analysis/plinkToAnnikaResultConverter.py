#!/usr/bin/env python3

# pLink 2.3.11 filtered_cross-linked_peptides Result file to MS Annika Result file converter
# 2024 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

import sys
import pandas as pd

__version = "1.0.0"
__date = "2024-06-20"

def get_crosslinker_type(protein_idA, protein_idB):
    if protein_idA.upper() == protein_idB.upper():
        return "Intra"
    else:
        return "Inter"

def get_proteins(crosslink):
    proteins = crosslink.split(",")[4]
    protein_1 = proteins.split("-")[0]
    protein_2 = proteins.split("-")[1]
    protein_1 = protein_1.split("(")[0].strip()
    protein_2 = protein_2.split("(")[0].strip()
    return [protein_1, protein_2]

def get_nr_proteins(xl_type):
    if xl_type == "Intra":
        return 1
    return 2

def get_sequences(crosslink):
    sequences = crosslink.split(",")[1]
    sequence_1 = sequences.split("-")[0]
    sequence_2 = sequences.split("-")[1]
    sequence_1 = sequence_1.split("(")[0].strip().upper()
    sequence_2 = sequence_2.split("(")[0].strip().upper()
    return [sequence_1, sequence_2]

def get_positions(crosslink):
    positions = crosslink.split(",")[1]
    position_1 = positions.split("-")[0]
    position_2 = positions.split("-")[1]
    position_1 = int(position_1.split("(")[1].split(")")[0])
    position_2 = int(position_2.split("(")[1].split(")")[0])
    return [position_1, position_2]

def get_score(csms):
    max_score = 0.0
    for csm in csms:
        score = float(csm.split(",")[6])
        if score > max_score:
            max_score = score
    return max_score

def create_annika_result(plink_filename, crosslinker = "DSS"):
    # load file
    plink_lines = []
    with open(plink_filename, "r", encoding = "utf-8") as f:
        plink_lines = f.readlines()
        f.close()

    if len(plink_lines) == 0:
        raise RuntimeError("Couldn't read file.")

    crosslinks = dict()
    current_xl = ""
    for line in plink_lines[2:]:
        if line.strip()[0] == ",":
            crosslinks[current_xl].append(line)
        else:
            current_xl = line
            crosslinks[current_xl] = list()
    print(f"read {len(crosslinks)} crosslinks from file.")
    nrows = len(crosslinks)

    # columns
    Checked = ["FALSE" for i in range(nrows)]
    Crosslinker = [crosslinker for i in range(nrows)]
    Crosslink_Type = [get_crosslinker_type(get_proteins(crosslink)[0], get_proteins(crosslink)[1]) for crosslink in crosslinks]
    CSMs = [len(crosslinks[crosslink]) for crosslink in crosslinks]
    Proteins = [get_nr_proteins(xl_type) for xl_type in Crosslink_Type]
    Sequence_A = [get_sequences(crosslink)[0] for crosslink in crosslinks]
    Accession_A = [get_proteins(crosslink)[0] for crosslink in crosslinks]
    Position_A = [get_positions(crosslink)[0] for crosslink in crosslinks]
    Sequence_B = [get_sequences(crosslink)[1] for crosslink in crosslinks]
    Accession_B = [get_proteins(crosslink)[1] for crosslink in crosslinks]
    Position_B = [get_positions(crosslink)[1] for crosslink in crosslinks]
    Protein_Descriptions_A = ["" for i in range(nrows)]
    Protein_Descriptions_B = ["" for i in range(nrows)]
    Best_CSM_Score = [get_score(crosslinks[crosslink]) for crosslink in crosslinks]
    In_protein_A = [0 for i in range(nrows)]
    In_protein_B = [0 for i in range(nrows)]
    Decoy = ["FALSE" for i in range(nrows)]
    Modifications_A = ["" for i in range(nrows)]
    Modifications_B = ["" for i in range(nrows)]
    Confidence = ["High" for i in range(nrows)]

    # create annika dataframe
    annika_df = pd.DataFrame({"Checked": Checked,
                              "Crosslinker": Crosslinker,
                              "Crosslink Type": Crosslink_Type,
                              "# CSMs": CSMs,
                              "# Proteins": Proteins,
                              "Sequence A": Sequence_A,
                              "Accession A": Accession_A,
                              "Position A": Position_A,
                              "Sequence B": Sequence_B,
                              "Accession B": Accession_B,
                              "Position B": Position_B,
                              "Protein Descriptions A": Protein_Descriptions_A,
                              "Protein Descriptions B": Protein_Descriptions_B,
                              "Best CSM Score": Best_CSM_Score,
                              "In protein A": In_protein_A,
                              "In protein B": In_protein_B,
                              "Decoy": Decoy,
                              "Modifications A": Modifications_A,
                              "Modifications B": Modifications_B,
                              "Confidence": Confidence})

    return annika_df

if __name__ == "__main__":

    filename = sys.argv[1]
    df = create_annika_result(filename)
    df.to_excel(filename + "_annikaFormat.xlsx")
