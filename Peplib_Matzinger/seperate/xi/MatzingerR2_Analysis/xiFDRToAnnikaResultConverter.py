#!/usr/bin/env python3

# xiFDR 2.2.1 Links Result file to MS Annika Result file converter
# 2024 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

import sys
import pandas as pd

__version = "1.0.1"
__date = "2024-06-28"

def get_crosslink_type(protein_idA, protein_idB):
    if protein_idA.upper() == protein_idB.upper():
        return "Intra"
    else:
        return "Inter"

def get_nr_proteins(xl_type):
    if xl_type == "Intra":
        return 1
    return 2

def get_peptide1(row):
    psmid = str(row["PSMIDs"]).split(";")[0]
    p1 = psmid.split("P1_")[1].split(" ")[0]
    if "." in p1:
        p1 = p1.split(".")[1]
    seq = ""
    for aa in p1:
        if aa.isupper():
            seq += aa
    if len(seq) == 0:
        raise RuntimeError("Couldn't parse peptide sequence!")
    return seq

def get_peptide2(row):
    psmid = str(row["PSMIDs"]).split(";")[0]
    p2 = psmid.split("P2_")[1].split(" ")[0]
    if "." in p2:
        p2 = p2.split(".")[1]
    seq = ""
    for aa in p2:
        if aa.isupper():
            seq += aa
    if len(seq) == 0:
        raise RuntimeError("Couldn't parse peptide sequence!")
    return seq

def get_positions(row):
    psmid = str(row["PSMIDs"]).split(";")[0]
    pos1 = psmid.split("P2_")[1].split(" ")[1]
    pos2 = psmid.split("P2_")[1].split(" ")[2]
    return [int(pos1), int(pos2)]

def create_annika_result(xifdr_filename, crosslinker = "DSS"):
    # load file
    df = pd.read_csv(xifdr_filename)

    #columns
    Crosslink_Type = []
    CSMs = []
    Sequence_A = []
    Accession_A = []
    Position_A = []
    Sequence_B = []
    Accession_B = []
    Position_B = []
    Best_CSM_Score = []

    nrows = 0
    for i, row in df.iterrows():
        if row["isDecoy"]:
            continue
        acc1 = str(row["Protein1"])
        acc2 = str(row["Protein2"])
        seq1 = get_peptide1(row)
        seq2 = get_peptide2(row)
        pos1, pos2 = get_positions(row)
        Crosslink_Type.append(get_crosslink_type(acc1, acc2))
        CSMs.append(row["count peptide pairs"])
        Sequence_A.append(seq1)
        Accession_A.append(acc1)
        Position_A.append(pos1)
        Sequence_B.append(seq2)
        Accession_B.append(acc2)
        Position_B.append(pos2)
        Best_CSM_Score.append(row["Score"])
        nrows += 1
    print(f"read {nrows} crosslinks from file.")

    # columns
    Checked = ["FALSE" for i in range(nrows)]
    Crosslinker = [crosslinker for i in range(nrows)]
    Proteins = [get_nr_proteins(xl_type) for xl_type in Crosslink_Type]
    Protein_Descriptions_A = ["" for i in range(nrows)]
    Protein_Descriptions_B = ["" for i in range(nrows)]
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
