#!/usr/bin/env python3

# MaxLynx (MaxQuant 2.6.2.0) crosslinkMsms Result file to MS Annika Result file converter
# 2024 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

import sys
import pandas as pd

__version = "1.0.0"
__date = "2024-06-25"

def get_crosslink_type(protein_idA, protein_idB):
    if protein_idA.upper() == protein_idB.upper():
        return "Intra"
    else:
        return "Inter"

def get_proteins(csms):
    return [str(csms[0]["Proteins1"]).strip(), str(csms[0]["Proteins2"]).strip()]

def get_nr_proteins(xl_type):
    if xl_type == "Intra":
        return 1
    return 2

def get_score(csms):
    max_score = 0.0
    for csm in csms:
        score = float(csm["Score"])
        if score > max_score:
            max_score = score
    return max_score

def create_annika_result(maxlynx_filename, crosslinker = "DSS"):
    # load file
    df = pd.read_csv(maxlynx_filename, sep = "\t")
    crosslinks = dict()
    for i, row in df.iterrows():
        if not pd.isna(row["Decoy"]):
            continue
        if str(row["Crosslink product type"]).strip() not in ["Intra-protein link", "Inter-protein link"]:
            continue
        seq1 = str(row["Sequence1"]).strip()
        seq2 = str(row["Sequence2"]).strip()
        pos1 = str(int(row["Peptide index of Crosslink 1"]))
        pos2 = str(int(row["Peptide index of Crosslink 2"]))
        crosslink = "-".join(sorted([seq1 + pos1, seq2 + pos2]))
        if crosslink in crosslinks:
            crosslinks[crosslink].append(row)
        else:
            crosslinks[crosslink] = [row]
    print(f"read {len(crosslinks)} crosslinks from file.")
    nrows = len(crosslinks)

    # columns
    Checked = ["FALSE" for i in range(nrows)]
    Crosslinker = [crosslinker for i in range(nrows)]
    Crosslink_Type = [get_crosslink_type(get_proteins(crosslinks[crosslink])[0], get_proteins(crosslinks[crosslink])[1]) for crosslink in crosslinks]
    CSMs = [len(crosslinks[crosslink]) for crosslink in crosslinks]
    Proteins = [get_nr_proteins(xl_type) for xl_type in Crosslink_Type]
    Sequence_A = [str(crosslinks[crosslink][0]["Sequence1"]).strip() for crosslink in crosslinks]
    Accession_A = [get_proteins(crosslinks[crosslink])[0] for crosslink in crosslinks]
    Position_A = [int(crosslinks[crosslink][0]["Peptide index of Crosslink 1"]) for crosslink in crosslinks]
    Sequence_B = [str(crosslinks[crosslink][0]["Sequence2"]).strip() for crosslink in crosslinks]
    Accession_B = [get_proteins(crosslinks[crosslink])[1] for crosslink in crosslinks]
    Position_B = [int(crosslinks[crosslink][0]["Peptide index of Crosslink 2"]) for crosslink in crosslinks]
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
