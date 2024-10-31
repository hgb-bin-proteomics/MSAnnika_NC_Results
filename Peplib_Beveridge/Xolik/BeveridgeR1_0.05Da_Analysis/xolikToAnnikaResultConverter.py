#!/usr/bin/env python3

# Xolik validated Result files to MS Annika Result file converter
# 2024 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

import re
import sys
import pandas as pd

__version = "1.0.0"
__date = "2024-10-30"

def get_peptides(row):
    pep1 = row["Peptide#1"].split(".")[1].strip()
    pep2 = row["Peptide#2"].split(".")[1].strip()
    pep1_pos = int(row["LinkSite#1"])
    pep2_pos = int(row["LinkSite#2"])
    pep1 = re.sub('[^a-zA-Z]+', '', pep1).upper()
    pep2 = re.sub('[^a-zA-Z]+', '', pep2).upper()
    return pep1, pep1_pos, pep2, pep2_pos

def create_annika_result(inter, intra, crosslinker = "DSS"):
    # load file
    df1 = pd.read_csv(inter)
    df2 = pd.read_csv(intra)
    print("Read INTERLINKS:")
    print(df1)
    print("Read INTRALINKS:")
    print(df2)
    crosslinks = dict()
    for df in [df1, df2]:
        for i, row in df.iterrows():
            if row["qValue"] > 0.01:
                continue
            pep1, pep1_pos, pep2, pep2_pos = get_peptides(row)
            xl_key = "-".join(sorted([f"{pep1}{pep1_pos}", f"{pep2}{pep2_pos}"]))
            if xl_key not in crosslinks:
                crosslinks[xl_key] = {"pep1": pep1, "pep1_pos": pep1_pos, "prot1": row["Protein#1"], "pep2": pep2, "pep2_pos": pep2_pos, "prot2": row["Protein#2"], "score": row["Score"], "csms": 1}
            else:
                if row["Score"] > crosslinks[xl_key]["score"]:
                    crosslinks[xl_key]["score"] = row["Score"]
                crosslinks[xl_key]["csms"] += 1
    print(f"Read {len(crosslinks)} crosslinks with q-value <= 0.01 from files {inter}, {intra}.")

    # columns
    Checked = ["FALSE" for crosslink in crosslinks]
    Crosslinker = [crosslinker for crosslink in crosslinks]
    Crosslink_Type = ["Intra" for crosslink in crosslinks]
    CSMs = [crosslinks[crosslink]["csms"] for crosslink in crosslinks]
    Proteins = [1 for crosslink in crosslinks]
    Sequence_A = [crosslinks[crosslink]["pep1"] for crosslink in crosslinks]
    Accession_A = [crosslinks[crosslink]["prot1"] for crosslink in crosslinks]
    Position_A = [crosslinks[crosslink]["pep1_pos"] for crosslink in crosslinks]
    Sequence_B = [crosslinks[crosslink]["pep2"] for crosslink in crosslinks]
    Accession_B = [crosslinks[crosslink]["prot2"] for crosslink in crosslinks]
    Position_B = [crosslinks[crosslink]["pep2_pos"] for crosslink in crosslinks]
    Protein_Descriptions_A = ["" for crosslink in crosslinks]
    Protein_Descriptions_B = ["" for crosslink in crosslinks]
    Best_CSM_Score = [crosslinks[crosslink]["score"] for crosslink in crosslinks]
    In_protein_A = [0 for crosslink in crosslinks]
    In_protein_B = [0 for crosslink in crosslinks]
    Decoy = ["FALSE" for crosslink in crosslinks]
    Modifications_A = ["" for crosslink in crosslinks]
    Modifications_B = ["" for crosslink in crosslinks]
    Confidence = ["High" for crosslink in crosslinks]

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

    filename1 = sys.argv[1]
    filename2 = sys.argv[2]
    df = create_annika_result(filename1, filename2)
    df.to_excel(filename1 + "_annikaFormat.xlsx")
