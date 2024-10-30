#!/usr/bin/env python3

# Kojak with Percolator validated Result file to MS Annika Result file converter
# 2024 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

import re
import sys
import pandas as pd

__version = "1.0.5"
__date = "2024-10-30"

def get_peptides(row):
    peps = row["peptide"].split("--")
    pep1 = peps[0]
    pep2 = peps[1]
    pep1_pos = int(pep1.split("(")[1].split(")")[0])
    pep2_pos = int(pep2.split("(")[1].split(")")[0])
    pep1 = re.sub('[^a-zA-Z]+', '', pep1).upper()
    pep2 = re.sub('[^a-zA-Z]+', '', pep2).upper()
    return pep1, pep1_pos, pep2, pep2_pos

def create_annika_result(inter, intra, crosslinker = "DSS"):
    # load file
    df1 = pd.read_csv(inter, sep = "\t", names = ["PSMId", "score", "q-value", "posterior_error_prob", "peptide", "proteinId1", "proteinId2"], header = 0)
    df2 = pd.read_csv(intra, sep = "\t", names = ["PSMId", "score", "q-value", "posterior_error_prob", "peptide", "proteinIds"], header = 0)
    print("Read INTERLINKS:")
    print(df1)
    print(df1["peptide"])
    print("Read INTRALINKS:")
    print(df2)
    print(df2["peptide"])
    crosslinks = dict()
    for df in [df1, df2]:
        for i, row in df.iterrows():
            if row["q-value"] > 0.01:
                continue
            pep1, pep1_pos, pep2, pep2_pos = get_peptides(row)
            xl_key = "-".join(sorted([f"{pep1}{pep1_pos}", f"{pep2}{pep2_pos}"]))
            if xl_key not in crosslinks:
                crosslinks[xl_key] = {"pep1": pep1, "pep1_pos": pep1_pos, "pep2": pep2, "pep2_pos": pep2_pos, "score": row["score"], "csms": 1}
            else:
                if row["score"] > crosslinks[xl_key]["score"]:
                    crosslinks[xl_key]["score"] = row["score"]
                crosslinks[xl_key]["csms"] += 1
    print(f"Read {len(crosslinks)} crosslinks with q-value <= 0.01 from files {inter}, {intra}.")

    # columns
    Checked = ["FALSE" for crosslink in crosslinks]
    Crosslinker = [crosslinker for crosslink in crosslinks]
    Crosslink_Type = ["Intra" for crosslink in crosslinks]
    CSMs = [crosslinks[crosslink]["csms"] for crosslink in crosslinks]
    Proteins = [1 for crosslink in crosslinks]
    Sequence_A = [crosslinks[crosslink]["pep1"] for crosslink in crosslinks]
    Accession_A = ["Cas9" for crosslink in crosslinks]
    Position_A = [crosslinks[crosslink]["pep1_pos"] for crosslink in crosslinks]
    Sequence_B = [crosslinks[crosslink]["pep2"] for crosslink in crosslinks]
    Accession_B = ["Cas9" for crosslink in crosslinks]
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
