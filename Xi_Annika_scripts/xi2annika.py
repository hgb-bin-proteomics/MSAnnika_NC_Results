import pandas as pd

def get_unique(df: pd.DataFrame) -> pd.DataFrame:
    unique_scans = dict()
    for i, row in df.iterrows():
        scan = row["ScanTitle"]
        score = row["match score"]
        if scan in unique_scans:
            if unique_scans[scan]["match score"] < score:
                unique_scans[scan] = row
        else:
            unique_scans[scan] = row
    if len(unique_scans.keys()) != len(df["ScanTitle"].unique()):
        raise RuntimeError("Couldn't determine unique scans!")
    else:
        print(f"Found {len(unique_scans)} unique scans.")
    return pd.concat(list(unique_scans.values()), ignore_index = True, axis = 1, names = df.columns.tolist()).T

def xiSearch2Annika(xi: pd.DataFrame, unique_scans = True) -> pd.DataFrame:
    if unique_scans:
        df = get_unique(xi)
    else:
        df = xi
    df.rename(columns = {"Run": "Spectrum File",
                         "Scan": "First Scan",
                         "PrecoursorCharge": "Charge",
                         "BasePeptide1": "Sequence A",
                         "Protein1": "Accession A",
                         "Link1": "Crosslinker Position A",
                         "BasePeptide2": "Sequence B",
                         "Protein2": "Accession B",
                         "Link2": "Crosslinker Position B",
                         "Fasta1": "Protein Descriptions A",
                         "Fasta2": "Protein Descriptions B",
                         "Start1": "A in protein",
                         "Start2": "B in protein",
                         "match score": "Combined Score",
                         "Pep1Score": "Score Alpha",
                         "Pep2Score": "Score Beta"},
                  inplace = True,
                  errors = "raise")
    df.dropna(subset = ["Peptide1", "Peptide2"], inplace = True)
    df["Accession A"] = df["Accession A"].apply(lambda x: str(x).replace("REV_", ""))
    df["Accession B"] = df["Accession B"].apply(lambda x: str(x).replace("REV_", ""))
    df["Crosslink Type"] = df.apply(lambda x: "Intra" if len(set(x["Accession A"].split(";")).intersection(set(x["Accession B"].split(";")))) > 0 else "Inter", axis = 1)
    df["Alpha T/D"] = df["Protein1decoy"].apply(lambda x: "T" if int(x) == 0 else "D")
    df["Beta T/D"] = df["Protein2decoy"].apply(lambda x: "T" if int(x) == 0 else "D")
    df["Modifications A"] = df.apply(lambda x: "Not Parsed", axis = 1)
    df["Modifications B"] = df.apply(lambda x: "Not Parsed", axis = 1)
    columns_to_keep = ["Spectrum File", "First Scan", "Sequence A", "Accession A",
                       "Crosslinker Position A", "Sequence B", "Accession B",
                       "Crosslinker Position B", "Protein Descriptions A",
                       "Protein Descriptions B", "A in protein", "B in protein",
                       "Combined Score", "Score Alpha", "Score Beta", "Crosslinker",
                       "Crosslink Type", "Alpha T/D", "Beta T/D", "Modifications A",
                       "Modifications B", "Charge"]
    df.drop(df.columns.difference(columns_to_keep), axis = 1, inplace = True)
    return df
