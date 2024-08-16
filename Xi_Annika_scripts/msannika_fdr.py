#!/usr/bin/env python3

# MS ANNIKA FDR VALIDATOR
# 2024 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

# version tracking
__version = "1.1.0"
__date = "2024-01-24"

# REQUIREMENTS
# pip install numpy
# pip install pandas
# pip install openpyxl

######################

"""
DESCRIPTION:
A script to group and validate results from MS Annika searches.
USAGE:
msannika_fdr.py f [f ...]
                  [-fdr FDR][--false_discovery_rate FDR]
                  [-h][--help]
                  [--version]
positional arguments:
  f                     MS Annika result files in Microsoft Excel format (.xlsx)
                        to process.
optional arguments:
  -fdr FDR, --false_discovery_rate FDR
                        False discovery rate to validate results for. Supports
                        both percentage input (e.g. 1) or fraction input (e.g.
                        0.01). By default not set and the input results will
                        just be grouped to crosslinks (if CSMs as input) or
                        nothing will be done (if crosslinks as input).
                        Default: None
  -h, --help            show this help message and exit
  --version             show program's version number and exit
"""

######################

import argparse
import numpy as np
import pandas as pd

from typing import List

class MSAnnika_CSM_Grouper:

    @staticmethod
    def get_crosslink_key(row: pd.Series) -> str:
        a = str(row["Sequence A"]) + str(row["Crosslinker Position A"])
        b = str(row["Sequence B"]) + str(row["Crosslinker Position B"])
        return "-".join(sorted([a, b]))

    @staticmethod
    def get_nr_proteins(row: pd.Series) -> int:
        proteins_str = str(row["Accession A"]).strip(";") + ";" + str(row["Accession B"]).strip(";")
        return len(proteins_str.split(";"))

    @staticmethod
    def get_best_csm_score(csms: List[pd.Series]) -> float:

        best_score = 0

        for csm in csms:
            if csm["Combined Score"] > best_score:
                best_score = csm["Combined Score"]

        return best_score

    @staticmethod
    def get_decoy_flag(row: pd.Series) -> bool:
        return True if "D" in row["Alpha T/D"] or "D" in row["Beta T/D"] else False

    @staticmethod
    def group(data: pd.DataFrame) -> pd.DataFrame:

        crosslinks = dict()

        for i, row in data.iterrows():
            crosslink = MSAnnika_CSM_Grouper.get_crosslink_key(row)
            if crosslink in crosslinks:
                crosslinks[crosslink].append(row)
            else:
                crosslinks[crosslink] = [row]

        rows = list()
        columns = ["CHECKED",
                   "Crosslinker",
                   "Crosslink Type",
                   "# CSMs",
                   "# Proteins",
                   "Sequence A",
                   "Accession A",
                   "Position A",
                   "Sequence B",
                   "Accession B",
                   "Position B",
                   "Protein Descriptions A",
                   "Protein Descriptions B",
                   "Best CSM Score",
                   "In protein A",
                   "In protein B",
                   "Decoy",
                   "Modifications A",
                   "Modifications B",
                   "Confidence"]

        for crosslink in crosslinks:
            row = pd.Series({"CHECKED": False,
                             "Crosslinker": crosslinks[crosslink][0]["Crosslinker"],
                             "Crosslink Type": crosslinks[crosslink][0]["Crosslink Type"],
                             "# CSMs": len(crosslinks[crosslink]),
                             "# Proteins": MSAnnika_CSM_Grouper.get_nr_proteins(crosslinks[crosslink][0]),
                             "Sequence A": crosslinks[crosslink][0]["Sequence A"],
                             "Accession A": crosslinks[crosslink][0]["Accession A"],
                             "Position A": crosslinks[crosslink][0]["Crosslinker Position A"],
                             "Sequence B": crosslinks[crosslink][0]["Sequence B"],
                             "Accession B": crosslinks[crosslink][0]["Accession B"],
                             "Position B": crosslinks[crosslink][0]["Crosslinker Position B"],
                             "Protein Descriptions A": crosslinks[crosslink][0]["Accession A"],
                             "Protein Descriptions B": crosslinks[crosslink][0]["Accession B"],
                             "Best CSM Score": MSAnnika_CSM_Grouper.get_best_csm_score(crosslinks[crosslink]),
                             "In protein A": crosslinks[crosslink][0]["A in protein"],
                             "In protein B": crosslinks[crosslink][0]["B in protein"],
                             "Decoy": MSAnnika_CSM_Grouper.get_decoy_flag(crosslinks[crosslink][0]),
                             "Modifications A": crosslinks[crosslink][0]["Modifications A"],
                             "Modifications B": crosslinks[crosslink][0]["Modifications B"],
                             "Confidence": "Unknown"})
            rows.append(row)

        return pd.concat(rows, ignore_index = True, axis = 1, names = columns).T

class MSAnnika_CSM_Validator:

    @staticmethod
    def get_class(row: pd.Series) -> str:
        return "Decoy" if "D" in row["Alpha T/D"] or "D" in row["Beta T/D"] else "Target"

    @staticmethod
    def get_cutoff(data: pd.DataFrame, fdr: float) -> float:

        data["Class"] = data.apply(lambda row: MSAnnika_CSM_Validator.get_class(row), axis = 1)
        data["Class_label"] = data.apply(lambda row: 0 if row["Class"] == "Target" else 1, axis = 1)
        labels = data["Class_label"].to_numpy()
        labels_sorted = labels[data["Combined Score"].to_numpy().argsort()]

        scores = sorted(data["Combined Score"].tolist())
        for i, score in enumerate(scores):
            if labels_sorted[i:].sum() / (labels_sorted[i:].shape[0] - labels_sorted[i:].sum()) < fdr:
                return score

        return scores[0]

    @staticmethod
    def validate(data: pd.DataFrame, fdr: float) -> pd.DataFrame:

        cutoff = MSAnnika_CSM_Validator.get_cutoff(data, fdr)
        df = data[data["Combined Score"] >= cutoff].copy()

        if "Confidence" not in df.columns:
            return df

        df["Confidence"] = df.apply(lambda row: "High", axis = 1)

        return df

class MSAnnika_Crosslink_Validator:

    @staticmethod
    def get_class(row: pd.Series) -> str:
        return "Decoy" if row["Decoy"] else "Target"

    @staticmethod
    def get_cutoff(data: pd.DataFrame, fdr: float) -> float:

        data["Class"] = data.apply(lambda row: MSAnnika_Crosslink_Validator.get_class(row), axis = 1)
        data["Class_label"] = data.apply(lambda row: 0 if row["Class"] == "Target" else 1, axis = 1)
        labels = data["Class_label"].to_numpy()
        labels_sorted = labels[data["Best CSM Score"].to_numpy().argsort()]

        scores = sorted(data["Best CSM Score"].tolist())
        for i, score in enumerate(scores):
            if labels_sorted[i:].sum() / (labels_sorted[i:].shape[0] - labels_sorted[i:].sum()) < fdr:
                return score

        return scores[0]

    @staticmethod
    def validate(data: pd.DataFrame, fdr: float) -> pd.DataFrame:

        cutoff = MSAnnika_Crosslink_Validator.get_cutoff(data, fdr)
        df = data[data["Best CSM Score"] >= cutoff].copy()
        df["Confidence"] = df.apply(lambda row: "High", axis = 1)

        return df

def main(argv = None) -> List[pd.DataFrame]:
    parser = argparse.ArgumentParser()
    parser.add_argument(metavar = "f",
                        dest = "files",
                        help = "Name/Path of the MS Annika result files to process.",
                        type = str,
                        nargs = "+")
    parser.add_argument("-fdr", "--false_discovery_rate",
                        dest = "fdr",
                        default = None,
                        help = "FDR for CSM/crosslink validation.",
                        type = float)
    parser.add_argument("--version",
                        action = "version",
                        version = __version)
    args = parser.parse_args(argv)

    if args.fdr is not None:
        print(f"Using {float(args.fdr) if float(args.fdr) < 1.0 else float(args.fdr) / 100.0} FDR.")

    result_list = list()

    for f, file in enumerate(args.files):
        df = pd.read_csv(file)

        if "Combined Score" in df.columns:
            crosslinks = MSAnnika_CSM_Grouper.group(df)
            crosslinks.to_csv(".csv".join(file.split(".csv")[:-1]) + "_crosslinks.csv", index = False)
            result_list.append(crosslinks)
            if args.fdr is not None:
                validated_csms = MSAnnika_CSM_Validator.validate(df, float(args.fdr) if float(args.fdr) < 1.0 else float(args.fdr) / 100.0)
                validated_csms.to_csv(".csv".join(file.split(".csv")[:-1]) + "_validated.csv", index = False)
                result_list.append(validated_csms)
                validated_crosslinks = MSAnnika_Crosslink_Validator.validate(crosslinks, float(args.fdr) if float(args.fdr) < 1.0 else float(args.fdr) / 100.0)
                validated_crosslinks.to_csv(".csv".join(file.split(".csv")[:-1]) + "_crosslinks_validated.csv", index = False)
                result_list.append(validated_crosslinks)
        else:
            if args.fdr is not None:
                validated_crosslinks = MSAnnika_Crosslink_Validator.validate(df, args.fdr if args.fdr < 1.0 else args.fdr / 100)
                validated_crosslinks.to_csv(".csv".join(file.split(".csv")[:-1]) + "_validated.csv", index = False)
                result_list.append(validated_crosslinks)
            else:
                print(f"Crosslink file without FDR given. Nothing to do. Skipping file {file}.")

        print(f"Processed {f + 1} files...")

    print("Done!")
    return result_list

if __name__ == "__main__":
    r = main()
