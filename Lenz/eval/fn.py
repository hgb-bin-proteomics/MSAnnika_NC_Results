import pandas as pd

def get_proteins():
    return pd.read_excel("41467_2021_23666_MOESM10_ESM.xlsx", header = 1)

def get_annika_1():
    return pd.read_excel("Crosslinks_1%_nodecoy_1.xlsx")

def get_annika_2():
    return pd.read_excel("Crosslinks_1%_nodecoy_2.xlsx")

def get_xi():
    xi = pd.read_csv("xi2annika_crosslinks_validated.csv")
    print(f"Read xi results with {xi[xi['Class']=='Target'].shape[0]} target hits and {xi[xi['Class']=='Decoy'].shape[0]} decoy hits.")
    xi = xi[xi["Class"]=="Target"]
    return xi

def get_interactions(df):
    interactions = dict()
    for i, row in df.iterrows():
        a = str(row["Protein1_UniProt_Accession"]).strip()
        b = str(row["Protein2_UniProt_Accession"]).strip()
        if a not in interactions:
            interactions[a] = {b}
        else:
            interactions[a].add(b)
        if b not in interactions:
            interactions[b] = {a}
        else:
            interactions[b].add(a)
    return interactions

def is_fp(a, b, interactions):
    for p1 in a:
        for p2 in b:
            if p2 in interactions and p1 in interactions[p2]:
                return False
            if p1 in interactions and p2 in interactions[p1]:
                return False
    return True

def test_df(df, interactions, colname_a = "Accession A", colname_b = "Accession B"):
    interlinks = 0
    intralinks = 0
    fp = 0
    for i, row in df.iterrows():
        a = {x.strip() for x in str(row[colname_a]).split(";")}
        b = {x.strip() for x in str(row[colname_b]).split(";")}
        if len(a.intersection(b)) > 0:
            intralinks += 1
        else:
            if is_fp(a, b, interactions):
                fp += 1
            interlinks += 1

    print(f"Crosslinks: {df.shape[0]}")
    print(f"Intra + Inter: {intralinks + interlinks}")
    print(f"Intralinks: {intralinks}")
    print(f"Interlinks: {interlinks}")
    print(f"FDR: {fp / (df.shape[0] - fp)}")
    print(f"FDR (inter): {fp / (interlinks - fp)}")
