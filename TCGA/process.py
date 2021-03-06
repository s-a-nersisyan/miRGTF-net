import sys

import pandas as pd
import numpy as np


sample_sheet_file = sys.argv[1]
cart_path = sys.argv[2]

sample_sheet = pd.read_csv(sample_sheet_file, sep="\t")
# Select only primary tumors
sample_sheet = sample_sheet.loc[sample_sheet["Sample Type"] == "Primary Tumor"]
case_ids = np.unique(sample_sheet["Case ID"])

df_gene = pd.DataFrame()
df_miRNA = pd.DataFrame()

gencode = pd.read_csv("GENCODE_22.tsv", sep="\t", index_col=0)
miRBase = pd.read_csv("miRBase_22_1.tsv", sep="\t", index_col=1)

for case_id in case_ids:
    files = sample_sheet.loc[sample_sheet["Case ID"] == case_id]
    # Files should contain two rows: one for genes, second for miRNA
    if len(files) != 2:
        continue

    # Load paths for paired gene / miRNA expression files
    try:
        gene_row = files.loc[files["Data Type"] == "Gene Expression Quantification"]
        miRNA_row = files.loc[files["Data Type"] == "Isoform Expression Quantification"]
        gene_path = gene_row["File ID"].iloc[0] + "/" + gene_row["File Name"].iloc[0]
        miRNA_path = miRNA_row["File ID"].iloc[0] + "/" + miRNA_row["File Name"].iloc[0]
    except:
        continue

    # First, gene expression table
    df = pd.read_csv(cart_path + "/" + gene_path, sep="\t", header=None)
    df_gene["EnsemblID"] = df[0]
    df_gene[case_id] = df[1]

    # The, miRNA expression table
    df = pd.read_csv(cart_path + "/" + miRNA_path, sep="\t")
    df["MIMAT"] = [m.replace("mature,", "") if "MIMAT" in m else None for m in df["miRNA_region"]]
    # Sum up expression of all miRNA isoforms
    df = df[["MIMAT", "reads_per_million_miRNA_mapped"]].groupby("MIMAT").sum()

    df_miRNA["miRNA"] = miRBase["miRNA"]
    df_miRNA[case_id] = miRBase.join(df).fillna(0)["reads_per_million_miRNA_mapped"]

# Convert EnsemblIDs to Gene Symbols
df_gene = df_gene.set_index("EnsemblID").join(gencode).set_index("Gene Symbol")
# Set correct index
df_miRNA = df_miRNA.set_index("miRNA")

# Convert to TPM and apply log2-transformation
for df, out_name in zip([df_gene, df_miRNA], ["gene_expression.tsv", "miRNA_expression.tsv"]):
    # Filter by median >= 1 and remove duplicate genes
    index_name = df.index.name
    df = df.loc[df.median(axis=1).sort_values(ascending=False).index]
    df = df.reset_index().drop_duplicates(subset=index_name, keep='first').set_index(index_name)
    df = df.div(df.sum(axis=0), axis=1)*10**6  # TPM transformation
    df = np.log2(df + 1)
    df = df.loc[df.median(axis=1) >= 1]  # Filter low-expressed entries

    df.to_csv(out_name, sep="\t")
    print(df)
