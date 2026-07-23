import GEOparse
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os

print("Bioinformatics Analysis Starting...")

os.makedirs("data_raw", exist_ok=True)
os.makedirs("results", exist_ok=True)

# 1. Download Microarray compatible data from NCBI GEO
print("Downloading GSE25724 dataset from NCBI...")
gse = GEOparse.get_GEO(geo="GSE25724", destdir="./data_raw")

# 2. Convert gene expression values into a Pandas DataFrame
print("Creating DataFrame...")
df = gse.pivot_samples('VALUE')

# 3. Find the top 10 Differentially Expressed (most variable) genes
print("Calculating top 10 most variable genes...")
top_10_genes = df.var(axis=1).sort_values(ascending=False).head(10).index
df_top10 = df.loc[top_10_genes]

# 4. Heatmap plotting with Seaborn
print("Drawing heatmap...")
plt.figure(figsize=(10, 8))
sns.heatmap(df_top10, cmap="viridis", annot=False)

plt.title("GSE25724 - Top 10 Most Variable Genes")
plt.xlabel("Samples")
plt.ylabel("Gene IDs")
plt.tight_layout()

# 5. Save the result visually to the results folder
target_file = "./results/task01_heatmap.png"
plt.savefig(target_file)
print(f"Done! Heatmap successfully saved: {target_file}")