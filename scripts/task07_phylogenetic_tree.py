from Bio import Entrez, SeqIO, Phylo
from Bio.Align import PairwiseAligner
from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor
import matplotlib.pyplot as plt
import os

print("--- Starting HBB Sequence Acquisition Pipeline ---")
# Set your academic email for Entrez access
Entrez.email = "dursunberkay.elci@agu.edu.tr"

# NCBI Accession numbers for HBB (Human, Chimp, Gorilla)
species_ids = {
    "Human": "NM_000518.5",
    "Chimp": "NM_001009044.1",
    "Gorilla": "NM_001009045.1"
}

records = []

# 1. Fetch sequences from NCBI
for species_name, accession_id in species_ids.items():
    print(f"Fetching {species_name} HBB sequence (ID: {accession_id})...")
    handle = Entrez.efetch(db="nucleotide", id=accession_id, rettype="fasta", retmode="text")
    seq_record = SeqIO.read(handle, "fasta")
    seq_record.id = species_name
    seq_record.description = f"HBB sequence for {species_name}"
    records.append(seq_record)
    handle.close()

print("\n--- Sequence Acquisition Summary ---")
for record in records:
    print(f"Species: {record.id:10} | Length: {len(record.seq):5} base pairs")

print("\n--- Calculating Evolutionary Distances ---")
# 2. Calculate Pairwise Distances
aligner = PairwiseAligner()
aligner.mode = 'global'

names = [rec.id for rec in records]
distances = []

for i in range(len(records)):
    row = []
    for j in range(i + 1):
        if i == j:
            row.append(0.0)
        else:
            # Hizalama skoru üzerinden uzaklık hesaplaması
            score = aligner.score(records[i].seq, records[j].seq)
            max_score = aligner.score(records[i].seq, records[i].seq)
            distance = 1 - (score / max_score)
            row.append(distance)
    distances.append(row)

# 3. Build the Phylogenetic Tree
print("Building the UPGMA Phylogenetic Tree...")
dist_matrix = DistanceMatrix(names, distances)
constructor = DistanceTreeConstructor()
tree = constructor.upgma(dist_matrix)

# 4. Print Tree to Console (ASCII)
print("\n--- Phylogenetic Tree (Terminal View) ---")
Phylo.draw_ascii(tree)

# 5. Save Tree as PNG
os.makedirs("results", exist_ok=True)
fig = plt.figure(figsize=(8, 5))
ax = fig.add_subplot(1, 1, 1)

# Ağacı çiz ve şekillendir
Phylo.draw(tree, axes=ax, do_show=False)
plt.title("Phylogenetic Tree of HBB Gene (Human, Chimp, Gorilla)")

target_file = "results/task07_phylogenetic_tree.png"
plt.savefig(target_file, dpi=300, bbox_inches='tight')
print(f"\nDone! High-resolution tree saved to: {target_file}")