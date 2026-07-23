from Bio.PDB import MMCIFParser
import urllib.request
import os
import warnings

# Suppress warnings for a clean terminal output
warnings.filterwarnings("ignore")

print("--- Starting 3D Protein Structure Analysis ---")

# 1. Setup directories
os.makedirs("data_raw/pdb", exist_ok=True)
os.makedirs("results", exist_ok=True)

# 2. Direct HTTPS download to bypass Biopython PDBList bugs
pdb_id = "1HHO"
cif_file = f"data_raw/pdb/{pdb_id}.cif"
url = f"https://files.rcsb.org/download/{pdb_id}.cif"

print(f"Downloading structure {pdb_id} directly from RCSB servers...")
try:
    urllib.request.urlretrieve(url, cif_file)
    print(f"Structure successfully downloaded to: {cif_file}")
except Exception as e:
    print(f"Network error occurred: {e}")
    exit()

# 3. Parse the structure
print("Parsing 3D coordinates and structure hierarchy...\n")
parser = MMCIFParser()
structure = parser.get_structure(pdb_id, cif_file)

# 4. Analyze the structure
model = list(structure.get_models())[0] # Select the first 3D model
chains = list(model.get_chains())
total_residues = 0
total_atoms = 0

print("--- Structure Summary ---")
print(f"Protein: {pdb_id}")
print(f"Number of Chains: {len(chains)}")

for chain in chains:
    residues = list(chain.get_residues())
    atoms = list(chain.get_atoms())
    total_residues += len(residues)
    total_atoms += len(atoms)
    print(f"  - Chain {chain.id}: {len(residues)} residues, {len(atoms)} atoms")

# 5. Save report
target_file = "results/task08_protein_structure_report.txt"
with open(target_file, "w") as f:
    f.write(f"3D PROTEIN STRUCTURE ANALYSIS: {pdb_id} (Human Hemoglobin)\n")
    f.write("="*55 + "\n")
    f.write(f"Total Chains: {len(chains)}\n")
    f.write(f"Total Residues (Amino Acids/Ligands/Water): {total_residues}\n")
    f.write(f"Total Atoms: {total_atoms}\n\n")
    f.write("Chain Breakdown:\n")
    for chain in chains:
        f.write(f"- Chain {chain.id}: {len(list(chain.get_residues()))} residues\n")

print(f"\nDone! Detailed structure report saved to: {target_file}")