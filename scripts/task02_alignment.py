import os
import subprocess
from Bio import Entrez, SeqIO, AlignIO, Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
import matplotlib.pyplot as plt

def fetch_sequences(accessions, email, output_file):
    """Fetches protein sequences from NCBI and saves them to a FASTA file."""
    Entrez.email = email
    records = []
    
    for acc in accessions:
        print(f"Fetching {acc}...")
        handle = Entrez.efetch(db="protein", id=acc, rettype="fasta", retmode="text")
        record = SeqIO.read(handle, "fasta")
        # Simplify the species name for the tree plot
        record.id = record.description.split("[")[1].replace("]", "") if "[" in record.description else acc
        records.append(record)
        handle.close()
        
    SeqIO.write(records, output_file, "fasta")
    print(f"Saved unaligned sequences to {output_file}")

def main():
    # 1. Configuration
    email = "your.email@student.agu.edu.tr" 
    
    # NCBI Protein Accession Codes for Hemoglobin Alpha Subunit
    accessions = [
        "NP_000508.1",      # Homo sapiens (Human)
        "NP_001009232.1",   # Pan troglodytes (Chimp)
        "NP_001188325.1",   # Macaca mulatta (Macaque)
        "NP_032246.2",      # Mus musculus (Mouse)
        "NP_036706.2",      # Rattus norvegicus (Rat)
        "NP_001003310.1"    # Canis lupus familiaris (Dog)
    ]
    
    unaligned_file = "results/hemoglobin_unaligned.fasta"
    aligned_file = "results/hemoglobin_aligned.fasta"
    tree_image = "results/phylogenetic_tree.png"
    
    os.makedirs("results", exist_ok=True)

    # 2. Download Sequences
    print("\n--- STEP 1: Downloading Data ---")
    fetch_sequences(accessions, email, unaligned_file)

    # 3. Align Sequences with MUSCLE using Python's subprocess
    print("\n--- STEP 2: Sequence Alignment ---")
    muscle_exe = os.path.abspath("scripts/muscle.exe") if os.name == "nt" else "muscle"
    
    if not os.path.exists(muscle_exe):
        print(f"ERROR: '{muscle_exe}' not found! Please check your scripts folder.")
        return

    print("Running MUSCLE alignment engine...")
    try:
        # Try modern MUSCLE v5 syntax first (-align and -output)
        subprocess.run([muscle_exe, "-align", unaligned_file, "-output", aligned_file], check=True, capture_output=True)
    except subprocess.CalledProcessError:
        try:
            # Fallback to older MUSCLE v3 syntax (-in and -out) if v5 fails
            subprocess.run([muscle_exe, "-in", unaligned_file, "-out", aligned_file], check=True, capture_output=True)
        except Exception as e:
            print(f"Alignment failed. Error details: {e}")
            return
            
    print(f"Alignment complete. Saved to {aligned_file}")

    # 4. Construct Phylogenetic Tree
    print("\n--- STEP 3: Building Phylogenetic Tree ---")
    alignment = AlignIO.read(aligned_file, "fasta")
    
    calculator = DistanceCalculator('identity')
    distance_matrix = calculator.get_distance(alignment)
    
    constructor = DistanceTreeConstructor()
    tree = constructor.nj(distance_matrix) # Neighbor-Joining Method
    
    # 5. Visualize and Save
    fig = plt.figure(figsize=(10, 6))
    axes = fig.add_subplot(1, 1, 1)
    
    # Styling the tree
    Phylo.draw(tree, axes=axes, do_show=False, branch_labels=None)
    plt.title("Evolutionary Tree of Hemoglobin Alpha", fontsize=16)
    
    plt.savefig(tree_image, dpi=300, bbox_inches="tight")
    print(f"\nSUCCESS! Tree visualization saved as '{tree_image}'.")

if __name__ == "__main__":
    main()