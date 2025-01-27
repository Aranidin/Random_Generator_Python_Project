import os
from Bio.Seq import Seq
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import matplotlib.pyplot as plt

# 1) Reading the DNA Sequence
def read_dna_from_text_file(file_path):

    print("Reading file: " + file_path)
    # Checks if the file exists at the specified path using the os.path.exists() function
    if not os.path.exists(file_path):
        raise FileNotFoundError("Reading file: " + file_path + "does not exist. Please check the file path.")

    # Opens the file at the specified file_path in read mode ('r')
    with open(file_path, 'r') as file:
        # Reads the file content, removes extra spaces, converts the sequence to uppercase
        sequence = file.read().strip().upper()

    valid_bases = {"A", "T", "C", "G"}
    # Converts the sequence into a set of unique characters, Checks if all characters in the sequence belong to the set of valid bases
    if not set(sequence).issubset(valid_bases):
        raise ValueError("Invalid DNA sequence: contains characters other than A, T, C, G.")

    return sequence


# 2) Calculating GC Content
def gc_content(sequence):

    gc_count = sequence.count('G') + sequence.count('C')
    total_bases = len(sequence)
    return round((gc_count / total_bases) * 100, 2)


# 3) Finding ORFs
# Locate all regions of the DNA that start with a start codon (ATG) and end with a stop codon (TAA, TAG, or TGA)
def find_orfs(sequence):

    start_codon = "ATG"
    stop_codons = {"TAA", "TAG", "TGA"}
    orfs = []

    for i in range(len(sequence)):
        if sequence[i:i + 3] == start_codon:
            for j in range(i + 3, len(sequence), 3):
                codon = sequence[j:j + 3]
                if codon in stop_codons:
                    orfs.append(sequence[i:j + 3])
                    break

    return orfs


# 4) Translating DNA to Protein
def translate_dna_to_protein(sequence):

    # Converts the input DNA sequence into a Biopython Seq object
    dna_seq = Seq(sequence)
    # Translates the Seq object into a protein sequence
    return str(dna_seq.translate(to_stop=True))


# 5) Protein Analysis
def analyze_protein(protein_sequence):

    # Creates a ProteinAnalysis object for the given protein sequence
    analysis = ProteinAnalysis(protein_sequence)
    return {
        "pI": round(analysis.isoelectric_point(), 2),
        "Molecular Weight": round(analysis.molecular_weight(), 2),
        "Amino Acid Percentages": analysis.get_amino_acids_percent()
    }


# 6) Visualizations

## GC Content (Bar chart)
def plot_gc_content(gc_percentage):
    """
    Visualize GC content as a bar chart.
    """
    plt.bar(["GC Content"], [gc_percentage])
    plt.title("GC Content of DNA Sequence")
    plt.ylabel("Percentage")
    plt.ylim(0, 100)
    plt.show()




## ORF Lengths (Histogram)
def plot_orf_lengths(orfs):

    lengths = [len(orf) for orf in orfs]
    if lengths:
        plt.hist(lengths, bins=range(0, max(lengths) + 30, 30), edgecolor="black", alpha=0.7)
        plt.title("Distribution of ORF Lengths")
        plt.xlabel("ORF Length (bases)")
        plt.ylabel("Frequency")
        plt.show()
    else:
        print("No ORFs to visualize.")


## Amino Acid Composition (Pie Chart)
def plot_amino_acid_composition(protein_sequence):

    # Counts the occurrences of each amino acid in the protein sequence using Counter
    from collections import Counter
    amino_acids = Counter(protein_sequence)
    # Separates amino acid labels and their counts for plotting
    labels, sizes = zip(*amino_acids.items())

    plt.pie(sizes, labels=labels, autopct="%1.1f%%", startangle=140)
    plt.title("Amino Acid Composition")
    plt.show()

## Processing of DNA file
def process_dna_file(file_path):
    dna_sequence = read_dna_from_text_file(file_path)
    gc_percentage = gc_content(dna_sequence)
    orfs = find_orfs(dna_sequence)
    protein_sequence = translate_dna_to_protein(dna_sequence)
    protein_properties = analyze_protein(protein_sequence)

    # Return all results as a dictionary
    return {
        "dna_sequence": dna_sequence,
        "gc_percentage": gc_percentage,
        "orfs": orfs,
        "protein_sequence": protein_sequence,
        "protein_analysis": protein_properties
    }


if __name__ == "__main__":
    # Step 1: Ask the user for the file path
    file_path = input("Enter the path to the DNA sequence file (e.g., dna_sequence.txt): ").strip()

    # Step 2: Read and validate the DNA sequence
    dna_sequence = read_dna_from_text_file(file_path)
    print("\nDNA Sequence Loaded Successfully!")

    # Step 3: Calculate the GC contents
    gc_percentage = gc_content(dna_sequence)
    print("\nGC Content: " + str(gc_percentage) + "%")
    plot_gc_content(gc_percentage)

    # Step 4: Find ORFs
    orfs = find_orfs(dna_sequence)
    print("\nIdentified ORFs:")
    if orfs:
        # Prints each ORF along with its index, visualizes the distribution of ORF lengths
        for index, orf in enumerate(orfs, 1):
            print("ORF " + str(index) + ": " + orf)
        plot_orf_lengths(orfs)
    else:
        print("No ORFs found in the sequence.")

    # Step 5: Translate DNA to protein
    protein_sequence = translate_dna_to_protein(dna_sequence)
    print("\nTranslated Protein Sequence:")
    print(protein_sequence)
    plot_amino_acid_composition(protein_sequence)

    # Step 6: Analyze the protein
    protein_properties = analyze_protein(protein_sequence)
    print("\nProtein Analysis:")
    print("  Isoelectric Point (pI): " + str(protein_properties['pI']))
    print("  Molecular Weight: " + str(protein_properties['Molecular Weight']) + " Da")
    print("  Amino Acid Composition (Percentages):")
    for aa, percentage in protein_properties["Amino Acid Percentages"].items():
        print("    " + aa + ": " + str(round(percentage * 100, 2)) + "%")
