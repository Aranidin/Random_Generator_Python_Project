# Checks if a file exists.
import os
from Bio.Seq import Seq
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import matplotlib.pyplot as plt

# Is used to count the frequency of amino acids in a protein sequence.
from collections import Counter


# 1) Reading the DNA Sequence

# Read a DNA sequence from a file (either FASTA or TXT format).
def read_dna_sequence(file_path):
    print("Reading file: " + file_path)

    # Check if the file exists. If not, raise a FileNotFoundError.
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"File {file_path} does not exist. Please check the file path.")

    # Open the file in read mode ('r') and read all lines.
    with open(file_path, "r") as file:
        lines = file.readlines()

        # Check if the file starts with '>', indicating it's in FASTA format.
        if lines[0].startswith(">"):
            print("Detected FASTA format. Ignoring header...")
            # Extract all lines that do not start with '>'.
            sequence_lines = [line.strip() for line in lines if not line.startswith(">")]
        else:
            print("Detected TXT format. Processing entire file...")

            # For TXT format, assume all lines contain the sequence.
            sequence_lines = [line.strip() for line in lines]

        # Join all lines into a single uppercase DNA sequence.
        sequence = "".join(sequence_lines).upper()

        # Validate the sequence by checking for invalid characters.
        valid_bases = {"A", "T", "C", "G"}
        invalid_chars = set(sequence) - valid_bases
        if invalid_chars:
            print(f"Warning: Invalid characters found in sequence: {invalid_chars}")

            # Remove any invalid characters from the sequence.
            sequence = ''.join(filter(lambda x: x in valid_bases, sequence))
            print("Cleaned sequence by removing invalid characters.")

        # Return the cleaned DNA sequence.
        return sequence


# 2) Calculating GC Content
def gc_content(sequence):
    gc_count = sequence.count('G') + sequence.count('C')
    return round((gc_count / len(sequence)) * 100, 2)


# 3) Finding ORFs

# Locate all regions of the DNA that start with a start codon (ATG) and end with a stop codon (TAA, TAG, or TGA).
def find_orfs(sequence):
    start_codon = "ATG"
    stop_codons = {"TAA", "TAG", "TGA"}
    orfs = []

    # Iterate through the sequence, checking every position for a start codon.
    for i in range(len(sequence)):
        if sequence[i:i + 3] == start_codon:
            # Once a start codon is found, look for the nearest stop codon in the same reading frame.
            for j in range(i + 3, len(sequence), 3):
                if sequence[j:j + 3] in stop_codons:
                    orfs.append(sequence[i:j + 3])
                    break
    if not orfs:
        print("No ORFs found in the provided DNA sequence.")

    return orfs


# 4) Translating DNA to Protein
def translate_dna_to_protein(sequence):
    # Trim the sequence to make length divisible by 3
    trimmed_sequence = sequence[:len(sequence) - (len(sequence) % 3)]

    # Convert the DNA sequence into a Biopython `Seq` object and translate it into a protein sequence.
    return str(Seq(trimmed_sequence).translate(to_stop=True))


# 5) Protein Analysis
def analyze_protein(protein_sequence):
    # Create a ProteinAnalysis object using the given protein sequence.
    analysis = ProteinAnalysis(protein_sequence)

    # Return a dictionary containing:
    # - The isoelectric point (rounded to two decimal places).
    # - The molecular weight (rounded to two decimal places).
    # - The amino acid composition as percentages.
    return {
        "pI": round(analysis.isoelectric_point(), 2),
        "Molecular Weight": round(analysis.molecular_weight(), 2),
        "Amino Acid Percentages": analysis.get_amino_acids_percent()
    }


# 6) Visualizations

## Amino Acid Composition (Pie Chart)
def plot_amino_acid_composition(protein_sequence):
    # Count the occurrences of each amino acid in the protein sequence.
    amino_acids = Counter(protein_sequence)

    sorted_amino_acids = sorted(amino_acids.items(), key=lambda x: x[1], reverse=True)
    labels, sizes = zip(*sorted_amino_acids)

    # Separate the amino acid labels and their counts for plotting.
    labels, sizes = zip(*sorted_amino_acids)

    # Create a pie chart showing the distribution of amino acids.
    full_names = {
        'A': 'Alanine', 'C': 'Cysteine', 'D': 'Aspartic Acid', 'E': 'Glutamic Acid', 'F': 'Phenylalanine',
        'G': 'Glycine', 'H': 'Histidine', 'I': 'Isoleucine', 'K': 'Lysine', 'L': 'Leucine',
        'M': 'Methionine', 'N': 'Asparagine', 'P': 'Proline', 'Q': 'Glutamine', 'R': 'Arginine',
        'S': 'Serine', 'T': 'Threonine', 'V': 'Valine', 'W': 'Tryptophan', 'Y': 'Tyrosine'
    }
    label_names = [f"{aa}: {full_names.get(aa, 'Unknown')} ({(count / sum(sizes)) * 100:.1f}%)" for aa, count in
                   sorted_amino_acids]
    plt.figure()
    plt.pie(sizes, labels=labels, autopct="%1.1f%%", startangle=140)
    plt.title("Amino Acid Composition")
    plt.legend(label_names, title="Amino Acids", loc="center left", bbox_to_anchor=(1, 0.5))
    plt.show()

# Main script execution.

if __name__ == "__main__":
    # Prompt the user for the path to the DNA file.
    file_path = input("Enter the path to the DNA sequence file (FASTA or TXT format): ").strip()

    # Step 1: Read the DNA sequence from the file.
    dna_sequence = read_dna_sequence(file_path)
    print("\nDNA Sequence Loaded Successfully!")

    # Step 2: Calculate and display GC content.
    gc_percentage = gc_content(dna_sequence)
    print("\nGC Content:", gc_percentage, "%")

    # Step 3: Identify and display ORFs in the DNA sequence.
    orfs = find_orfs(dna_sequence)
    if not orfs:
        print("No ORFs were found in the provided sequence.")
    else:
        print("\nIdentified ORFs:")
        for index, orf in enumerate(orfs, 1):
            print(f"ORF {index}: {orf}")

    # Step 4: Translate the DNA sequence into a protein sequence.
    protein_sequence = translate_dna_to_protein(dna_sequence)
    print("\nTranslated Protein Sequence:", protein_sequence)

    # Step 5: Visualize the amino acid composition of the protein.
    plot_amino_acid_composition(protein_sequence)

    # Step 6: Analyze and display protein properties.
    protein_properties = analyze_protein(protein_sequence)
    print("\nProtein Analysis:")
    print("  Isoelectric Point (pI):", protein_properties['pI'])
    print("  Molecular Weight:", protein_properties['Molecular Weight'], "Da")
    print("  Amino Acid Composition (Percentages):")
    for aa, percentage in protein_properties["Amino Acid Percentages"].items():
        print(f"    {aa}: {round(percentage * 100, 2)}%")