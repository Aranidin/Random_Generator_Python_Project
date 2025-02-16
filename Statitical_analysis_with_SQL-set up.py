## Imports

import sqlite3
import os # Check if the file exists
from Bio.Seq import Seq
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import matplotlib.pyplot as plt

# Is used to count the frequency of amino acids in a protein sequence.
from collections import Counter


## Database Setup

# Connect to the SQLite database named 'Python_Project.db'
conn = sqlite3.connect('Python_Project.db')

# Create a cursor object to execute SQL commands
cursor = conn.cursor()

# Create 'dna_sequences' Table
# This table stores DNA sequences with their descriptions and timestamps
cursor.execute('''CREATE TABLE IF NOT EXISTS dna_sequences (
    id INTEGER PRIMARY KEY AUTOINCREMENT,  -- Unique ID for each DNA sequence
    sequence TEXT,  -- Stores the DNA sequence
    description TEXT,  -- Stores description of the sequence
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP  -- Timestamp of record creation
)''')

# Create 'gc_content' Table
# Stores GC content percentages, linked to the DNA sequences
cursor.execute('''CREATE TABLE IF NOT EXISTS gc_content (
    id INTEGER PRIMARY KEY AUTOINCREMENT,  -- Unique ID for each GC content record
    sequence_id INTEGER,  -- Links to DNA sequence via foreign key
    gc_percentage REAL,  -- Calculated GC percentage
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,  -- Timestamp of record creation
    FOREIGN KEY(sequence_id) REFERENCES dna_sequences(id) -- Links to 'dna_sequences' table
)''')

# Create 'orfs' Table
# Stores Open Reading Frames (ORFs) with their start and end positions
cursor.execute('''
    CREATE TABLE IF NOT EXISTS orfs (
        id INTEGER PRIMARY KEY AUTOINCREMENT, -- Unique ID for each ORF (auto-incrementing)
        sequence_id INTEGER,  -- Foreign key linking to a DNA sequence
        orf_sequence TEXT,  -- The nucleotide sequence representing the ORF
        start_position INTEGER,  -- Start position of the ORF in the sequence
        end_position INTEGER,  -- End position of the ORF in the sequence
        created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,  -- Timestamp of record creation
        FOREIGN KEY(sequence_id) REFERENCES dna_sequences(id)  -- Links to 'dna_sequences' table
)''')

# Create 'protein_sequences' Table
# Stores protein sequences derived from DNA sequences
cursor.execute('''
    CREATE TABLE IF NOT EXISTS protein_sequences (
        id INTEGER PRIMARY KEY AUTOINCREMENT,  -- Unique ID for each record
        sequence_id INTEGER,  -- Foreign key linking to DNA sequence
        protein_sequence TEXT,  -- Protein sequence from translation
        created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,  -- Timestamp of record creation
        FOREIGN KEY(sequence_id) REFERENCES dna_sequences(id)  -- Links to 'dna_sequences' table
)''')

# Create 'protein_analysis' Table
# Stores protein analysis data including pI and molecular weight
cursor.execute('''CREATE TABLE IF NOT EXISTS protein_analysis (
    id INTEGER PRIMARY KEY AUTOINCREMENT,  -- Unique ID for protein analysis
    sequence_id INTEGER,  -- Links to DNA sequence via foreign key
    pi REAL,  -- Isoelectric point
    molecular_weight REAL,  -- Molecular weight of the protein
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,  -- Timestamp of creation
    FOREIGN KEY(sequence_id) REFERENCES dna_sequences(id) -- Links to 'dna_sequences' table
)''')

# Create 'visualizations' Table
# Stores information about generated visualizations
cursor.execute('''CREATE TABLE IF NOT EXISTS visualizations (
    id INTEGER PRIMARY KEY AUTOINCREMENT,  -- Unique ID for each visualization
    sequence_id INTEGER,  -- Links to DNA sequence via foreign key
    chart_type TEXT,  -- Type of chart (e.g., Pie Chart)
    image_path TEXT,  -- Path to the saved image file
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,  -- Timestamp of record creation
    FOREIGN KEY(sequence_id) REFERENCES dna_sequences(id) -- Links to 'dna_sequences' table
)''')

# Save all table creation operations to the database
conn.commit()

## Functions for Saving Data

# Saves a DNA sequence and description to the database
def save_dna_sequence(sequence, description):
    # Executes SQL command to insert DNA sequence and description
    cursor.execute(
        "INSERT INTO dna_sequences (sequence, description) VALUES (?, ?)",
        (sequence, description)  # Passes values into placeholders
    )
    conn.commit()  # Commits transaction to save changes
    return cursor.lastrowid  # Returns the ID of the inserted sequence

# Saves GC content to the database
def save_gc_content(sequence_id, gc_percentage):
    # Inserts sequence ID and GC content into the 'gc_content' table
    cursor.execute(
        "INSERT INTO gc_content (sequence_id, gc_percentage) VALUES (?, ?)",
        (sequence_id, gc_percentage)
    )
    conn.commit()

# Saves ORFs with start and end positions to the database
def save_orfs(sequence_id, orfs):
    for orf_sequence, start, end in orfs:
        # Check if ORF already exists in the database
        cursor.execute('''
            SELECT id FROM orfs  -- Query to find existing ORF ID
            WHERE sequence_id = ? AND orf_sequence = ? AND start_position = ? AND end_position = ?
        ''', (sequence_id, orf_sequence, start, end)) # Use parameters to prevent SQL injection

        # Fetch the first matching record (if any)
        existing_orf = cursor.fetchone()

        # If no matching ORF found, insert a new one
        if not existing_orf:
            cursor.execute('''
                INSERT INTO orfs (sequence_id, orf_sequence, start_position, end_position)
                VALUES (?, ?, ?, ?)  -- Insert ORF details into the table
            ''', (sequence_id, orf_sequence, start, end))

    # Commit the changes to save them into the database
    conn.commit()

    # Print a message showing how many new ORFs were processed
    print(f"Saved {len(orfs)} new ORFs to the database.")

# Saves Protein Sequence to the Database
# Saves Protein Sequence to the Database
def save_protein_sequence(sequence_id, protein_sequence):
    cursor.execute(
        "INSERT INTO protein_sequences (sequence_id, protein_sequence) VALUES (?, ?)",
        (sequence_id, protein_sequence)
    )
    conn.commit()  # Saves the changes to the database
    print(f"Saved protein sequence for sequence ID {sequence_id} to the database.")


# Saves protein analysis results to the database
def save_protein_analysis(sequence_id, pi, molecular_weight):
    cursor.execute(
        "INSERT INTO protein_analysis (sequence_id, pi, molecular_weight) VALUES (?, ?, ?)",
        (sequence_id, pi, molecular_weight)
    )
    conn.commit()

# Saves chart details
def save_visualization(sequence_id, chart_type, image_path):
    cursor.execute(
        "INSERT INTO visualizations (sequence_id, chart_type, image_path) VALUES (?, ?, ?)",
        (sequence_id, chart_type, image_path)  # Adds chart metadata
    )
    conn.commit()

# Closes the database connection
def close_db():
    conn.close()  # Terminates the database connection
    print("Database connection closed.")



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

                    # Append the ORF (sequence and start-end positions) to the list
                    orfs.append((sequence[i:j + 3], i + 1, j + 3))
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

## Main script execution

# Step 1: Read DNA Sequence
file_path = input("Enter the path to the DNA sequence file (FASTA or TXT): ").strip()
dna_sequence = read_dna_sequence(file_path)
print("\nDNA Sequence Loaded Successfully!")

# Step 2: Save to Database & Show GC Content
sequence_id = save_dna_sequence(dna_sequence, "User-provided sequence")  # Saves sequence to the database
gc_percentage = gc_content(dna_sequence)  # Calculates GC content
print("\nGC Content:", gc_percentage, "%")  # Displays GC content
save_gc_content(sequence_id, gc_percentage)  # Stores GC content to the database

# Step 3: Save to Database & Find ORFs
orfs = find_orfs(dna_sequence)
if not orfs:
    print("No ORFs were found in the provided sequence.")
else:
    print("\nIdentified ORFs:")
    for index, (orf_sequence, start, end) in enumerate(orfs, 1):
        print(f"ORF {index}: Sequence: {orf_sequence}, Start: {start}, End: {end}")  # Displays each ORF
    save_orfs(sequence_id, orfs)  # Saves ORFs with start and end positions to the database
    print(f"Saved {len(orfs)} ORFs to the database with start and end positions.")

# Step 4: Translate to Protein & Store Protein Analysis
protein_sequence = translate_dna_to_protein(dna_sequence)
print("\nTranslated Protein Sequence:", protein_sequence)
protein_properties = analyze_protein(protein_sequence)
print("\nProtein Analysis:")
print("  Isoelectric Point (pI):", protein_properties['pI'])
print("  Molecular Weight:", protein_properties['Molecular Weight'], "Da")
save_protein_sequence(sequence_id, protein_sequence) # Save protein sequence to the database
save_protein_analysis(sequence_id, protein_properties['pI'],
                      protein_properties['Molecular Weight'])  # Stores protein analysis to the database

# Step 5: Plot & Save Amino Acid Composition
plot_amino_acid_composition(protein_sequence)
save_visualization(sequence_id, "Pie Chart",
                   "amino_acid_composition.png")  # Saves chart details to the database

# Step 6: Close Database
close_db()
