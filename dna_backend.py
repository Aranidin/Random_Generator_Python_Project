import os
from Bio.Seq import Seq
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import matplotlib.pyplot as plt
from dna_features_viewer import GraphicFeature, CircularGraphicRecord, GraphicRecord
import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import json
import sqlite3
from collections import Counter
# Is used to count the frequency of amino acids in a protein sequence.
from collections import Counter
## Database Setup

# Connect to the SQLite database named 'Python_Project.db'
conn = sqlite3.connect('Python_Project.db', check_same_thread=False)

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
def read_dna_from_text_file(file_path):
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
        print("Extracted Sequence:", sequence)  # Debugging

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
# Locate all regions of the DNA that start with a start codon (ATG) and end with a stop codon (TAA, TAG, or TGA)
def find_orfs(sequence):

    start_codon = "ATG"
    stop_codons = {"TAA", "TAG", "TGA"}
    orfs = []

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
    """# Converts the input DNA sequence into a Biopython Seq object
    dna_seq = Seq(sequence)
    # Translates the Seq object into a protein sequence
    return str(dna_seq.translate(to_stop=True))"""


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
# Linear DNA with annotations (Dinara)
# Function for interactive linear DNA visualization
def create_dna_figure(sequence, window_start=0, window_size=100, features_dict=None):
    # Calculate window end 
    window_end = min(window_start + window_size, len(sequence))
    visible_sequence = sequence[window_start:window_end]
    
    # Calculate positions for sequence display
    x_positions = list(range(window_start + 1, window_end + 1))
    y_positions = [0] * len(visible_sequence)
    
    fig = make_subplots(rows=1, cols=1, shared_xaxes=True)
    
    # DNA sequence bases
    for i, base in enumerate(visible_sequence):
        color = {
            'A': '#2CD05E',  # Green
            'T': '#CC2C2C',  # Red
            'G': '#FFAA00',  # Orange
            'C': '#0066FF'   # Blue
        }.get(base.upper(), '#808080')
        
        fig.add_trace(
            go.Scatter(
                x=[x_positions[i]],
                y=[y_positions[i]],
                text=[base],
                mode="text",
                textfont=dict(
                    size=18,
                    family="Courier New",
                    color=color
                ),
                hoverinfo='text',
                hovertext=f'Position: {window_start + i + 1}<br>Base: {base}',
                showlegend=False
            )
        )
    
    # Add feature annotations
    if features_dict and isinstance(features_dict, dict):
        colors = ['#FF9999', '#99FF99', '#9999FF', '#FFFF99', '#FF99FF', '#99FFFF']
        
        for idx, (name, feat_seq) in enumerate(features_dict.items()):
            if not isinstance(feat_seq, str):
                continue
                
            # Find all occurrences of the feature sequence
            current_pos = 0
            while True:
                # Find the next occurrence of the feature sequence
                start = sequence.find(feat_seq, current_pos)
                if start == -1:  # No more occurrences found
                    break
                    
                end = start + len(feat_seq)
                current_pos = start + 1  # Move to next position for next search
                
                # Only show features that are visible in current window
                if not (end < window_start or start > window_end):
                    color = colors[idx % len(colors)]
                    
                    # Adjust feature boundaries to window
                    visible_start = max(start, window_start)
                    visible_end = min(end, window_end)
                    
                    # Add feature box
                    fig.add_trace(
                        go.Scatter(
                            x=[visible_start + 1, visible_start + 1, visible_end + 1, visible_end + 1, visible_start + 1],
                            y=[0.3, 0.7, 0.7, 0.3, 0.3],
                            fill="toself",
                            fillcolor=color,
                            line=dict(color=color),
                            name=f"{name} ({start + 1}-{end})",
                            hoverinfo='text',
                            hovertext=f'Feature: {name}<br>Start: {start + 1}<br>End: {end}<br>Length: {len(feat_seq)} bp<br>Sequence: {feat_seq}'
                        )
                    )
                    
                    # Add feature label if there's enough space
                    if visible_end - visible_start > len(name):
                        fig.add_trace(
                            go.Scatter(
                                x=[(visible_start + visible_end)/2 + 1],
                                y=[0.5],
                                mode='text',
                                text=[name],
                                textposition='middle center',
                                textfont=dict(
                                    size=10,
                                    color='black'
                                ),
                                showlegend=False,
                                hoverinfo='skip'
                            )
                        )
    
    # Update layout
    fig.update_layout(
        plot_bgcolor='white',
        paper_bgcolor='white',
        height=200,
        margin=dict(l=50, r=50, t=50, b=50),
        xaxis=dict(
            showgrid=False,
            zeroline=False,
            title='Base Position',
            range=[window_start + 1, window_end + 1]
        ),
        yaxis=dict(
            showgrid=False,
            zeroline=False,
            showticklabels=False,
            range=[-0.2, 1.0]
        ),
        showlegend=True,
        legend=dict(
            orientation="h",
            yanchor="bottom",
            y=1.02,
            xanchor="right",
            x=1
        ),
        title=dict(
            text=f'DNA Sequence Viewer (Position {window_start + 1} to {window_end})',
            x=0.5,
            xanchor='center'
        )
    )
    
    return fig
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
    #plt.show()
    aa_plot_path = os.path.join("static/images", "amino_acid_composition.png")
    plt.savefig(aa_plot_path, bbox_inches='tight', dpi=100, format="png")
    plt.close()
    return "amino_acid_composition.png"
##Processing dna file
def process_dna_file(file_path):
    print("Processing DNA file:", file_path)  # Debugging

    dna_sequence = read_dna_from_text_file(file_path)
    print("Extracted DNA Sequence:", dna_sequence)  # Debugging

    try:
        gc_percentage = gc_content(dna_sequence)
        print("GC Content:", gc_percentage)
    except Exception as e:
        print("Error in gc_content:", e)
        return {"error": f"gc_content error: {e}"}

    try:
        orfs = find_orfs(dna_sequence)
        print("ORFs Found:", orfs)
    except Exception as e:
        print("Error in find_orfs:", e)
        return {"error": f"find_orfs error: {e}"}

    try:
        protein_sequence = translate_dna_to_protein(dna_sequence)
        print("Protein Sequence:", protein_sequence)
    except Exception as e:
        print("Error in translate_dna_to_protein:", e)
        return {"error": f"translate_dna_to_protein error: {e}"}

    try:
        protein_properties = analyze_protein(protein_sequence)
        print("Protein Properties:", protein_properties)
    except Exception as e:
        print("Error in analyze_protein:", e)
        return {"error": f"analyze_protein error: {e}"}
    ### ADDITION: Save Results to Database ###
    try:
        # Save DNA Sequence
        sequence_id = save_dna_sequence(dna_sequence, "User_Provided Sequence")

        # Save GC Content
        save_gc_content(sequence_id, gc_percentage)

        # Save ORFs
        if orfs:
            save_orfs(sequence_id, orfs)

        # Save Protein Sequence
        save_protein_sequence(sequence_id, protein_sequence)

        # Save Protein Analysis
        save_protein_analysis(sequence_id,
                              protein_properties['pI'],
                              protein_properties['Molecular Weight'])

        # Save Amino Acid Composition Chart Path
        save_visualization(sequence_id, "Amino Acid Composition", aa_plot_path)

        print(f"✅ Results for File Sequence {sequence_id} saved to database.")

    except Exception as e:
        print(f"❌ Database Save Error: {e}")

    ###

    # Return all results as a dictionary
    return {
        "dna_sequence": dna_sequence,
        "gc_percentage": gc_percentage,
        "orfs": orfs,
        "protein_sequence": protein_sequence,
        "protein_analysis": protein_properties
    }

##Processing sequence
def processing_sequence(dna_sequence, features_dict=None, upload_path_lin = None, upload_path_circ= None):

    try:
        gc_percentage = gc_content(dna_sequence)
        print("GC Content:", gc_percentage)
    except Exception as e:
        print("Error in gc_content:", e)
        return {"error": f"gc_content error: {e}"}

    try:
        orfs = find_orfs(dna_sequence)
        print("ORFs Found:", orfs)
    except Exception as e:
        print("Error in find_orfs:", e)
        return {"error": f"find_orfs error: {e}"}

    try:
        protein_sequence = translate_dna_to_protein(dna_sequence)
        print("Protein Sequence:", protein_sequence)
    except Exception as e:
        print("Error in translate_dna_to_protein:", e)
        return {"error": f"translate_dna_to_protein error: {e}"}
    
    try:
        protein_properties = analyze_protein(protein_sequence)
        aa_plot_path = plot_amino_acid_composition(protein_sequence)
        print("Protein Properties:", protein_properties)
    except Exception as e:
        print("Error in analyze_protein:", e)
        return {"error": f"analyze_protein error: {e}"}
    # Preprocessing for visualization and annotation (Dinara)
    try:
        # Create DNA visualization
        #dna_fig = create_dna_figure(dna_sequence)
        dna_fig = create_dna_figure(dna_sequence, features_dict=features_dict)
        # Convert plot to JSON for embedding
        dna_visualization = dna_fig.to_json()
    except Exception as e:
        print("Error in dna_visualization:", e)
        return {"error": f"dna_visualization error: {e}"}

    ### STORE RESULTS INTO DATABASE (ADDITIONALLY) ###
    try:
        # Create a unique entry for the DNA sequence
        sequence_id = save_dna_sequence(dna_sequence, "User-Provided Sequence")

        # Save GC Content
        save_gc_content(sequence_id, gc_percentage)

        # Save ORFs with start & stop positions
        if orfs:
            save_orfs(sequence_id, orfs)

        # Save Protein Sequence
        save_protein_sequence(sequence_id, protein_sequence)

        # Save Protein Analysis
        save_protein_analysis(sequence_id,
                              protein_properties['pI'],
                              protein_properties['Molecular Weight'])

        # Save Amino Acid Composition Chart Path
        save_visualization(sequence_id, "Amino Acid Composition", aa_plot_path)

        # Save DNA Visualization Path
        save_visualization(sequence_id, "DNA Plot", "dna_visualization.json")

        print(f"✅ Results for Sequence {sequence_id} saved to database.")

    except Exception as e:
        print("Error saving to database:", e)

    # Return all results as a dictionary

    return {
        "dna_sequence": dna_sequence,
        "gc_percentage": gc_percentage,
        "orfs": orfs,
        "protein_sequence": protein_sequence,
        "protein_analysis": protein_properties,
        "dna_plot": dna_visualization, # Added for interactive visualization,
        "aa_plot": aa_plot_path

    }


if __name__ == "__main__":
    file_path = input("Enter the path to the DNA sequence file (FASTA or TXT): ").strip()
    dna_sequence = read_dna_from_text_file(file_path)
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

    #   Step 6: Close Database
    close_db()

