import os
from Bio.Seq import Seq
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import matplotlib.pyplot as plt

import plotly.graph_objects as go
from plotly.subplots import make_subplots
import json
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

    # Counts the occurrences of each amino acid in the protein sequence using Counter
    from collections import Counter
    amino_acids = Counter(protein_sequence)
    # Separates amino acid labels and their counts for plotting
    labels, sizes = zip(*amino_acids.items())

    plt.pie(sizes, labels=labels, autopct="%1.1f%%", startangle=140)
    plt.title("Amino Acid Composition")
    plt.show()

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

    # Return all results as a dictionary
    return {
        "dna_sequence": dna_sequence,
        "gc_percentage": gc_percentage,
        "orfs": orfs,
        "protein_sequence": protein_sequence,
        "protein_analysis": protein_properties
    }

##Processing sequence
def processing_sequence(dna_sequence, features_dict=None):

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
    # Return all results as a dictionary
    return {
        "dna_sequence": dna_sequence,
        "gc_percentage": gc_percentage,
        "orfs": orfs,
        "protein_sequence": protein_sequence,
        "protein_analysis": protein_properties,
        "dna_plot": dna_visualization # Added fo visualization
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
