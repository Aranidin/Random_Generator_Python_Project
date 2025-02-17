from flask import Flask, render_template, request, flash, jsonify, send_from_directory, redirect
from dna_features_viewer import GraphicFeature, CircularGraphicRecord, GraphicRecord
from Bio import SeqIO
from Bio.Seq import Seq
from io import StringIO
import os
import requests
import matplotlib
import matplotlib.pyplot as plt
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import json
import numpy as np
from Bio.Restriction import *
import primer_design_module as pr_d

from dna_backend import (
    read_dna_from_text_file,
    gc_content,
    find_orfs,
    translate_dna_to_protein,
    analyze_protein,
    plot_gc_content,
    plot_orf_lengths,
    plot_amino_acid_composition,
    process_dna_file,
    processing_sequence,
    create_dna_figure,
)
import main_primer
from Bio.SeqUtils import MeltingTemp as mt


matplotlib.use('Agg') # for import of matplotlib non-interactive plots

app = Flask(__name__)
app.config["UPLOAD_FOLDER"] = "uploads" #folder for file uploads
app.secret_key = "thinkofsomethingsecret" #secret key
app.config["UPLOAD_FOLDER"] = "static/images"
os.makedirs(app.config["UPLOAD_FOLDER"], exist_ok=True)  # Ensure folder exists

# load annotations of DNA sequences
with open('snap_gene_features.txt') as json_file:
    snap_gene_feat = json.load(json_file)

# folders for static linear and circular plots
linear_path = os.path.join(app.config["UPLOAD_FOLDER"], "dna_linear_record.png")
linear_path_circ = os.path.join(app.config["UPLOAD_FOLDER"], "dna_circular_record.png")

@app.route('/', methods=['GET', 'POST'])
def home():
    record_script, record_div = None, None  #  prevent errors if no plot
    dna_input = None
    results = {}  # Initialize results dictionary to store processing outputs
    
    if request.method == 'POST':
        print("Post request received") #Debugging

        #get form data 
        email = request.form.get("email")
        accession_number = request.form.get("access_n")
        fasta_file = request.files.get("fasta_file")
        print(f"Uploaded file: {fasta_file}")  # Debugging

        if email and accession_number:
            try:
                #Fetch DNA sequence from NCBI
                dna_input = fetch_sequence_from_ncbi(accession_number, email)
                print(f"Fetched sequence: {dna_input[:50]}...") #Displays part of the sequence

                #Validate DNA sequnce
                if not is_valid_dna_sequence(dna_input):
                    flash("Invalid DNA sequence, only A, T, C, G, U, and N are allowed!", "error")
                elif len(dna_input) > 1000:
                    flash("Sequence length exceeds 1000 bases, please reduce it!", "error")
                else:
                    #perform sequence analysis
                    results = processing_sequence(dna_input)
                    # Annotate DNA Features
                    features = []
                    for name, seq in snap_gene_feat.items():
                        start = dna_input.find(seq)
                        if start != -1:
                            end = start + len(seq)
                            features.append(GraphicFeature(
                            start=start,
                            end=end,
                            strand=+1,
                            color="#ffd700",
                            label=name
                            ))
                    # Make linear static DNA Plot with annotations 
                    # (Adapted from DnaFeaturesViewer README)
                    record = GraphicRecord(sequence_length=len(dna_input), features=features)
                    fig, (ax, ax2) = plt.subplots(2,1,figsize=(10, 3), sharex=True, gridspec_kw={"height_ratios": [4, 1]})
                    record.plot(ax=ax)

                    # Visualize local GC content  (Adapted from DnaFeaturesViewer README)
                    gc_local = lambda s: 100.0 * len([c for c in s if c in "GC"]) / 50
                    xx = np.arange(len(dna_input) - 50)
                    yy = [gc_local(dna_input[x : x + 50]) for x in xx]
                    ax2.fill_between(xx + 25, yy, alpha=0.3)
                    ax2.set_ylim(bottom=0)
                    ax2.set_ylabel("GC(%)")
                    linear_path = os.path.join(app.config["UPLOAD_FOLDER"], "dna_linear_record.png")
                    plt.savefig(linear_path, dpi=300, bbox_inches="tight", format="png")
                    plt.close(fig)

                    # Make circular DNA Plot with annotations
                    # (Adapted from DnaFeaturesViewer README)
                    record_circ = CircularGraphicRecord(sequence_length=len(dna_input), features=features)
                    fig_circ, ax_circ= plt.subplots(figsize=(10,2))
                    record_circ.plot(ax=ax_circ)
                    linear_path_circ = os.path.join(app.config["UPLOAD_FOLDER"], "dna_circular_record.png")
                    plt.savefig(linear_path_circ, dpi=200, bbox_inches="tight", format="png")
                    plt.close(fig_circ)

                    
                    results['linear_plot'] = 'dna_linear_record.png'
                    results['dna_sequence'] = dna_input
                    results['features'] = features
                    results['circular_plot'] = 'dna_circular_record.png'
                    return render_template("results.html", **results)
            except Exception as e:
                flash(f"Error fetching NCBI sequence: {e}", "error")                      

        #check for file upload 
        if fasta_file and fasta_file.filename !="":
            #Debugging : Log file type and size
            fasta_file.stream.seek(0) # Ensures pointer is at the beginning
            file_size = len(fasta_file.read())
            print(f"File type: {fasta_file.content_type}, File size: {file_size}") #Debigging: File type and length

            if file_size == 0:
                print("Uploaded File is empty!")
                flash("Uploaded file is empty!", "error")
            else:
               fasta_file.stream.seek(0) #reset position before saving
               file_path = os.path.join(app.config["UPLOAD_FOLDER"], fasta_file.filename)
               fasta_file.save(file_path)
               print(f"File received: {fasta_file.filename}")  # Debugging log to check if file is received
               
               #validate file content (fallback)
               try:  
                #Using utility function from dna_backend
                flash("FASTA file upload and parsing successful!", "validation") 
                dna_input = read_fasta(file_path)
                results = process_dna_file(file_path)
                features = []
                for name, seq in snap_gene_feat.items():
                    start = dna_input.find(seq)
                    if start != -1:
                        end = start + len(seq)
                        features.append(GraphicFeature(
                        start=start,
                        end=end,
                        strand=+1,
                        color="#ffd700",
                        label=name
                        ))
                # Make linear static DNA Plot with annotations 
                # (Adapted from DnaFeaturesViewer README)
                record = GraphicRecord(sequence_length=len(dna_input), features=features)
                fig, (ax, ax2) = plt.subplots(2,1,figsize=(10, 3), sharex=True, gridspec_kw={"height_ratios": [4, 1]})
                record.plot(ax=ax)

                # Visualize local GC content  (Adapted from DnaFeaturesViewer README)
                gc_local = lambda s: 100.0 * len([c for c in s if c in "GC"]) / 50
                xx = np.arange(len(dna_input) - 50)
                yy = [gc_local(dna_input[x : x + 50]) for x in xx]
                ax2.fill_between(xx + 25, yy, alpha=0.3)
                ax2.set_ylim(bottom=0)
                ax2.set_ylabel("GC(%)")
                linear_path = os.path.join(app.config["UPLOAD_FOLDER"], "dna_linear_record.png")
                plt.savefig(linear_path, dpi=300, bbox_inches="tight", format="png")
                plt.close(fig)

                # Make circular DNA Plot with annotations
                # (Adapted from DnaFeaturesViewer README)
                record_circ = CircularGraphicRecord(sequence_length=len(dna_input), features=features)
                fig_circ, ax_circ= plt.subplots(figsize=(10,3))
                record_circ.plot(ax=ax_circ)
                linear_path_circ = os.path.join(app.config["UPLOAD_FOLDER"], "dna_circular_record.png")
                plt.savefig(linear_path_circ, dpi=300, bbox_inches="tight", format="png")
                plt.close(fig_circ)

                    
                results['linear_plot'] = 'dna_linear_record.png'
                results['dna_sequence'] = dna_input
                results['features'] = features
                results['circular_plot'] = 'dna_circular_record.png'
                return(render_template("results.html", **results))
               except Exception as e:
                   flash(f"Error reading FASTA file: {e}", "error")
                   #flash("Invalid FASTA file content. Please check your file and try again!", "error")


        #manual DNA sequence input
        if not dna_input:
            dna_input = request.form.get('textarea')  # Get input from the text area
            print(dna_input)  # Print the input in the terminal for debugging
        # Validate the DNA sequence input
        if dna_input:
            features = []
            for name, seq in snap_gene_feat.items():
                start = dna_input.find(seq)
                if start != -1:
                    end = start + len(seq)
                    features.append(GraphicFeature(
                    start=start,
                    end=end,
                    strand=+1,  # +1 for forward strand
                    color="#ffd700",
                    label=name
                    ))
            if not is_valid_dna_sequence(dna_input):
                print("Invalid DNA sequence, only A, T, C, G, U, and N are allowed!", "error")
            elif len(dna_input) > 1000:
                flash("Sequence length exceeds 1000 bases, please reduce it!", "error")
            else:
                # Process the DNA sequence (whether from file or manual input)
                try:
                    # Perform sequence analysis (e.g., GC content, ORF detection, etc.)
                    results = processing_sequence(dna_input)
                    # Annotate DNA Features
                    features = []
                    for name, seq in snap_gene_feat.items():
                        start = dna_input.find(seq)
                        if start != -1:
                            end = start + len(seq)
                            features.append(GraphicFeature(
                            start=start,
                            end=end,
                            strand=+1,
                            color="#ffd700",
                            label=name
                            ))
                    # Make linear static DNA Plot with annotations 
                    # (Adapted from DnaFeaturesViewer README)
                    record = GraphicRecord(sequence_length=len(dna_input), features=features)
                    fig, (ax, ax2) = plt.subplots(2,1,figsize=(10, 3), sharex=True, gridspec_kw={"height_ratios": [4, 1]})
                    record.plot(ax=ax)

                    # Visualize local GC content  (Adapted from DnaFeaturesViewer README)
                    gc_local = lambda s: 100.0 * len([c for c in s if c in "GC"]) / 50
                    xx = np.arange(len(dna_input) - 50)
                    yy = [gc_local(dna_input[x : x + 50]) for x in xx]
                    ax2.fill_between(xx + 25, yy, alpha=0.3)
                    ax2.set_ylim(bottom=0)
                    ax2.set_ylabel("GC(%)")
                    linear_path = os.path.join(app.config["UPLOAD_FOLDER"], "dna_linear_record.png")
                    plt.savefig(linear_path, dpi=300, bbox_inches="tight", format="png")
                    plt.close(fig)

                    # Make circular DNA Plot with annotations
                    # (Adapted from DnaFeaturesViewer README)
                    record_circ = CircularGraphicRecord(sequence_length=len(dna_input), features=features)
                    fig_circ, ax_circ= plt.subplots(figsize=(10,3))
                    record_circ.plot(ax=ax_circ)
                    linear_path_circ = os.path.join(app.config["UPLOAD_FOLDER"], "dna_circular_record.png")
                    plt.savefig(linear_path_circ, dpi=300, bbox_inches="tight", format="png")
                    plt.close(fig_circ)

                    
                    results['linear_plot'] = 'dna_linear_record.png'
                    results['dna_sequence'] = dna_input
                    results['features'] = features
                    results['circular_plot'] = 'dna_circular_record.png'
                    return render_template("results.html", **results)
                except Exception as e:
                    flash(f"Error processing DNA sequence: {e}", "error")
        # Validate inputs
        if email:
            if not is_valid_email(email):
                flash("Invalid email address, try again!", "error")
            else:
                flash("Valid email address!", "validation")

        if accession_number:
            if not is_valid_accession_number(accession_number):
                flash("Invalid accession number, try again!", "error")
            # else:
            #     flash("Valid accession number!", "validation")
        
        if dna_input:
            if not is_valid_dna_sequence(dna_input):
                flash("Invalid DNA sequence, only A, T, C, G, U, and N are allowed!", "error")
            #max input of 1000 bases
            if len(dna_input) > 1000:
                flash("Only sequences with a max. length of 1000 allowed!")
            else:
                #flash("Success!", "validation")
                try:
                    #results = processing_sequence(dna_input)
                    results = processing_sequence(dna_input)

                    # Create DNA visualization
                    dna_fig = create_dna_figure(dna_input, features_dict=snap_gene_feat)
                    # Convert plot to JSON for embedding
                    plot_json = dna_fig.to_json()
                    results['dna_plot'] = plot_json
                    results = processing_sequence(dna_input)
                    # Annotate DNA Features
                    features = []
                    for name, seq in snap_gene_feat.items():
                        start = dna_input.find(seq)
                        if start != -1:
                            end = start + len(seq)
                            features.append(GraphicFeature(
                            start=start,
                            end=end,
                            strand=+1,
                            color="#ffd700",
                            label=name
                            ))
                    # Make linear static DNA Plot with annotations 
                    # (Adapted from DnaFeaturesViewer README)
                    record = GraphicRecord(sequence_length=len(dna_input), features=features)
                    fig, (ax, ax2) = plt.subplots(2,1,figsize=(10, 3), sharex=True, gridspec_kw={"height_ratios": [4, 1]})
                    record.plot(ax=ax)

                    # Visualize local GC content  (Adapted from DnaFeaturesViewer README)
                    gc_local = lambda s: 100.0 * len([c for c in s if c in "GC"]) / 50
                    xx = np.arange(len(dna_input) - 50)
                    yy = [gc_local(dna_input[x : x + 50]) for x in xx]
                    ax2.fill_between(xx + 25, yy, alpha=0.3)
                    ax2.set_ylim(bottom=0)
                    ax2.set_ylabel("GC(%)")
                    linear_path = os.path.join(app.config["UPLOAD_FOLDER"], "dna_linear_record.png")
                    plt.savefig(linear_path, dpi=300, bbox_inches="tight", format="png")
                    plt.close(fig)

                    # Make circular DNA Plot with annotations
                    # (Adapted from DnaFeaturesViewer README)
                    record_circ = CircularGraphicRecord(sequence_length=len(dna_input), features=features)
                    fig_circ, ax_circ= plt.subplots(figsize=(10,3))
                    record_circ.plot(ax=ax_circ)
                    linear_path_circ = os.path.join(app.config["UPLOAD_FOLDER"], "dna_circular_record.png")
                    plt.savefig(linear_path_circ, dpi=300, bbox_inches="tight", format="png")
                    plt.close(fig_circ)

                    
                    results['linear_plot'] = 'dna_linear_record.png'
                    results['dna_sequence'] = dna_input
                    results['features'] = features
                    results['circular_plot'] = 'dna_circular_record.png'
                    return render_template("results.html", **results, features = snap_gene_feat)
                except Exception as e:
                    flash(f"Error processing DNA sequence: {e}", "error")

    
    return render_template('index.html', record_script=record_script, record_div=record_div)

@app.route('/circ_plot', methods=['POST'])
def generate_circ_plot():
    data = request.get_json()
    dna_input = data.get('sequence', '')
    
    if not dna_input:
        return jsonify({"error": "No DNA sequence provided"})
    
    try:
        # Generate features for the plot
        features = []
        for name, seq in snap_gene_feat.items():
            start = dna_input.find(seq)
            if start != -1:
                end = start + len(seq)
                features.append(GraphicFeature(
                    start=start,
                    end=end,
                    strand=+1,
                    color="#ffd700",
                    label=name
                ))
        
        # Create the circular plot
        record = CircularGraphicRecord(sequence_length=len(dna_input), features=features)
        fig, ax = plt.subplots(figsize=(8, 8))  # Made square and larger for better circular visualization
        ax.set_title("Circular DNA Map")
        record.plot(ax=ax)
        
        # Save to the same directory as linear plot
        image_path = os.path.join(app.config["UPLOAD_FOLDER"], "dna_circular_record.png")
        plt.savefig(image_path, dpi=300, bbox_inches="tight", format="png")
        plt.close(fig)
        
        return jsonify({"success": True})
    except Exception as e:
        return jsonify({"error": str(e)})
@app.route('/lin_plot', methods=['POST'])
def generate_lin_plot():
    data = request.get_json()
    dna_input = data.get('sequence', '')
    
    if not dna_input:
        return jsonify({"error": "No DNA sequence provided"})
    
    try:
        # Generate features for the plot
        features = []
        for name, seq in snap_gene_feat.items():
            start = dna_input.find(seq)
            if start != -1:
                end = start + len(seq)
                features.append(GraphicFeature(
                    start=start,
                    end=end,
                    strand=+1,  # +1 for forward strand
                    color="#ffd700",
                    label=name
                ))
        
        # Create the linear plot
        record = GraphicRecord(sequence_length=len(dna_input), features=features)
        fig, ax = plt.subplots(figsize=(5, 5))
        record.plot(ax=ax)
        image_path = os.path.join(app.config["UPLOAD_FOLDER"], "dna_linear_record.png")
        plt.savefig(image_path, dpi=300, bbox_inches="tight", format="png")
        plt.close(fig)
        return jsonify({"image_path": image_path})
    except Exception as e:
        return jsonify({"error": str(e)})

@app.route('/static/images/<filename>')
def get_image(filename):
    return send_from_directory(app.config["UPLOAD_FOLDER"], filename)
@app.route('/help')
def help_page():
    return render_template("help.html")
@app.route('/help/biology-basics')
def biology_basics():
    return render_template("biology_basics.html")

@app.route('/help/example')
def help_examples():
    return render_template("help_example.html")

@app.route('/help/functions')
def help_functions():
    return render_template("functions.html")


@app.route('/save_features', methods=['POST'])
def save_features():
    features = request.json
    try:
        with open('features.json', 'w') as f:
            json.dump(features, f)
        return jsonify({"status": "success"})
    except Exception as e:
        return jsonify({"status": "error", "message": str(e)})

#Function to read a fasta file and parse DNA sequence
def read_fasta(file_path):
    record = SeqIO.read(file_path, "fasta")
    return str(record.seq)

#Validation functions for NCBI input
def is_valid_email(email):
    if "@" in email and "." in email and len(email) > 2:
        return True
    return False
#acc number must be min 4 long and contain digits from the third position
def is_valid_accession_number(accession_number):
    if len(accession_number) >= 4:
        return True
    return False
#only allow the base Nucleotides, need to look at fasta files for upper, lower case  
def is_valid_dna_sequence(dna_input):
    valid_bases = ["A", "T", "C", "G", "N", "U", "a", "t", "c", "g", "n", "u"]
    dna_input = dna_input.replace("\n", "").replace("\r", "").replace(" ", "")  # Remove unwanted characters
    for i in dna_input:
        if i not in valid_bases:
            return False
    return True

def fetch_sequence_from_ncbi(accession_number, email):
    #NCBI Entrez E-utilities URL
    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"

    params = {
        "db": "nucleotide",  #search in NCBI nucleotide database
        "term": accession_number,  # Accessionnumber to search for
        "retmode" : "xml",  #XML format
        "email" : email  #email for Entrez usage tracking from NCBI
    }
    response = requests.get(url, params=params)
    if response.status_code == 200:
        #Parse XML response to get the ID
        from xml.etree import ElementTree as ET
        tree = ET.ElementTree(ET.fromstring(response.text))
        root = tree.getroot()
        id_list = root.find(".//IdList")

        if id_list is not None:
            #get the first ID from the list
            sequence_id = id_list.find("Id").text
            #get sequence using ID
            fetch_url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
            fetch_params = {
                "db": "nucleotide",
                "id": sequence_id,
                "rettype": "fasta",
                "retmode": "text",
                "email": email
            }
            fetch_response = requests.get(fetch_url, params = fetch_params)
            if fetch_response.status_code == 200:
                record = SeqIO.read(StringIO(fetch_response.text), "fasta")
                return str(record.seq) #return sequence as string
            else:
                raise Exception("Failed to detch sequence data!")
        else:
            raise Exception("Accession number not found!")
    else:
        raise Exception("Failed to search for accession number!")


#to be able to go back and forth in the visualized sequence
# + visualize annotated features

@app.route('/update_dna_view', methods=['POST'])
def update_dna_view():
    try:
        data = request.get_json()
        sequence = data.get('sequence', '')
        window_start = data.get('window_start', 0)
        window_size = data.get('window_size', 150)
        
        if not sequence:
            return jsonify({"error": "No sequence provided"})
        
        # Create updated visualization
        fig = create_dna_figure(
            sequence=sequence,
            window_start=window_start,
            window_size=window_size,
            features_dict=snap_gene_feat
        )
        
        return jsonify(fig.to_dict())
    except Exception as e:
        return jsonify({"error": str(e)})

@app.route('/load_features', methods=['GET'])
def load_features():
    try:
        with open('snap_gene_features.txt') as json_file:
            features = json.load(json_file)
        return jsonify(features)
    except FileNotFoundError:
        return jsonify({})
    except Exception as e:
        return jsonify({"error": str(e)})


@app.route('/primer_design', methods=['POST'])
def primer_design():
    try:
        dna_sequence = request.form['dna_sequence']
        forward_primer = request.form['forward_primer']
        reverse_primer = request.form['reverse_primer']

        # Validate primers
        if not all(base.upper() in 'ATCG' for base in forward_primer):
            flash('Forward primer contains invalid bases. Use only A, T, C, G.', 'error')
            return redirect(request.referrer)

        if not all(base.upper() in 'ATCG' for base in reverse_primer):
            flash('Reverse primer contains invalid bases. Use only A, T, C, G.', 'error')
            return redirect(request.referrer)

        # Check if primers exist in sequence
        if forward_primer.upper() not in dna_sequence.upper():
            flash('Forward primer not found in the DNA sequence.', 'error')
            return redirect(request.referrer)

        if reverse_primer.upper() not in dna_sequence.upper():
            flash('Reverse primer not found in the DNA sequence.', 'error')
            return redirect(request.referrer)

        # Create primer design object
        primer_design = pr_d.primer_design(Seq(dna_sequence), forward_primer, reverse_primer)
        
        # Get primer analysis results
        analysis_results = primer_design.primer_analysis()
        amplicon = primer_design.amplicon
        
        # Perform restriction analysis
        lysis = pr_d.lysis_analysis(amplicon)

        # Prepare analysis results
        primer_results = {
            'forward_primer': {
                'sequence': forward_primer,
                'length': len(forward_primer),
                'gc_content': gc_content(forward_primer),
                'tm': round(mt.Tm_NN(forward_primer), 1)
            },
            'reverse_primer': {
                'sequence': reverse_primer,
                'length': len(reverse_primer),
                'gc_content': gc_content(reverse_primer),
                'tm': round(mt.Tm_NN(reverse_primer), 1)
            },
            'amplicon': str(amplicon),
            'amplicon_length': len(amplicon),
            'restriction_analysis': {
                'single_cut': lysis.single_cut_enzymes,
                'multi_cut': lysis.multi_cut_enzymes
            }
        }

        return render_template('primer_results.html', 
                             results=primer_results, 
                             dna_sequence=dna_sequence)

    except Exception as e:
        flash(f'Error in primer analysis: {str(e)}', 'error')
        return redirect(request.referrer)

### Database ###
import sqlite3
import csv
import zipfile
from flask import send_file, flash, redirect

@app.route('/download_results', methods=['GET'])
def download_results():
    try:
        # Ensure database connection
        conn = sqlite3.connect('Python_Project.db')
        cursor = conn.cursor()

        # Query all data from the database
        cursor.execute("""
            SELECT d.id, d.sequence, d.description, d.created_at, 
                   g.gc_percentage, 
                   o.orf_sequence, o.start_position, o.end_position, 
                   p.protein_sequence, 
                   pa.pi, pa.molecular_weight
            FROM dna_sequences d
            LEFT JOIN gc_content g ON d.id = g.sequence_id
            LEFT JOIN orfs o ON d.id = o.sequence_id
            LEFT JOIN protein_sequences p ON d.id = p.sequence_id
            LEFT JOIN protein_analysis pa ON d.id = pa.sequence_id
        """)
        records = cursor.fetchall()

        # Create CSV file
        csv_path = "static/downloads/dna_analysis_results.csv"
        with open(csv_path, mode='w', newline='') as csv_file:
            writer = csv.writer(csv_file)
            writer.writerow(["Sequence ID", "DNA Sequence", "Description", "Created At",
                             "GC Content (%)", "ORF Sequence", "ORF Start", "ORF End",
                             "Protein Sequence", "pI", "Molecular Weight"])
            for row in records:
                writer.writerow(row)

        # Create TXT file
        txt_path = "static/downloads/dna_analysis_results.txt"
        with open(txt_path, 'w') as txt_file:
            txt_file.write("DNA Analysis Results\n")
            txt_file.write("="*50 + "\n\n")
            for record in records:
                txt_file.write(f"Sequence ID: {record[0]}\n")
                txt_file.write(f"DNA Sequence: {record[1]}\n")
                txt_file.write(f"Description: {record[2]}\n")
                txt_file.write(f"Created At: {record[3]}\n")
                txt_file.write(f"GC Content (%): {record[4]}\n")
                txt_file.write(f"ORF Sequence: {record[5]} (Start: {record[6]}, End: {record[7]})\n")
                txt_file.write(f"Protein Sequence: {record[8]}\n")
                txt_file.write(f"Isoelectric Point (pI): {record[9]}\n")
                txt_file.write(f"Molecular Weight: {record[10]} Da\n")
                txt_file.write("\n" + "-"*50 + "\n\n")

        # Zip the results
        zip_path = "static/downloads/dna_analysis_results.zip"
        with zipfile.ZipFile(zip_path, 'w') as zipf:
            zipf.write(csv_path, arcname="dna_analysis_results.csv")
            zipf.write(txt_path, arcname="dna_analysis_results.txt")

        return send_file(zip_path, as_attachment=True, download_name="dna_analysis_results.zip")

    except Exception as e:
        flash(f"Error generating download: {e}", "error")
        return redirect('/')

    finally:
        cursor.close()
        conn.close()
###

if __name__ == '__main__':
    app.run(debug=True, threaded=False)
