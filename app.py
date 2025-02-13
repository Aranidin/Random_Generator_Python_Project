from flask import Flask, render_template, request, flash, jsonify, send_from_directory
from dna_features_viewer import GraphicFeature, CircularGraphicRecord 
import base64
from io import BytesIO
from Bio import SeqIO
from io import StringIO
import os
import requests
import matplotlib
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import json
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

app = Flask(__name__)
app.config["UPLOAD_FOLDER"] = "uploads" #folder for file uploads
app.secret_key = "thinkofsomethingsecret" #secret key
app.config["UPLOAD_FOLDER"] = "static/images"
os.makedirs(app.config["UPLOAD_FOLDER"], exist_ok=True)  # Ensure folder exists
# load annotations of DNA sequences
with open('snap_gene_features.txt') as json_file:
    snap_gene_feat = json.load(json_file)
#dna_input = None
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
            if not is_valid_dna_sequence(dna_input):
                flash("Invalid DNA sequence, only A, T, C, G, U, and N are allowed!", "error")
            elif len(dna_input) > 1000:
                flash("Sequence length exceeds 1000 bases, please reduce it!", "error")
            else:
                # Process the DNA sequence (whether from file or manual input)
                try:
                    # Perform sequence analysis (e.g., GC content, ORF detection, etc.)
                    results = processing_sequence(dna_input)
                    # add circular plot
                    image_path = generate_and_save_circular_plot(dna_input)
                    return render_template("results.html", **results, image_path=image_path)
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
            else:
                flash("Valid accession number!", "validation")
        
        if dna_input:
            if not is_valid_dna_sequence(dna_input):
                flash("Invalid DNA sequence, only A, T, C, G, U, and N are allowed!", "error")
            #max input of 1000 bases
            if len(dna_input) > 1000:
                flash("Only sequences with a max. length of 1000 allowed!")
            else:
                flash("Success!", "validation")
                try:
                    

                    results = processing_sequence(dna_input)
                    # Create DNA visualization
                    dna_fig = create_dna_figure(dna_input)
                    # Convert plot to JSON for embedding
                    plot_json = dna_fig.to_json()
                    results['dna_plot'] = plot_json
                    return render_template("results.html", **results)
                except Exception as e:
                    flash(f"Error processing DNA sequence: {e}", "error")

    
    return render_template('index.html', record_script=record_script, record_div=record_div)
def generate_and_save_circular_plot(dna_input):
    #  Generate a circular DNA plot and save it as a PNG file
    features=[]
    for name, seq in snap_gene_feat:
        start = dna_input.find(seq)
        if start != -1:
            end = start + len(seq)
            features.append(GraphicFeature(start=start, end=end, label=name, color="purple"))
    record = CircularGraphicRecord(sequence_length=len(dna_input), features=features)
    fig, _ = record.plot(figure_width=5)

    image_filename = "circular_plot.png"
    image_path = os.path.join(app.config["UPLOAD_FOLDER"], image_filename)
    fig.savefig(image_path, format="png")  # Save the figure
    return f"static/images/{image_filename}"  # Return relative path for HTML

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
    return "<h1>Functions</h1><p>Understand the different features and tools available.</p>"


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
@app.route('/load_features', methods=['GET'])
def load_features():
    try:
        with open('features.json') as json_file:
            features = json.load(json_file)
        return jsonify(features)
    except FileNotFoundError:
        return jsonify({})  # Return empty dict if file doesn't exist
    except Exception as e:
        return jsonify({"error": str(e)})
    
@app.route('/update_dna_view', methods=['POST'])
def update_dna_view():
    data = request.get_json()
    sequence = data['sequence']
    window_start = data['window_start']
    window_size = data['window_size']
    features = data.get('features', {})  # Get features if provided, empty dict if not
    
    fig = create_dna_figure(sequence, window_start, window_size, features)
    return jsonify(fig.to_dict())

if __name__ == '__main__':
    app.run(debug=True, threaded=False)
