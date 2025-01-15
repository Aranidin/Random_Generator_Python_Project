from flask import Flask, render_template, request, flash
from dna_features_viewer import GraphicFeature, GraphicRecord
from bokeh.embed import components, file_html
from bokeh.resources import CDN
from Bio import SeqIO
import os
app = Flask(__name__)
app.config["UPLOAD_FOLDER"] = "uploads" #folder for file uploads
app.secret_key = "thinkofsomethingsecret" #secret key

@app.route('/', methods=['GET', 'POST'])
def home():
    record_script, record_div = None, None  # Initialize variables to prevent errors if no plot

    if request.method == 'POST':
        #NCBI form
        email = request.form.get("email")
        accession_number = request.form.get("access_n")
        fasta_file = request.files.get("fasta_file")

        #check for file upload
        if fasta_file:
            file_path = os.path.join(app.config["UPLOAD_FOLDER"], fasta_file.filename)
            fasta_file.save(file_path)

            #read sequence from fasta file
            dna_input = read_fasta(file_path)
            flash("Fasta file upload successful!", "validation")

        #manual DNA sequence input
        if not dna_input:
            dna_input = request.form.get('textarea')  # Get input from the text area
            print(dna_input)  # Print the input in the terminal for debugging

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

        return render_template('index.html')
        print(f"Email: {email}, Accession Number: {accession_number}, DNA Input: {dna_input}")

        # Example features for the DNA sequence visualization
        features = [
            GraphicFeature(start=0, end=20, strand=+1, color="#ffd700", label="Small feature"),
            GraphicFeature(start=20, end=500, strand=+1, color="#ffcccc", label="Gene 1"),
            GraphicFeature(start=400, end=700, strand=-1, color="#cffccc", label="Gene 2"),
            GraphicFeature(start=600, end=900, strand=+1, color="#ccccff", label="Gene 3")
        ]
        record = GraphicRecord(sequence_length=1000, features=features)
        record_p = record.plot_with_bokeh(figure_width=5)
        htm = file_html(record_p, CDN, "my plot")
        # Generate the script and div for embedding: Does not work!
        record_script, record_div = components(record_p)
    
    return render_template('index.html', record_script=record_script, record_div=record_div)

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
    if len(accession_number) >= 4 and accession_number[2:].isdigit():
        return True
    return False
#only allow the base Nucleotides, need to look at fasta files for upper, lower case  
def is_valid_dna_sequence(dna_input):
    valid_bases = ["A", "T", "C", "G", "N", "U"]
    for i in dna_input:
        if i not in valid_bases:
            return False
    return True

if __name__ == '__main__':
    app.run(debug=True)
