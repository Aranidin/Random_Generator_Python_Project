from flask import Flask, render_template, request, flash
from dna_features_viewer import GraphicFeature, GraphicRecord
from bokeh.embed import components, file_html
from bokeh.resources import CDN
app = Flask(__name__)

#secret key
app.secret_key = "thinkofsomethingsecret"

@app.route('/', methods=['GET', 'POST'])
def home():
    record_script, record_div = None, None  # Initialize variables to prevent errors if no plot

    if request.method == 'POST':
        #see if NCBI form was submitted
        email = request.form.get("email")
        accession_number = request.form.get("acess_n")
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
                flash("Valid accession number", "validation")
        
        if dna_input:
            if not is_valid_dna_sequence(dna_input):
                flash("Invalid DNA sequence, only A, T, C, G, U, and N are allowed", "error")
            else: 
                flash("DNA sequence is valid!", "validation")

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
#Validation functions for NCBI input

def is_valid_email(email):
    if "@" in email and "." in email and len(email) > 2:
        return True
    return False
#something here doesn't work yet
def is_valid_accession_number(accession_number):
    if accession_number[3:].isdigit():
        return True
    return False
    
def is_valid_dna_sequence(dna_input):
    valid_bases = ["A", "T", "C", "G", "N", "U"]
    for i in dna_input:
        if i not in valid_bases:
            return False
    return True

if __name__ == '__main__':
    app.run(debug=True)
