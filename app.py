from flask import Flask, render_template, request
from dna_features_viewer import GraphicFeature, GraphicRecord
from bokeh.embed import components, file_html
from bokeh.resources import CDN
app = Flask(__name__)

@app.route('/', methods=['GET', 'POST'])
def home():
    record_script, record_div = None, None  # Initialize variables to prevent errors if no plot

    if request.method == 'POST':
        dna_input = request.form.get('textarea')  # Get input from the text area
        print(dna_input)  # Print the input in the terminal for debugging

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

if __name__ == '__main__':
    app.run(debug=True)
