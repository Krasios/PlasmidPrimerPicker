from flask import Flask,request,render_template,flash
from werkzeug.utils import secure_filename
from parser import parse
import tempfile
import os

#from suggester import runSuggester


UPLOAD_DIR = tempfile.gettempdir()
ALLOWED_FILETYPES = set(['gb'])

app = Flask(__name__)
app.secret_key = os.urandom(24)

def valid_file(file):
    return '.' in file and file.split('.')[-1].lower() in ALLOWED_FILETYPES
 
@app.route('/',methods=['GET','POST'])
def select_insert():
    global PARAMS
    if request.method == 'POST':
        if 'plasmid' not in request.files or request.files['plasmid'].filename == '':
            flash('Missing Input File')
            return render_template('webapp.html')


        seqFile = secure_filename(request.files['plasmid'].filename)
        if '.' not in seqFile or seqFile.split('.')[-1].lower() != 'gb':
            flash('Invalid Input File')
            return render_template('webapp.html')
        request.files['plasmid'].save(os.path.join(UPLOAD_DIR,seqFile))
        #circular = True if request.form['circular'] == 'True' else False
        #results = {'identity': 5, 'orientation': 10}
        
        data = parse(seqFile)
        if len(data) == 4:
            results = {'primers': [data[0], data[1], data[2]], 'attributes': data[3]}
        else:
            results = {'primers': ["N/A", "N/A", "N/A"], 'attributes': "N/A"}
        return render_template('results.html', results=results)

    return render_template('webapp.html')

app.run(debug=True, port=5000)
