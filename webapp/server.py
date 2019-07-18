
import os
import sys
from flask import Flask
from flask import jsonify

#### Import lib/snpdb.py

if "__file__" in globals():
    PROJECT_ROOT = os.path.join(os.path.dirname(__file__), "..")
else:
    PROJECT_ROOT = ".."
sys.path.insert(0, os.path.join(PROJECT_ROOT, "lib"))

from snpdb import SNPDB

DB_PATH = os.path.join(PROJECT_ROOT, "results/snpdb.sqlite3")
db = SNPDB(DB_PATH)

app = Flask(__name__,
            static_url_path = "",
            static_folder = os.path.join(PROJECT_ROOT, "webapp/static"))

@app.route('/')
def root():
     return app.send_static_file('index.html')

@app.route('/snp/<chromosome>/<int:position>')
def get_snp(chromosome, position):
    res = db.get_snp(chromosome, position)
    return jsonify(res)


