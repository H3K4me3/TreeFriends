
import os
import re
import sys
from flask import Flask
from flask import jsonify
from flask import request
from flask import render_template

#### Import lib/snpdb.py

if "__file__" in globals():
    PROJECT_ROOT = os.path.join(os.path.dirname(__file__), "..")
else:
    PROJECT_ROOT = ".."
sys.path.insert(0, os.path.join(PROJECT_ROOT, "lib"))

from snpdb import SNPDB
import ParsimonyInfer

DB_PATH = os.path.join(PROJECT_ROOT, "results/snpdb.sqlite3")
db = SNPDB(DB_PATH)

app = Flask(__name__,
            static_url_path = "",
            static_folder = os.path.join(PROJECT_ROOT, "webapp/static"))

@app.route('/')
def application():
    chromosome = request.args.get('chromosome')
    position = request.args.get('position')
    rsid = request.args.get("rsid")
    if rsid is not None:
        rsid = int(rsid)
        data = get_snp_byrsnum(rsid)
        if not data['found']:
            return render_template("application.html", rsid_num=rsid, chromosome=None, position=None, data=data)
        return render_template("application.html", rsid_num=rsid, chromosome=data['results']['chromosome'],
            position=data['results']['position'], data=data)
    chromosome = chromosome if chromosome is None else str(chromosome)
    position = position if position is None else int(position)

    return render_template("application.html",
        rsid_num=None, chromosome=chromosome, position=position, data=get_snp(chromosome, position))

## A restful API
@app.route('/snp/<chromosome>/<int:position>')
def api_get_snp(chromosome, position):
    return jsonify(get_snp(chromosome, position))

def get_snp(chromosome, position):
    res = db.get_snp(chromosome, position)
    if res is None:
        return res
    assert isinstance(res, dict)
    edge_changes = ParsimonyInfer.stat_edge_changes(ParsimonyInfer.mkNodeTuple(res))
    edge_changes = edge_changes._asdict()
    merged = dict()
    merged['results'] = res
    merged['edge_changes'] = edge_changes
    return merged

@app.route("/rsid/<int:rsid_num>")
def api_get_snp_byrsnum(rsid_num):
    return jsonify(get_snp_byrsnum(rsid_num))

def get_snp_byrsnum(rsid_num):
    assert isinstance(rsid_num, int)
    res = list(db.get_rsidnum_batch([rsid_num]))
    assert len(res) == 1
    res = res[0]
    res = dict(res)
    merged = dict()
    merged['rsid_num'] = rsid_num
    merged['found'] = False
    ## rsid not found
    if res['chromosome'] is None and res['position'] is None:
        return merged
    merged['found'] = True
    edge_changes = ParsimonyInfer.stat_edge_changes(ParsimonyInfer.mkNodeTuple(res))
    edge_changes = edge_changes._asdict()
    merged['results'] = res
    merged['edge_changes'] = edge_changes
    return merged
