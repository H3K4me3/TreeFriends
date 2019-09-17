
import os
import sys
import subprocess
import sqlite3

SCRIPT_DIR = os.path.dirname(os.path.realpath(__file__))
PROJECT_ROOT = os.path.join(SCRIPT_DIR, "..")
CSV_FILE = os.path.join(PROJECT_ROOT, "results", "rsid_position.csv")
DB_PATH = os.path.join(PROJECT_ROOT, "results/snpdb.sqlite3")

SH_CMD = """
    curl http://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/snp151.txt.gz | gunzip -c | awk '{if ($12 == "single") print $2 "," $4 "," $5 }'
"""

def download_csv():
    if os.path.exists(CSV_FILE):
        os.remove(CSV_FILE)
    csv_file = open(CSV_FILE, "w")
    proc = subprocess.Popen(SH_CMD, shell=True, stdout=csv_file, stderr=sys.stderr)
    proc.wait()
    if proc.returncode > 0:
        print("===== Download process failed with code {} =====".format(proc.returncode))
        os.remove(CSV_FILE)
    else:
        print("===== File downloaded to {} =====".format(CSV_FILE))

sql = """
    drop table if exists rsid_position;
    create table rsid_position(seqnames TEXT, start INT, rsid TEXT);
    .mode csv
    .import rsid_position.csv rsid_position
"""

def csv_iter():
    with open(CSV_FILE, "r") as csv_file:
        for line in csv_file:
            seqnames, start, rsid = line.rstrip().split(",")
            start = int(start)
            yield seqnames, start, rsid

def insert_csv_to_db():
    conn = sqlite3.connect(DB_PATH)
    conn.execute(
        " drop index if exists rsid_position_idx_byposition "
    )
    conn.execute(
        " drop index if exists rsid_position_idx_byrsid "
    )
    conn.execute(
        " drop table if exists rsid_position "
    )
    conn.execute(
        " create table rsid_position(seqnames TEXT, start INT, rsid TEXT) "
    )
    conn.executemany(
        " INSERT INTO rsid_position(seqnames, start, rsid) VALUES (?, ?, ?) ", csv_iter()
    )
    conn.execute(
        " CREATE UNIQUE INDEX rsid_position_idx_byposition ON rsid_position(seqnames, start) "
    )
    conn.execute(
        " CREATE UNIQUE INDEX rsid_position_idx_byrsid ON rsid_position(rsid) "
    )
    conn.commit()
    conn.close()

if __name__ == "__main__":
    download_csv() 
    insert_csv_to_db()
