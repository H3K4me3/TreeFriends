
import os
import sys
import subprocess
import sqlite3
from datetime import datetime

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

def csv_iter():
    with open(CSV_FILE, "r") as csv_file:
        line_count = 0
        for line in csv_file:
            seqnames, start, rsid = line.rstrip().split(",")
            rsid_num = int(rsid[2:])
            start = int(start)
            if line_count % 1000000 == 0:
                print("{}\trs{}\t{}\t{}\t{}".format(line_count, rsid_num, seqnames, start, datetime.now().ctime()))
            line_count += 1
            yield rsid_num, seqnames, start

def insert_csv_to_db():
    conn = sqlite3.connect(DB_PATH)
    conn.execute(
        " drop index if exists idx_rsid_position_on_rsid_num "
    )
    conn.execute(
        " drop table if exists rsid_position "
    )
    conn.execute(
        " create table rsid_position(rsid_num INT, seqnames TEXT, start INT) "
    )
    conn.executemany(
        " INSERT INTO rsid_position(rsid_num, seqnames, start) VALUES (?, ?, ?) ", csv_iter()
    )
    conn.execute(
        # Non-unique index
        " CREATE INDEX idx_rsid_position_on_rsid_num ON rsid_position(rsid_num) "
    )
    conn.commit()
    conn.close()

if __name__ == "__main__":
    download_csv() 
    insert_csv_to_db()
