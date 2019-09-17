
import os
import sys
import subprocess

SCRIPT_DIR = os.path.dirname(os.path.realpath(__file__))
PROJECT_ROOT = os.path.join(SCRIPT_DIR, "..")
CSV_FILE = os.path.join(PROJECT_ROOT, "results", "rsid_position.csv")

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

if __name__ == "__main__":
    download_csv()
