#!/usr/bin/env python3

import os
import sys
import sqlite3
import argparse

SCRIPT_DIR = os.path.dirname(__file__)
DB_PATH = os.path.join(SCRIPT_DIR, "../results/snpdb.sqlite3")

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("chr")
    parser.add_argument("start", type = int)
    parser.add_argument("end",   type = int)
    args = parser.parse_args()
    return args.chr, args.start, args.end

class SNPDB:
    """
    Interface for manipulating the database
    """
    conn = None

    def __init__(self):
        self.load_db()

    def load_db(self):
        if self.conn is None:
            self.conn = sqlite3.connect(DB_PATH)
            #self.conn.row_factory = sqlite3.Row
        return True

    def subset_range(self, seqname, lower, upper):
        c = self.conn

        sql = """
        select * from snp where
        seqnames == ? and start >= ? and start <= ?
        """
        cursor = c.execute(sql, (seqname, lower, upper))
        return cursor

    @staticmethod
    def format_res(cursor, format = "tsv", out = sys.stdout):
        if format != "tsv":
            raise NotImplementedError("Not implemented for format {}".format(format))
        headers = [x[0] for x in cursor.description]
        out.write("\t".join(headers))
        out.write("\n")
        for row in cursor:
            out.write("\t".join((str(x) for x in row)))
            out.write("\n")
        out.flush()
        return True

def main():
    if not os.path.isfile(DB_PATH):
        sys.exit("The DB file does not exist")
    db = SNPDB()
    chr, start, end = parse_args()
    cursor = db.subset_range(chr, start, end)
    SNPDB.format_res(cursor)


if __name__ == "__main__":
    main()

# Local Variables:
# compile-command: "./snpdb.py chr10 1 11000"
# End:
