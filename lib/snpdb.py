import sqlite3

class SNPDB:
    conn = None
    def __init__(self, dbpath):
        self.conn = sqlite3.connect(dbpath)
        self.register_expected_REF_AC()
        self.createview_suspicious_stat()
    def register_pyfunc(self, name, nargs, func):
        self.conn.create_function(name, nargs, func)
        return self
    def register_expected_REF_AC(self):
        conn = self.conn
        cursor = conn.cursor()
        cursor.execute("select AN from snp limit 1")
        AN = cursor.fetchone()[0]
        if AN != 125568:
            raise Exception("AN value is not expected: {}".format(AN))
        cursor.close()
        def expected_REF_AC(AC_A, AC_C, AC_G, AC_T):
            if AC_A is None:
                AC_A = 0
            if AC_C is None:
                AC_C = 0
            if AC_G is None:
                AC_G = 0
            if AC_T is None:
                AC_T = 0
            return AN - (AC_A + AC_C + AC_G + AC_T)
        self.register_pyfunc("expected_REF_AC", 4, expected_REF_AC)
        return self
    def createview_suspicious_stat(self):
        conn = self.conn
        conn.execute(
            """
            CREATE TEMP VIEW suspicious_stat AS
            SELECT
                seqnames,
                start,
                REF,
                "hg",
                AC_A,
                AC_C,
                AC_G,
                AC_T,

                expected_REF_AC(AC_A, AC_C, AC_G, AC_T) AS REF_AC
            FROM snp
            """
        )
    def query(self, string):
        cursor = self.conn.cursor()
        cursor.execute(string)
        return cursor
