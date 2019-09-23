import sqlite3

class SNPDB:
    conn = None
    def __init__(self, dbpath):
        ## Use check_same_thread=False -- we do not have writing operations
        self.conn = sqlite3.connect(dbpath, check_same_thread=False)
        self.register_expected_REF_AC()
        self.createview_suspicious_stat()
        self.createview_tree_nodes()
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
    ## only available for this session
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
    def createview_tree_nodes(self):
        conn = self.conn
        conn.execute(
            """
            CREATE TEMP VIEW tree_nodes AS
            SELECT
              seqnames as chromosome,
              start as position,
              REF as ref,
              NS as sample_number,
              AN as allele_number,
              hg, panTro, gorGor, ponAbe, rheMac,
              `hg-panTro` as ancestry1,
              `hg-pan-gor` as ancestry2,
              `ponAbe-else` as ancestry3,
              root as ancestry4
            FROM snp
            """
        )


    def query(self, string, binding = None):
        cursor = self.conn.cursor()
        if binding is None:
            cursor.execute(string)
        else:
            cursor.execute(string, binding)
        return cursor

    def get_available_chromosomes(self):
        cursor = self.query(
            """
            select distinct seqnames from snp
            """
        )
        ans = []
        for row in cursor:
            ans.append(row[0])
        return ans

    def get_snp(self, chromosome, loc):
        cursor = self.query(
            """
            select * from tree_nodes where chromosome=? and position=?
            """,
            (chromosome, loc)
        )
        ## Convert to dict
        cursor.row_factory = sqlite3.Row
        res = cursor.fetchone()
        if res is None:
            return None
        return dict(res)
    
    def get_snp_batch(self, coordinates):
        cursor = self.conn.cursor()
        cursor.row_factory = sqlite3.Row
        for chromosome, position in coordinates:
            cursor.execute(
                """
                WITH loc(chromosome, position) AS (VALUES(?, ?))
                SELECT * FROM loc LEFT JOIN tree_nodes ON
                    loc.chromosome=tree_nodes.chromosome AND loc.position=tree_nodes.position
                """,
                (chromosome, position)
            )
            res = cursor.fetchone()
            if cursor.fetchone() is not None:
                raise Exception("Unexpected: there should be only one row")
            yield dict(res)
    
    def get_rsidnum_batch(self, rsidnum):
        allowed_chromosome = ["chr" + str(i) for i in range(1,23)]
        allowed_chromosome.append("chrX")
        cursor = self.conn.cursor()
        cursor.row_factory = sqlite3.Row
        sql = """
        WITH
            query_rsid_num(rsid_num) AS (VALUES(?)),
            rsid_position_res AS (
                SELECT
                    printf('rs%d', query_rsid_num.rsid_num) AS rsid,
                    seqnames AS chromosome,
                    start AS position
                FROM
                    query_rsid_num
                LEFT JOIN
                    rsid_position
                ON
                    query_rsid_num.rsid_num=rsid_position.rsid_num
            ),
            final_res AS (
                SELECT *
                FROM
                    rsid_position_res
                LEFT JOIN
                    tree_nodes
                ON
                    rsid_position_res.chromosome=tree_nodes.chromosome AND
                    rsid_position_res.position=tree_nodes.position
            )

        SELECT * FROM final_res WHERE chromosome IN ({}) OR chromosome IS NULL
        """.format(",".join(map(lambda x: "'{}'".format(x), allowed_chromosome)))
        for id in rsidnum:
            cursor.execute(sql, (id,))
            res = cursor.fetchall()
            if len(res) == 1:
                yield dict(res[0])
            else:
                res = [ each for each in res if each['ref'] is not None ]
                if len(res) != 1:
                    raise Exception("Something unexpected happened with rs{}".format(id))
                yield dict(res[0])
            
        

