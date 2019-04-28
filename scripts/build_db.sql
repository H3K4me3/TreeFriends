
.separator "\t"

DROP TABLE if EXISTS snp;

CREATE TABLE snp (
       "seqnames"   TEXT,
       "start"      INT,
       "end"        INT,
       "strand"     TEXT,
       "REF"        TEXT,
       "ALL"        TEXT,
       "NS"         INT,
       "AN"         INT,
       "AC_A"       INT,
       "AC_C"       INT,
       "AC_G"       INT,
       "AC_T"       INT,
       "hg"         TEXT,
       "panTro"     TEXT,
       "gorGor"     TEXT,
       "ponAbe"     TEXT,
       "rheMac"     TEXT,
       "root"       TEXT,
       "ponAbe-else"  TEXT,
       "hg-pan-gor"   TEXT,
       "hg-panTro"    TEXT,
       /* "end" should be same as "start" */
       PRIMARY KEY (seqnames, start)
);

/* I will assume that index have been created for the primary key.
 * So I will not create an index on SNP location.                  */

.import /dev/stdin snp

.exit
