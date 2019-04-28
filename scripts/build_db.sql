
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
       "hg-panTro"    TEXT
);

.import /dev/stdin snp
.exit
