
using SQLite
using JuliaDB

const DB_FILE = "results/snpdb.sqlite3"

if !isfile(DB_FILE)
    throw("Database not available")
end

const db_conn = SQLite.DB(DB_FILE)

SQLite.tables(db_conn)


## Function to calculate statistics
function has_suspicious_AC(AC_A::Union{Missing, Int},
                           AC_C::Union{Missing, Int},
                           AC_G::Union{Missing, Int},
                           AC_T::Union{Missing, Int})::Bool
    AN = 125568
    convert_missing(x)::Int = if ismissing(x) return 0 else x end
    sum = convert_missing(AC_A) + convert_missing(AC_C) +
          convert_missing(AC_G) + convert_missing(AC_T)
    sum >= AN
end

SQLite.register(db_conn, has_suspicious_AC)

SQLite.execute!(
    db_conn,
    """
    CREATE TEMP VIEW snp_stat AS
    SELECT
        seqnames,
        start,
        REF,
        "hg",
        "panTro",
        "gorGor",
        "ponAbe",
        "rheMac",

        "root",
        "ponAbe-else",
        "hg-pan-gor",
        "hg-panTro",

        has_suspicious_AC(AC_A, AC_C, AC_G, AC_T) AS suspicious_ac
    FROM snp
    """
)

SQLite.tables(db_conn)

t = collect(SQLite.Query(db_conn, "select * from snp_stat where suspicious_AC > 0")) |> DataFrame
t

