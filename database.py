import pandas as pd
import sqlite3
from functools import reduce


def create_tables(db_file):

    conn = sqlite3.connect(db_file)

    conn.execute("""
        CREATE TABLE "tblGBReferences" (
            "RefID" INTEGER PRIMARY KEY,
            "Authors" TEXT,
            "Title" TEXT,
            "Journal" TEXT,
            "PMID" TEXT,
            "Year" TEXT
        )
    """)

    conn.execute("""
        CREATE TABLE "tblGBIsolates" (
            "Accession" TEXT PRIMARY KEY,
            "Country" TEXT,
            "RecordYear" INTEGER,
            "IsolateYear" INTEGER,
            "Host" TEXT,
            "Specimen" TEXT,
            "IsolateName" TEXT,
            "Comment" TEXT
        )
    """)

    conn.execute("""
        CREATE TABLE "tblGBRefSeqLink" (
            "RefID" INTEGER,
            "Accession" TEXT,
            FOREIGN KEY (RefID) REFERENCES tblGBReferences (RefID),
            FOREIGN KEY (Accession) REFERENCES tblGBIsolates (Accession)
        )
    """)

    conn.execute("""
        CREATE TABLE "tblSequences" (
            "Accession" TEXT,
            "Gene" TEXT,
            "CDS_NAME" TEXT,
            "AASeq" TEXT,
            "NumAA" INTEGER,
            "AA_start" INTEGER,
            "AA_stop" INTEGER,
            "NASeq" TEXT,
            "NumNA" INTEGER,
            "NA_start" INTEGER,
            "NA_stop" INTEGER,
            "PcntMatch" REAL,
            "HSPLength" INTEGER,
            FOREIGN KEY (Accession) REFERENCES tblGBIsolates (Accession)
        )
    """)

    conn.execute("""
        CREATE TABLE "tblLiteratures" (
            "LitID" INTEGER PRIMARY KEY,
            "Authors" TEXT,
            "Title" TEXT,
            "Journal" TEXT,
            "PMID" TEXT,
            "Year" TEXT
        )
    """)

    conn.execute("""
        CREATE TABLE "tblLitTextReview" (
            "LitID" INTEGER,
            "Viruses" TEXT,
            "NumSeqs" TEXT,
            "Host" TEXT,
            "SampleYr" TEXT,
            "Country" TEXT,
            "GenBank" TEXT,
            "SeqMethod" TEXT,
            "CloneMethod" TEXT,
            "IsolateType" TEXT,
            "Gene" TEXT,
            FOREIGN KEY (LitID) REFERENCES tblLiteratures (LitID)
        )
    """)

    conn.execute("""
        CREATE TABLE "tblLitGBRefLink" (
            "LitID" INTEGER,
            "RefID" TEXT,
            FOREIGN KEY (LitID) REFERENCES tblLiteratures (LitID),
            FOREIGN KEY (RefID) REFERENCES tblGBReferences (RefID)
        )
    """)

    conn.close()


def create_database(
        virus_obj,
        references, features, genes,
        pubmed, pubmed_genbank):

    virus_obj.DB_FILE.unlink(missing_ok=True)

    create_tables(virus_obj.DB_FILE)

    # GenBank Tables
    tblGBReferencess = references[[
        'RefID', 'Authors', 'Title', 'Journal', 'PMID', 'Year']]
    dump_table(virus_obj.DB_FILE,
               'tblGBReferences', tblGBReferencess)

    create_ref_link(virus_obj, references)

    features['Specimen'] = features['isolate_source']
    features['Virus'] = features['organism']

    if 'Comment' not in features.columns:
        features['Comment'] = ""

    tblGBIsolates = features[[
        'Accession', 'Country', 'RecordYear',
        'IsolateYear', 'Host', 'Specimen', 'IsolateName',
        'Comment']]

    for i, row in tblGBIsolates.iterrows():
        tblGBIsolates.at[i, 'Host'] = row['Host'] if row['Host'] else ('Not applicable' if row['Comment'] else 'Not available')
        tblGBIsolates.at[i, 'Specimen'] = row['Specimen'] if row['Specimen'] else ('Not applicable' if row['Comment'] else 'Not available')
        tblGBIsolates.at[i, 'IsolateYear'] = row['IsolateYear'] if row['IsolateYear'] else ('Not applicable' if row['Comment'] else 'Not available')
        tblGBIsolates.at[i, 'Country'] = row['Country'] if row['Country'] else ('Not applicable' if row['Comment'] else 'Not available')


    dump_table(
        virus_obj.DB_FILE,
        'tblGBIsolates', tblGBIsolates)

    genes['PcntMatch'] = genes['pcnt_id']
    genes['HSPLength'] = genes['align_len']

    tblGBSequences = genes[[
        'Accession', 'Gene', 'CDS_NAME',
        'AASeq', 'NumAA', 'AA_start', 'AA_stop',
        'NASeq', 'NumNA', 'NA_start', 'NA_stop',
        'PcntMatch', 'HSPLength',
    ]]
    dump_table(
        virus_obj.DB_FILE,
        'tblSequences',
        tblGBSequences)

    # PubMed Tables
    tblLiteratures = pubmed[[
        'LitID',
        'Authors',
        'Title',
        'Journal',
        'PMID',
        'Year'
    ]]
    dump_table(
        virus_obj.DB_FILE,
        'tblLiteratures', tblLiteratures)

    tblLitTextReview = pubmed[[
        'LitID',
        'Viruses',
        'NumSeqs',
        'Host',
        'SampleYr',
        'Country',
        'GenBank',
        'SeqMethod',
        'CloneMethod',
        'IsolateType',
        'Gene'
    ]]
    dump_table(
        virus_obj.DB_FILE,
        'tblLitTextReview',
        tblLitTextReview)

    tblLitRefLink = []
    for pubmed, genbank_list in pubmed_genbank:
        for g in genbank_list:
            tblLitRefLink.append((pubmed['LitID'], g['RefID']))

    tblLitRefLink = list(set(tblLitRefLink))
    tblLitRefLink = [
        {
            'LitID': i,
            'RefID': j
        }
        for i, j in tblLitRefLink
    ]

    dump_table(
        virus_obj.DB_FILE,
        'tblLitGBRefLink',
        pd.DataFrame(tblLitRefLink))

    creat_views(virus_obj.DB_FILE)

    # get_table_schema_sql(virus_obj.DB_FILE)


def create_ref_link(virus_obj, ref):
    ref_link = []
    for i, row in ref.iterrows():
        accessions = row['accession']
        accessions = [
            i.strip() for i in accessions.split(',') if i.strip()]
        for acc in accessions:
            ref_link.append({
                'RefID': row['RefID'],
                'Accession': acc
            })

    dump_table(
        virus_obj.DB_FILE,
        'tblGBRefSeqLink',
        pd.DataFrame(ref_link))


def creat_views(db_file):

    vReferenceRecord = """
        CREATE VIEW vGBRecords AS
        SELECT a.*, c.*
        FROM tblGBReferences a, tblGBRefSeqLink b, tblGBIsolates c
        WHERE a.RefID = b.RefID
        AND b.Accession = c.Accession;
    """
    run_create_view(db_file, vReferenceRecord)

    vLitAccessionLink = """
        CREATE VIEW vLitAccessionLink AS
        SELECT
            d.*,
            a.*
        FROM
            tblGBIsolates a,
            tblGBRefSeqLink b,
            tblLitGBRefLink c,
            tblLiteratures d
        WHERE
            a.Accession = b.Accession
            AND
            b.RefID = c.RefID
            AND
            c.LitID = d.LitID
        ;
    """
    run_create_view(db_file, vLitAccessionLink)


def dump_table(db_file, table_name, table):
    conn = sqlite3.connect(str(db_file))

    table.to_sql(
        table_name,
        conn,
        if_exists='append',
        index=False)


def run_create_view(db_file, sql):
    conn = sqlite3.connect(str(db_file))

    cursor = conn.cursor()

    cursor.execute(sql)
    conn.commit()
    conn.close()


def load_table(db_file, table_name):
    conn = sqlite3.connect(str(db_file))
    return pd.read_sql(f'SELECT * FROM {table_name};', conn)


def load_tables(db_file, table_name_list):
    tables = []
    for i in table_name_list:
        tables.append(load_table(db_file, i))

    return tables


def split_table(table_config, dataframe):
    tables = []
    for table_name, columns in table_config:
        table = dataframe[columns]
        tables.append((table_name, table))

    return tables


def merge_table(tables, key):
    return reduce(lambda x, y: pd.merge(x, y, on=key, how='inner'), tables)


def get_table_schema_sql(db_file):

    conn = sqlite3.connect(db_file)

    cursor = conn.execute(f"SELECT name FROM sqlite_master WHERE type='table' AND name NOT LIKE 'sqlite_%';")

    for i in cursor.fetchall():
        table_name = i[0]

        cursor = conn.execute(f"SELECT sql FROM sqlite_master WHERE type='table' AND name='{table_name}';")

        schema = cursor.fetchone()[0]

        print(schema)

    conn.close()
