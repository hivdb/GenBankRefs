import pandas as pd
import sqlite3
from functools import reduce


def create_tables(db_file):

    conn = sqlite3.connect(db_file)

    conn.execute("""
        CREATE TABLE "tblGBRefs" (
            "RefID" INTEGER PRIMARY KEY,
            "Authors" TEXT,
            "Title" TEXT,
            "Journal" TEXT,
            "PMID" TEXT,
            "Year" TEXT
        )
    """)

    conn.execute("""
        CREATE TABLE "tblIsolates" (
            "Accession" TEXT PRIMARY KEY,
            "Country" TEXT,
            "RecordYear" INTEGER,
            "IsolateYear" INTEGER,
            "Host" TEXT,
            "Specimen" TEXT,
            "IsolateName" TEXT,
            "SeqLength" INTEGER,
            "IsolateType" TEXT
        )
    """)

    conn.execute("""
        CREATE TABLE "tblGBRefLink" (
            "RefID" INTEGER,
            "Accession" TEXT,
            FOREIGN KEY (RefID) REFERENCES tblGBRefs (RefID),
            FOREIGN KEY (Accession) REFERENCES tblIsolates (Accession)
        )
    """)

    conn.execute("""
        CREATE TABLE "tblSequences" (
            "SeqID" INTEGER PRIMARY KEY,
            "Accession" TEXT,
            "Gene" TEXT,
            "CDS_NAME" TEXT,
            "AA_seq" TEXT,
            "AA_length" INTEGER,
            "AA_start" INTEGER,
            "AA_stop" INTEGER,
            "NA_Seq" TEXT,
            "NA_length" INTEGER,
            "NA_start" INTEGER,
            "NA_stop" INTEGER,
            "PcntMatch" REAL,
            "HSPLength" INTEGER,
            FOREIGN KEY (Accession) REFERENCES tblSequences (Accession)
        )
    """)

    conn.execute("""
        CREATE TABLE "tblSequenceQA" (
            "SeqID" INTEGER,
            "AA_num_ins" INTEGER,
            "AA_num_del" INTEGER,
            "AA_blast_failed" INTEGER,
            "NA_num_ins" INTEGER,
            "NA_num_del" INTEGER,
            "NA_blast_failed" INTEGER,
            "num_N" INTEGER,
            "translation_issue" INTEGER,
            FOREIGN KEY (SeqID) REFERENCES tblSequences (SeqID)
        )
    """)

    conn.execute("""
        CREATE TABLE "tblPublications" (
            "PubID" INTEGER PRIMARY KEY,
            "Authors" TEXT,
            "Title" TEXT,
            "Journal" TEXT,
            "PMID" TEXT,
            "Year" TEXT
        )
    """)

    conn.execute("""
        CREATE TABLE "tblPublicationData" (
            "PubID" INTEGER,
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
            FOREIGN KEY (PubID) REFERENCES tblPublications (PubID)
        )
    """)

    conn.execute("""
        CREATE TABLE "tblGBPubRefLink" (
            "PubID" INTEGER,
            "RefID" TEXT,
            FOREIGN KEY (PubID) REFERENCES tblPublications (PubID),
            FOREIGN KEY (RefID) REFERENCES tblGBRefs (RefID)
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
    tblGBRefss = references[[
        'RefID', 'Authors', 'Title', 'Journal', 'PMID', 'Year']]
    dump_table(virus_obj.DB_FILE,
               'tblGBRefs', tblGBRefss)

    create_ref_link(virus_obj, references)

    features['Specimen'] = features['isolate_source']
    features['Virus'] = features['organism']

    if 'IsolateType' not in features.columns:
        features['IsolateType'] = ""

    tblIsolates = features[[
        'Accession', 'Country', 'RecordYear',
        'IsolateYear', 'Host', 'Specimen', 'IsolateName',
        'SeqLength',
        'IsolateType']]

    for i, row in tblIsolates.iterrows():
        tblIsolates.at[i, 'Host'] = row['Host'] if row['Host'] else (
            'Not applicable' if row['IsolateType'] else 'Not available')
        tblIsolates.at[i, 'Specimen'] = row['Specimen'] if row['Specimen'] else (
            'Not applicable' if row['IsolateType'] else 'Not available')
        tblIsolates.at[i, 'IsolateYear'] = row['IsolateYear'] if row['IsolateYear'] else (
            'Not applicable' if row['IsolateType'] else 'Not available')
        tblIsolates.at[i, 'Country'] = row['Country'] if row['Country'] else (
            'Not applicable' if row['IsolateType'] else 'Not available')

    dump_table(
        virus_obj.DB_FILE,
        'tblIsolates', tblIsolates)

    genes['PcntMatch'] = genes['pcnt_id']

    tblGBSequences = genes[[
        'Accession', 'Gene', 'CDS_NAME',
        'AA_seq', 'AA_length', 'AA_start', 'AA_stop',
        'NA_seq', 'NA_length', 'NA_start', 'NA_stop',
        'PcntMatch',
    ]]
    dump_table(
        virus_obj.DB_FILE,
        'tblSequences',
        tblGBSequences)

    tblSequenceQA = genes[
        (genes['AA_num_ins'] != 0) |
        (genes['AA_num_del'] != 0) |
        (genes['AA_blast_failed'] == 1) |
        (genes['NA_num_ins'] != 0) |
        (genes['NA_num_del'] != 0) |
        (genes['NA_blast_failed'] == 1) |
        (genes['num_N'] != 0) |
        (genes['translation_issue'] != 0)
    ]

    tblSequenceQA = tblSequenceQA[[
        'SeqID',
        'AA_num_ins', 'AA_num_del', 'AA_blast_failed',
        'NA_num_ins', 'NA_num_del', 'NA_blast_failed',
        'num_N', 'translation_issue'
    ]]

    dump_table(
        virus_obj.DB_FILE,
        'tblSequenceQA',
        tblSequenceQA)

    # PubMed Tables
    tblPublications = pubmed[[
        'PubID',
        'Authors',
        'Title',
        'Journal',
        'PMID',
        'Year'
    ]]
    dump_table(
        virus_obj.DB_FILE,
        'tblPublications', tblPublications)

    tblPublicationData = pubmed[[
        'PubID',
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
        'tblPublicationData',
        tblPublicationData)

    tblPubRefLink = []
    for pubmed, genbank_list in pubmed_genbank:
        for g in genbank_list:
            tblPubRefLink.append((pubmed['PubID'], g['RefID']))

    tblPubRefLink = list(set(tblPubRefLink))
    tblPubRefLink = [
        {
            'PubID': i,
            'RefID': j
        }
        for i, j in tblPubRefLink
    ]

    dump_table(
        virus_obj.DB_FILE,
        'tblGBPubRefLink',
        pd.DataFrame(tblPubRefLink))

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
        'tblGBRefLink',
        pd.DataFrame(ref_link))


def creat_views(db_file):

    vReferenceRecord = """
        CREATE VIEW vGBRecords AS
        SELECT a.*, c.*
        FROM tblGBRefs a, tblGBRefLink b, tblIsolates c
        WHERE a.RefID = b.RefID
        AND b.Accession = c.Accession;
    """
    run_create_view(db_file, vReferenceRecord)

    vPubAccessionLink = """
        CREATE VIEW vPubAccessionLink AS
        SELECT
            d.*,
            a.*
        FROM
            tblIsolates a,
            tblGBRefLink b,
            tblGBPubRefLink c,
            tblPublications d
        WHERE
            a.Accession = b.Accession
            AND
            b.RefID = c.RefID
            AND
            c.PubID = d.PubID
        ;
    """
    run_create_view(db_file, vPubAccessionLink)

    vPubDataAccessionLink = """
        CREATE VIEW vPubDataAccessionLink AS
        SELECT
            d.*,
            a.*
        FROM
            tblIsolates a,
            tblGBRefLink b,
            tblGBPubRefLink c,
            (SELECT * FROM tblPublications i LEFT JOIN tblPublicationData j ON i.PubID = j.PubID) d
        WHERE
            a.Accession = b.Accession
            AND
            b.RefID = c.RefID
            AND
            c.PubID = d.PubID
        ;
    """
    run_create_view(db_file, vPubDataAccessionLink)

    vNotClinicalAccession = """
    SELECT      
        Comment,      
        COUNT(*) AS count,     
        SUM(COUNT(*)) OVER () AS total_count 
    FROM 
        tblIsolates 
    WHERE 
        Comment IS NOT "" 
    GROUP BY 
        Comment;"""
    
    run_create_view(db_file, vNotClinicalAccession)

    vIsolateHost = """
    SELECT 
        Host, COUNT(*) AS count 
    FROM 
        tblIsolates 
    GROUP BY 
        Host
    """
    run_create_view(db_file, vIsolateHost)

    vIsolateSpecimen = """
    SELECT 
        Specimen, COUNT(*) AS count 
    FROM 
        tblIsolates 
    GROUP BY 
        Specimen
    """
    run_create_view(db_file, vIsolateSpecimen)
    
    vIsolateCountry = """
    SELECT 
        Country, COUNT(*) AS count 
    FROM 
        tblIsolates 
    GROUP BY 
        Country
    """
    run_create_view(db_file, vIsolateCountry)
    
    vIsolateYear = """
    SELECT 
        IsolateYear, COUNT(*) AS count 
    FROM 
        tblIsolates 
    GROUP BY 
        IsolateYear
    """
    run_create_view(db_file, vIsolateYear)
    
    vMissingData = """
    SELECT 
    (SELECT COUNT(*) 
     FROM tblIsolates 
     WHERE Host IN ("Not Applicable", "Not Available", "Not available", "Not applicable")) AS host_count,
    (SELECT COUNT(*) 
     FROM tblIsolates 
     WHERE Specimen IN ("Not Applicable", "Not Available", "Not available", "Not applicable")) AS specimen_count,
    (SELECT COUNT(*) 
     FROM tblIsolates 
     WHERE Country IN ("Not Applicable", "Not Available", "Not available", "Not applicable")) AS country_count,
    (SELECT COUNT(*) 
     FROM tblIsolates 
     WHERE IsolateYear IN ("Not Applicable", "Not Available", "Not available", "Not applicable")) AS isolate_yr_count;
    """

    run_create_view(db_file, vMissingData)

    vSubmissionSetsSupplementingGB = """
    SELECT 
        COUNT(DISTINCT CASE WHEN tblIsolates.Host LIKE '%*%' THEN tblGBRefLink.RefID END) AS submission_hosts,
        COUNT(DISTINCT CASE WHEN tblIsolates.Host LIKE '%*%' THEN tblIsolates.Accession END) AS numseq_hosts,
        COUNT(DISTINCT CASE WHEN tblIsolates.Specimen LIKE '%*%' THEN tblGBRefLink.RefID END) AS submission_specimens,
        COUNT(DISTINCT CASE WHEN tblIsolates.Specimen LIKE '%*%' THEN tblIsolates.Accession END) AS numseq_specimens,
        COUNT(DISTINCT CASE WHEN tblIsolates.Country LIKE '%*%' THEN tblGBRefLink.RefID END) AS submission_country,
        COUNT(DISTINCT CASE WHEN tblIsolates.Country LIKE '%*%' THEN tblIsolates.Accession END) AS numseq_country,
        COUNT(DISTINCT CASE WHEN tblIsolates.IsolateYear LIKE '%*%' THEN tblGBRefLink.RefID END) AS submission_IsolateYear,
        COUNT(DISTINCT CASE WHEN tblIsolates.IsolateYear LIKE '%*%' THEN tblIsolates.Accession END) AS numseq_IsolateYear
    FROM 
        tblGBRefLink
    JOIN 
        tblIsolates ON tblGBRefLink.Accession = tblIsolates.Accession;

    """

    run_create_view(db_file, vSubmissionSetsSupplementingGB)

    vChordTable = """
    SELECT DISTINCT 
        tblIsolates.Host, 
        tblIsolates.Country, 
        tblIsolates.IsolateYear, 
        tblSequences.Gene,
        COUNT(*) AS count_rows ,
        SUM(COUNT(*)) OVER () AS Total_count,
        ROUND((COUNT(*) * 100.0 / SUM(COUNT(*)) OVER ()), 2) AS "%"
    FROM 
        tblIsolates
    JOIN 
        tblSequences ON tblIsolates.Accession = tblSequences.Accession
    GROUP BY 
        tblIsolates.Host, tblIsolates.Country, tblIsolates.IsolateYear, tblSequences.Gene 
    ORDER BY 
        count_rows DESC, tblIsolates.Host, tblIsolates.Country, tblSequences.Gene;
    """
    # ROUND((COUNT(*) * 100.0 / SUM(COUNT(*)) OVER ()), 2) AS "%" takes too long to calculate
    run_create_view(db_file, vChordTable)

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
