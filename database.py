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

    vGBRefIsolates = """
        CREATE VIEW vGBRefIsolates AS
        SELECT a.*, c.*
        FROM tblGBRefs a, tblGBRefLink b, tblIsolates c
        WHERE a.RefID = b.RefID
        AND b.Accession = c.Accession;
    """
    run_create_view(db_file, vGBRefIsolates)

    vNonClinicalIsolate = """
        CREATE VIEW vNonClinicalIsolate AS
        SELECT
            *
        FROM
            tblIsolates
        WHERE
            IsolateType IS NOT ""
    ;
    """
    run_create_view(db_file, vNonClinicalIsolate)

    vSubmissionPub = """
        CREATE VIEW vSubmissionPub AS
        SELECT
            a.*, c.*
        FROM
            tblGBRefs a,
            tblGBPubRefLink b,
            tblPublications c
        WHERE
            a.RefID = b.RefID
            AND b.PubID = c.PubID;
    """
    run_create_view(db_file, vSubmissionPub)

    vAccessionPub = """
        CREATE VIEW vAccessionPub AS
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
    run_create_view(db_file, vAccessionPub)

    vAccessionPubData = """
        CREATE VIEW vAccessionPubData AS
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
    run_create_view(db_file, vAccessionPubData)

    # vNumAccessions = """
    # CREATE VIEW vNumAccessions AS
    # SELECT
    #     COUNT(DISTINCT Accession) AS NumAccessions
    # FROM
    #     tblIsolates;
    # """
    # run_create_view(db_file, vNumAccessions)

    # vNumSubmissionSets = """
    # CREATE VIEW vNumSubmissionSets AS
    # SELECT
    #     COUNT(DISTINCT RefID) AS NumSubmissionSets
    # FROM
    #     tblGBRefs;
    # """
    # run_create_view(db_file, vNumSubmissionSets)

    vIsolateMissingData = """
    CREATE VIEW vIsolateMissingData AS
    SELECT
        *
    FROM tblIsolates
    WHERE
        IsolateType IS ""
    AND
        (
        Host IN
        ("Not Applicable", "Not Available", "Not available", "Not applicable")
    OR
        Specimen IN
        ("Not Applicable", "Not Available", "Not available", "Not applicable")
    OR
        Country IN
        ("Not Applicable", "Not Available", "Not available", "Not applicable")
    OR
        IsolateYear IN
        ("Not Applicable", "Not Available", "Not available", "Not applicable")
    );
    """
    run_create_view(db_file, vIsolateMissingData)

    vSubmissionNotMatch = """
    CREATE VIEW vSubmissionNotMatch AS
    SELECT
        *
    FROM
        tblGBRefs
    WHERE
        RefID NOT IN (
            SELECT RefID from tblGBPubRefLink
        )
    """
    run_create_view(db_file, vSubmissionNotMatch)

    vPubNotMatch = """
    CREATE VIEW vPubNotMatch AS
    SELECT
        *
    FROM
        tblGBRefs
    WHERE
        RefID NOT IN (
            SELECT RefID from tblGBPubRefLink
        )
    """
    run_create_view(db_file, vPubNotMatch)

    # vNumMatchedSubmission = """
    # CREATE VIEW vNumMatchedSubmission AS
    # SELECT
    #     (SELECT COUNT(DISTINCT RefID) FROM tblGBPubRefLink) as num_submission_set,
    #     (SELECT COUNT(DISTINCT PubID) FROM tblGBPubRefLink) as num_publication,
    #     (SELECT COUNT(DISTINCT RefID) FROM tblGBRefs WHERE
    #         RefID NOT IN (SELECT DISTINCT RefID FROM tblGBPubRefLink)
    #     ) as num_submission_not_match
    # ;
    # """
    # run_create_view(db_file, vNumMatchedSubmission)

    vNumSuppliedIsolateDataByPubMed = """
    CREATE VIEW vNumSuppliedIsolateDataByPubMed AS
    SELECT
        COUNT(
            DISTINCT
            CASE WHEN tblIsolates.Host LIKE '%*%'
                THEN tblGBRefLink.RefID
            END) AS submission_host,
        COUNT(
            DISTINCT
            CASE WHEN tblIsolates.Host LIKE '%*%'
                THEN tblIsolates.Accession
            END) AS isolate_host,
        COUNT(
            DISTINCT
            CASE WHEN tblIsolates.Specimen LIKE '%*%'
                THEN tblGBRefLink.RefID
            END) AS submission_specimen,
        COUNT(
            DISTINCT
            CASE WHEN tblIsolates.Specimen LIKE '%*%'
                THEN tblIsolates.Accession
            END) AS isolate_specimen,
        COUNT(
            DISTINCT
            CASE WHEN tblIsolates.Country LIKE '%*%'
                THEN tblGBRefLink.RefID
            END) AS submission_country,
        COUNT(
            DISTINCT
            CASE WHEN tblIsolates.Country LIKE '%*%'
                THEN tblIsolates.Accession
            END) AS isolate_country,
        COUNT(
            DISTINCT CASE WHEN tblIsolates.IsolateYear LIKE '%*%'
                THEN tblGBRefLink.RefID
            END) AS submission_IsolateYear,
        COUNT(
            DISTINCT CASE WHEN tblIsolates.IsolateYear LIKE '%*%'
                THEN tblIsolates.Accession
            END) AS isolate_IsolateYear
    FROM
        tblGBRefLink
    JOIN
        tblIsolates ON tblGBRefLink.Accession = tblIsolates.Accession
    ;
    """

    run_create_view(db_file, vNumSuppliedIsolateDataByPubMed)

    tblIsolateOrig = """
    CREATE VIEW tblIsolateOrig AS
    SELECT
        Accession,
        CASE
            WHEN IsolateYear LIKE '%*' THEN ''
            ELSE IsolateYear
        END AS IsolateYear,
        CASE
            WHEN Host LIKE '%*' THEN ''
            ELSE Host
        END AS Host,
        CASE
            WHEN Specimen LIKE '%*' THEN ''
            ELSE Specimen
        END AS Specimen,
        CASE
            WHEN Country LIKE '%*' THEN ''
            ELSE Country
        END AS Country
    FROM tblIsolates;
    """
    run_create_view(db_file, tblIsolateOrig)

    vIsolateMetadataSummary = """
    CREATE VIEW vIsolateMetadataSummary AS
    WITH temp_selection AS (
        SELECT *
        FROM
            tblIsolateOrig
            JOIN tblSequences ON tblIsolateOrig.Accession = tblSequences.Accession
            JOIN tblGBRefLink ON tblIsolateOrig.Accession = tblGBRefLink.Accession
            JOIN tblGBPubRefLink ON tblGBRefLink.RefID = tblGBPubRefLink.RefID
        WHERE
            IsolateType == ''
    )
    SELECT DISTINCT
        Host,
        Country,
        IsolateYear,
        CASE
            WHEN IsolateYear BETWEEN 1900 AND 1990 THEN '<1990'
            WHEN IsolateYear BETWEEN 1991 AND 2000 THEN '1991-2000'
            WHEN IsolateYear BETWEEN 2001 AND 2010 THEN '2001-2010'
            WHEN IsolateYear BETWEEN 2011 AND 2020 THEN '2011-2020'
            WHEN IsolateYear BETWEEN 2021 AND 2025 THEN '2021-2025'
            ELSE ''
        END AS IsolateYr,
        Gene,
        COUNT(DISTINCT Accession) AS "#" ,
        (SELECT COUNT(DISTINCT Accession) from temp_selection) AS Total,
        ROUND(
            (
                COUNT(DISTINCT Accession) * 100.0 /
                (SELECT COUNT(DISTINCT Accession) from temp_selection)
            ), 2
        ) AS "%"
    FROM
        temp_selection
    GROUP BY
        Host,
        Country,
        IsolateYr,
        Gene
    ORDER BY
        "#" DESC,
        Host,
        Country,
        IsolateYr,
        Gene;
    """
    run_create_view(db_file, vIsolateMetadataSummary)


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
