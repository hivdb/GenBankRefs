import pandas as pd
import sqlite3
from functools import reduce
import json


def create_tables(db_file):

    conn = sqlite3.connect(db_file)

    conn.execute("""
        CREATE TABLE "tblGBSubmissionSets" (
            "RefID" INTEGER PRIMARY KEY,
            "Authors" TEXT,
            "FirstAuthorSurname" TEXT,
            "Title" TEXT,
            "Journal" TEXT,
            "PMID" TEXT,
            "Year" TEXT,
            "ShortName" TEXT
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
            "NonClinical" TEXT
        )
    """)

    conn.execute("""
        CREATE TABLE "tblGBRefLink" (
            "RefID" INTEGER,
            "Accession" TEXT,
            FOREIGN KEY (RefID) REFERENCES tblGBSubmissionSets (RefID),
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
        CREATE TABLE "tblIndels" (
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
            "FirstAuthorSurname" TEXT,
            "Title" TEXT,
            "Journal" TEXT,
            "PMID" TEXT,
            "Year" TEXT,
            "ShortName" TEXT,
            "ref_source" TEXT
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
            "Method" TEXT,
            FOREIGN KEY (PubID) REFERENCES tblPublications (PubID),
            FOREIGN KEY (RefID) REFERENCES tblGBSubmissionSets (RefID)
        )
    """)

    conn.close()


def create_database(virus_obj, references, features, genes, pubmed,
                    pubmed_genbank):

    virus_obj.DB_FILE.unlink(missing_ok=True)

    create_tables(virus_obj.DB_FILE)

    # GenBank Tables
    tblGBSubmissionSets = references[[
        'RefID',
        'Authors',
        'FirstAuthorSurname',
        'Title',
        'Journal',
        'PMID',
        'Year',
        'ShortName',
    ]]
    fill_in_table(virus_obj.DB_FILE, 'tblGBSubmissionSets', tblGBSubmissionSets)

    create_ref_link(virus_obj, references)

    features['Specimen'] = features['isolate_source']
    features['Virus'] = features['organism']

    if 'NonClinical' not in features.columns:
        features['NonClinical'] = ""

    tblIsolates = features[[
        'Accession', 'Country', 'RecordYear', 'IsolateYear', 'Host',
        'Specimen', 'IsolateName', 'SeqLength', 'NonClinical'
    ]]

    for i, row in tblIsolates.iterrows():
        tblIsolates.at[i, 'Host'] = row['Host'] if row['Host'] else (
            'Not applicable' if row['NonClinical'] else 'Not available')
        tblIsolates.at[i,
                       'Specimen'] = row['Specimen'] if row['Specimen'] else (
                           'Not applicable'
                           if row['NonClinical'] else 'Not available')
        tblIsolates.at[
            i, 'IsolateYear'] = row['IsolateYear'] if row['IsolateYear'] else (
                'Not applicable' if row['NonClinical'] else 'Not available')
        tblIsolates.at[i, 'Country'] = row['Country'] if row['Country'] else (
            'Not applicable' if row['NonClinical'] else 'Not available')

    fill_in_table(virus_obj.DB_FILE, 'tblIsolates', tblIsolates)

    genes['PcntMatch'] = genes['pcnt_id']

    tblGBSequences = genes[[
        'Accession',
        'Gene',
        'CDS_NAME',
        'AA_seq',
        'AA_length',
        'AA_start',
        'AA_stop',
        'NA_seq',
        'NA_length',
        'NA_start',
        'NA_stop',
        'PcntMatch',
    ]]
    fill_in_table(virus_obj.DB_FILE, 'tblSequences', tblGBSequences)

    tblIndels = genes[(genes['AA_num_ins'] != 0) | (genes['AA_num_del'] != 0) |
                      (genes['AA_blast_failed'] == 1) |
                      (genes['NA_num_ins'] != 0) | (genes['NA_num_del'] != 0) |
                      (genes['NA_blast_failed'] == 1) | (genes['num_N'] != 0) |
                      (genes['translation_issue'] != 0)]

    tblIndels = tblIndels[[
        'SeqID', 'AA_num_ins', 'AA_num_del', 'AA_blast_failed', 'NA_num_ins',
        'NA_num_del', 'NA_blast_failed', 'num_N', 'translation_issue'
    ]]

    fill_in_table(virus_obj.DB_FILE, 'tblIndels', tblIndels)

    # PubMed Tables
    tblPublications = pubmed[[
        'PubID',
        'Authors',
        'FirstAuthorSurname',
        'Title',
        'Journal',
        'PMID',
        'Year',
        'ShortName',
        'ref_source'
    ]]
    fill_in_table(virus_obj.DB_FILE, 'tblPublications', tblPublications)

    tblPublicationData = pubmed[[
        'PubID', 'Viruses', 'NumSeqs', 'Host', 'SampleYr', 'Country',
        'GenBank', 'SeqMethod', 'CloneMethod', 'IsolateType', 'Gene'
    ]]
    fill_in_table(virus_obj.DB_FILE, 'tblPublicationData', tblPublicationData)

    tblPubRefLink = []
    for pubmed, genbank_list in pubmed_genbank:
        for g in genbank_list:
            tblPubRefLink.append(
                (pubmed['PubID'], g['RefID'], g['match_method']))

    tblPubRefLink = list(set(tblPubRefLink))
    tblPubRefLink = [{'PubID': i, 'RefID': j, 'Method': k} for i, j, k in tblPubRefLink]

    fill_in_table(virus_obj.DB_FILE, 'tblGBPubRefLink',
                  pd.DataFrame(tblPubRefLink))

    creat_views(virus_obj.DB_FILE)

    dump_db_tables(virus_obj.DB_FILE, virus_obj.db_dump_folder)

    # get_table_schema_sql(virus_obj.DB_FILE)


def create_ref_link(virus_obj, ref):
    ref_link = []
    for i, row in ref.iterrows():
        accessions = row['accession']
        accessions = [i.strip() for i in accessions.split(',') if i.strip()]
        for acc in accessions:
            ref_link.append({'RefID': row['RefID'], 'Accession': acc})

    fill_in_table(virus_obj.DB_FILE, 'tblGBRefLink', pd.DataFrame(ref_link))


def creat_views(db_file):

    vGBRefIsolates = """
        CREATE VIEW vGBRefIsolates AS
        SELECT a.*, c.*
        FROM tblGBSubmissionSets a, tblGBRefLink b, tblIsolates c
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
            NonClinical IS NOT ""
    ;
    """
    run_create_view(db_file, vNonClinicalIsolate)

    vSubmissionPub = """
        CREATE VIEW vSubmissionPub AS
        SELECT
            a.*, c.*
        FROM
            tblGBSubmissionSets a,
            tblGBPubRefLink b,
            tblPublications c
        WHERE
            a.RefID = b.RefID
            AND b.PubID = c.PubID;
    """
    run_create_view(db_file, vSubmissionPub)

    vGPMatched = """
        CREATE VIEW vGPMatched AS
        WITH temp_isolate AS (
            SELECT
                iso.*,
                seq.Gene
            FROM
                tblIsolates iso,
                tblSequences seq
            WHERE
                iso.Accession = seq.Accession
            GROUP BY
                iso.Accession, seq.Gene
        ),
        temp_num_isolate AS (
            SELECT
                a.*,
                COUNT(DISTINCT c.ACCESSION) AS num_Isolate,
                'https://www.ncbi.nlm.nih.gov/nuccore?cmd=Search&doptcmdl=Summary&term=' || GROUP_CONCAT('"' || c.ACCESSION || '"%5BACCN%5D', "%20OR%20") AS gb_search
            FROM
                tblGBSubmissionSets a,
                tblGBRefLink b,
                temp_isolate c
            WHERE
                a.RefID = b.RefID
                AND b.Accession = c.Accession
            GROUP BY
                a.RefID
        ),
        temp_num_gene AS (
            SELECT
                a.*,
                c.Gene,
                COUNT(c.Gene) AS num_Gene
            FROM
                tblGBSubmissionSets a,
                tblGBRefLink b,
                temp_isolate c
            WHERE
                a.RefID = b.RefID
                AND b.Accession = c.Accession
            GROUP BY
                a.RefID, c.Gene
            ORDER BY
                c.Gene
        ),
        temp_num_gene_str AS (
            SELECT
                a.*,
                GROUP_CONCAT(a.Gene || ' (' || a.num_Gene || ')', ', ') AS Genes
            FROM
                temp_num_gene a
            GROUP BY
                a.RefID
        ),
        temp_GBRef AS (
            SELECT
                a.*,
                b.num_Isolate,
                b.gb_search,
                c.Genes
            FROM
                tblGBSubmissionSets a
                LEFT JOIN temp_num_isolate b ON a.RefID = b.RefID
                LEFT JOIN temp_num_gene_str c ON a.RefID = c.RefID
        )

        SELECT
            a.RefID as RefID,
            a.ShortName || ". " || a.Title || "; " || a.Authors || "; " || a.Journal || "; ["|| a.Genes || "](" || a.gb_search || ")"
            AS SubmissionSet,

            c.PubID as PubID,
            c.ShortName || ". " || c.Title || "; " || c.Authors || "; " || c.Journal AS Publication,

            CASE
                WHEN c.ShortName IS NOT NULL THEN c.ShortName
                ELSE a.ShortName
            END as ShortName
        FROM
            temp_GBRef a,
            tblGBPubRefLink b,
            tblPublications c
        WHERE
            a.RefID = b.RefID
            AND b.PubID = c.PubID
            -- AND LOWER(a.Title) NOT LIKE '%patent%'
            -- AND LOWER(a.Title) NOT LIKE '%direct submission%'
            -- AND LOWER(a.Title) NOT LIKE '%construct%'
            -- AND LOWER(a.Journal) NOT LIKE '%patent%'
            -- AND LOWER(a.Journal) NOT LIKE '%direct submission%'
            -- AND LOWER(a.Journal) NOT LIKE '%construct%'
        UNION
        SELECT
            a.RefID as RefID,
            a.ShortName || ". " || a.Title || "; " || a.Authors || "; " || a.Journal || "; ["|| a.Genes || "](" || a.gb_search || ")"
            AS SubmissionSet,

            '' as PubID,
            '' AS Publication,

            a.ShortName as ShortName
        FROM
            temp_GBRef a
        WHERE
            a.RefID NOT IN (SELECT RefID FROM tblGBPubRefLink)
            -- AND LOWER(a.Title) NOT LIKE '%patent%'
            -- AND LOWER(a.Title) NOT LIKE '%direct submission%'
            -- AND LOWER(a.Title) NOT LIKE '%construct%'
            -- AND LOWER(a.Journal) NOT LIKE '%patent%'
            -- AND LOWER(a.Journal) NOT LIKE '%direct submission%'
            -- AND LOWER(a.Journal) NOT LIKE '%construct%'
        UNION
        SELECT
            '' as RefID,
            '' AS SubmissionSet,

            c.PubID as PubID,
            c.ShortName || ". " || c.Title || "; " || c.Authors || "; " || c.Journal AS Publication,

            c.ShortName as ShortName
        FROM
            tblPublications c
        WHERE
            c.PubID NOT IN (SELECT PubID FROM tblGBPubRefLink)
        ;
    """

    run_create_view(db_file, vGPMatched)

    vGPMatchedMatch = """
        CREATE VIEW vGPMatchedMatch AS
        SELECT * FROM vGPMatched
        WHERE RefID != '' AND PubID != ''
        ;
    """
    run_create_view(db_file, vGPMatchedMatch)

    vGPMatchedNoRef = """
        CREATE VIEW vGPMatchedNoRef AS
        SELECT * FROM vGPMatched
        WHERE RefID == '' AND PubID != ''
        ;
    """
    run_create_view(db_file, vGPMatchedNoRef)

    vGPMatchedNoPub = """
        CREATE VIEW vGPMatchedNoPub AS
        SELECT * FROM vGPMatched
        WHERE RefID != '' AND PubID == ''
        ;
    """
    run_create_view(db_file, vGPMatchedNoPub)

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
    #     tblGBSubmissionSets;
    # """
    # run_create_view(db_file, vNumSubmissionSets)

    vIsolateMissingData = """
    CREATE VIEW vIsolateMissingData AS
    SELECT
        *
    FROM tblIsolates
    WHERE
        NonClinical IS ""
    AND
        (
        Host IN
        ("Not Applicable", "Not Available", "Not available", "Not applicable")
    OR
        Country IN
        ("Not Applicable", "Not Available", "Not available", "Not applicable")
    OR
        IsolateYear IN
        ("Not Applicable", "Not Available", "Not available", "Not applicable")
    );
    """
    # OR
    # Specimen IN
    # ("Not Applicable", "Not Available", "Not available", "Not applicable")
    run_create_view(db_file, vIsolateMissingData)

    vSubmissionNotMatch = """
    CREATE VIEW vSubmissionNotMatch AS
    SELECT
        *
    FROM
        tblGBSubmissionSets
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
        tblGBSubmissionSets
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
    #     (SELECT COUNT(DISTINCT RefID) FROM tblGBSubmissionSets WHERE
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

    vIsolateOrig = """
    CREATE VIEW vIsolateOrig AS
    SELECT
        Accession,
        CASE
            WHEN IsolateYear LIKE '%*' THEN 'Not Available'
            ELSE IsolateYear
        END AS IsolateYear,
        CASE
            WHEN Host LIKE '%*' THEN 'Not Available'
            ELSE Host
        END AS Host,
        CASE
            WHEN Specimen LIKE '%*' THEN 'Not Available'
            ELSE Specimen
        END AS Specimen,
        CASE
            WHEN Country LIKE '%*' THEN 'Not Available'
            ELSE Country
        END AS Country,
        NonClinical
    FROM tblIsolates;
    """
    run_create_view(db_file, vIsolateOrig)

    vIsolateMetadataSummary = """
    CREATE VIEW vIsolateMetadataSummary AS
    WITH temp_selection AS (
        SELECT
            match.ShortName,
            iso.Accession,
            CASE
                WHEN iso.Host != 'Not available' THEN iso.Host
                WHEN pData.Host != '' THEN pData.Host
                ELSE iso.Host
            END AS Host,
            CASE
                WHEN iso.Country != 'Not available' THEN iso.Country
                WHEN pData.Country != '' THEN pData.Country
                ELSE iso.Country
            END AS Country,
            CASE
                WHEN iso.IsolateYear != 'Not available' THEN iso.IsolateYear
                WHEN pData.SampleYr != '' THEN pData.SampleYr
                ELSE iso.IsolateYear
            END AS IsolateYear,
            seq.Gene
        FROM
            tblIsolates iso
            LEFT JOIN tblSequences seq ON iso.Accession = seq.Accession
            JOIN tblGBRefLink ON iso.Accession = tblGBRefLink.Accession
            JOIN vGPMatched match ON tblGBRefLink.RefID = match.RefID
            LEFT JOIN tblPubLicationData pData ON match.PubID = pData.PubID
        WHERE
            NonClinical == ''
        ORDER BY
            match.ShortName,
            iso.Accession,
            seq.Gene
    )
    SELECT
        temp_selection.ShortName AS Publications,

        RTRIM(RTRIM(Host, '*'), ' ') AS Host,
        RTRIM(RTRIM(Country, '*'), ' ') As Country,
        RTRIM(RTRIM(IsolateYear, '*'), ' ') AS IsolateYear,
        -- CASE
        --     WHEN IsolateYear BETWEEN 1900 AND 1990 THEN '<1990'
        --     WHEN IsolateYear BETWEEN 1991 AND 2000 THEN '1991-2000'
        --     WHEN IsolateYear BETWEEN 2001 AND 2010 THEN '2001-2010'
        --     WHEN IsolateYear BETWEEN 2011 AND 2020 THEN '2011-2020'
        --     WHEN IsolateYear BETWEEN 2021 AND 2025 THEN '2021-2025'
        --     ELSE ''
        -- END AS IsolateYr,

        Gene,
        GROUP_CONCAT(Accession) AS Accessions,
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
        ShortName,
        Host,
        Country,
        IsolateYear,
        Gene
    ORDER BY
        Host,
        Country,
        IsolateYear,
        Gene,
        "#" DESC
        ;
    """
    run_create_view(db_file, vIsolateMetadataSummary)

    vMatchNotByPMID = """
        SELECT
            gbr.RefID,
            GROUP_CONCAT(gbrfl.Accession) AS CombinedAccessions,
            gbr.Authors,
            gbr.Title,
            gbr.Journal,
            pub.PubID,
            pub.PMID AS RefPMID,
            pub.Authors AS RefAuthor,
            pub.Title AS RefTitle,
            pub.Journal AS RefJournal,
            pubd.GenBank AS RefAccession
        FROM
            tblGBPubRefLink gbprl
        JOIN
            tblGBSubmissionSets gbr ON gbprl.RefID = gbr.RefID
        JOIN
            tblPublications pub ON gbprl.pubID = pub.pubID
        JOIN
            tblPublicationData pubd ON gbprl.pubID = pubd.pubID
        JOIN
            tblGBRefLink gbrfl ON gbrfl.refID = gbr.refID
        WHERE
            gbr.PMID IS NULL OR gbr.PMID = ''
        GROUP BY
            gbr.RefID, gbr.PMID, gbr.Authors, gbr.Title, gbr.Journal, pub.Authors, pub.Title, pub.Journal;
    """
    run_create_view(db_file, vMatchNotByPMID)

    vSubmissionPubLinkedSeqData = """
    CREATE VIEW vSubmissionPubLinkedSeqData AS
        SELECT
            match.SubmissionSet,  
            match.Publication,  
            seq.Accession AS IsolateAccession,  
        CASE
            WHEN iso.Host != 'Not available' THEN iso.Host
            WHEN pData.Host IS NOT NULL AND pData.Host != '' THEN pData.Host
            ELSE iso.Host
        END AS Host,
        CASE
            WHEN iso.Country != 'Not available' THEN iso.Country
            WHEN pData.Country IS NOT NULL AND pData.Country != '' THEN pData.Country
            ELSE iso.Country
        END AS Country,
        CASE
            WHEN iso.IsolateYear != 'Not available' THEN iso.IsolateYear
            WHEN pData.SampleYr IS NOT NULL AND pData.SampleYr != '' THEN pData.SampleYr
            ELSE iso.IsolateYear
        END AS IsolateYear,
        seq.Gene
    FROM
        tblIsolates iso
        LEFT JOIN tblSequences seq ON iso.Accession = seq.Accession
        JOIN tblGBRefLink ON iso.Accession = tblGBRefLink.Accession
        JOIN vGPMatched match ON tblGBRefLink.RefID = match.RefID
        JOIN tblPubLicationData pData  
    WHERE
        NonClinical = ''
    
    """

    run_create_view(db_file, vSubmissionPubLinkedSeqData)

    vSubmissionPubLinkedSeqDataAgg = """
    CREATE VIEW IF NOT EXISTS vSubmissionPubLinkedSeqDataAgg AS
    SELECT
        COUNT(DISTINCT SubmissionSet) AS NumSubmissions,
        SubmissionSet,
        GROUP_CONCAT(DISTINCT Publication) AS Publications,
        REPLACE(COALESCE(Host, ''), '*', '') AS Host,
        REPLACE(COALESCE(Country, ''), '*', '') AS Country,
        REPLACE(COALESCE(IsolateYear, ''), '*', '') AS IsolateYear,
        GROUP_CONCAT(DISTINCT Gene) AS Genes,
        GROUP_CONCAT(DISTINCT IsolateAccession) AS Accessions,  
        COUNT(DISTINCT IsolateAccession) AS NumAccessions
    FROM
        vSubmissionPubLinkedSeqData
    GROUP BY
        SubmissionSet;
    """
    run_create_view(db_file, vSubmissionPubLinkedSeqDataAgg)



def fill_in_table(db_file, table_name, table):
    conn = sqlite3.connect(str(db_file))

    table.to_sql(table_name, conn, if_exists='append', index=False)


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

    cursor = conn.execute(
        f"SELECT name FROM sqlite_master WHERE type='table' AND name NOT LIKE 'sqlite_%';"
    )

    for i in cursor.fetchall():
        table_name = i[0]

        cursor = conn.execute(
            f"SELECT sql FROM sqlite_master WHERE type='table' AND name='{table_name}';"
        )

        schema = cursor.fetchone()[0]

        print(schema)

    conn.close()


def dump_db_tables(db_path, db_dump_folder):
    tables = [
        # 'tblGBSubmissionSets',
        'vIsolateMissingData',
        # 'vSubMissionNotMatch',
        'vIsolateMetadataSummary',
        # 'vNonClinicalIsolate',

        # 'tblPublications',
        # 'tblSubmissionPub'

        # 'vNumSuppliedIsolateDataByPubMed',
        'vGPMatchedMatch',
        'vGPMatchedNoRef',
        'vGPMatchedNoPub',
    ]

    for t in tables:
        json_file_path = db_dump_folder / f"{t}.json"
        dump_table_to_json(json_file_path, db_path, t)


def dump_table_to_json(json_file_path, db_path, table_name):
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()

    cursor.execute(f"SELECT * FROM {table_name}")
    rows = cursor.fetchall()

    column_names = [description[0] for description in cursor.description]

    data = [dict(zip(column_names, row)) for row in rows]

    with open(json_file_path, 'w', encoding='utf-8') as json_file:
        json.dump(data, json_file, indent=4, ensure_ascii=False)
    # print(
    #     f"Data from view '{table_name}' has been exported to {json_file_path}")

    conn.commit()
    conn.close()


def dump_view_to_excel(virus_obj, db_path):
    conn = sqlite3.connect(db_path)

    query = "SELECT * FROM vMatchNotByPMID"

    # Execute the query and fetch data into a pandas DataFrame
    df = pd.read_sql_query(query, conn)

    # Ensure output directory exists
    output_dir = f"OutputData/{virus_obj}"

    # Save the DataFrame to an Excel file
    excel_path = f"{output_dir}/{virus_obj}_matched_except_PMID.xlsx"
    df.to_excel(excel_path, index=False)

    # Close the database connection
    conn.close()

    return excel_path
