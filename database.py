import pandas as pd
import sqlite3
from functools import reduce


def create_database(
        virus_obj,
        references, features, genes,
        pubmed, pubmed_genbank):

    virus_obj.DB_FILE.unlink(missing_ok=True)

    # GenBank Tables
    tblGBReferences = references[[
        'RefID', 'Authors', 'Title', 'Journal', 'PMID', 'Year']]
    dump_table(virus_obj.DB_FILE,
               'tblGBReference', tblGBReferences)

    create_ref_link(virus_obj, references)

    features['Specimen'] = features['isolate_source']
    features['Virus'] = features['organism']

    tblIsolates = features[[
        'Accession', 'Country', 'Description', 'RecordYear',
        'IsolateYear', 'Host', 'Specimen', 'IsolateName', 'Virus']]
    dump_table(
        virus_obj.DB_FILE,
        'tblIsolates', tblIsolates)

    genes['PcntMatch'] = genes['pcnt_id']
    genes['HSPLength'] = genes['align_len']

    tblSequences = genes[[
        'Accession', 'Gene', 'CDS_NAME',
        'AASeq', 'NumAA', 'AA_start', 'AA_stop',
        'NASeq', 'NumNA', 'NA_start', 'NA_stop',
        'PcntMatch', 'HSPLength',
    ]]
    dump_table(
        virus_obj.DB_FILE,
        'tblSequences',
        tblSequences)

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
        'tblLitRefLink',
        pd.DataFrame(tblLitRefLink))

    creat_views(virus_obj.DB_FILE)


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
        'tblRefLink',
        pd.DataFrame(ref_link))


def creat_views(db_file):

    vReferenceRecord = """
        CREATE VIEW vReferenceRecord AS
        SELECT a.*
        FROM tblGBReference a, tblRefLink b, tblIsolates c
        WHERE a.RefID = b.RefID
        AND b.Accession = c.Accession;
    """
    run_create_view(db_file, vReferenceRecord)

    vLitAccessionLink = """
        CREATE VIEW vLitAccessionLink AS
        SELECT
            distinct a.*
        FROM
            tblIsolates a,
            tblRefLink b,
            tblLitRefLink c
        WHERE
            a.Accession = b.Accession
            AND
            b.RefID = c.RefID;
    """
    run_create_view(db_file, vLitAccessionLink)


def dump_table(db_file, table_name, table, index=False):
    conn = sqlite3.connect(str(db_file))

    if index:
        table.to_sql(
            table_name,
            conn,
            if_exists='replace',
            index_label='ID')
    else:
        table.to_sql(
            table_name,
            conn,
            if_exists='replace',
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
