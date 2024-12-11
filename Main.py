import importlib.util
from Bio import SeqIO
from Bio import Entrez
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import importlib
import database

from datetime import datetime
import pandas as pd
from pathlib import Path

from GenBankFunctions import (filter_by_taxonomy,
                              pooled_blast)

from DataFrameLogic import (process_authors_titles,
                            combine_refs_and_features,
                            compare_output_files)

from Utilities import (extract_year_from_journal,
                       extract_year_from_date_fields,
                       process_author_field)
from compare_pubmed_genbank import compare_pubmed_genbank

# Silences the warnings that occur when empty cells are replaced with 'NA'
pd.set_option('future.no_silent_downcasting', True)
Entrez.email = "rshafer.stanford.edu"
timestamp = datetime.now().strftime('%m_%d')

# This is called by the main function to allow users to select the virus


def select_virus():
    viruses_list = (
        ('Orthonairovirus haemorrhagiae', 'CCHF'),
        ('Henipavirus nipahense', 'Nipah'),
    )
    for index, (cname, name) in enumerate(viruses_list):
        print(f"{index + 1}.", cname, f"({name})")

    virus_id = input('Please select a virus by ID: ')
    assert virus_id.isdigit(), 'Virus not found'
    assert int(virus_id) <= len(viruses_list), 'Virus not found'

    return viruses_list[int(virus_id) - 1][-1]

# This is called by the main function to allow users to determine whether
# BLAST is run


def select_run_blast(default=None):
    if default is not None:
        return default

    result = input('Run blast? [y/n]: ')
    result = result.lower()
    assert (result in ['y', 'n']), "Please use y/n."

    return 1 if result == 'y' else 0


def load_virus_obj(virus):
    """
        virus_conf folder contains specific configuration (virus config) for a virus.
        If this file doesn't exist, the default.py file will be used.

        A virus config contains the input and output file path, and functions
        for cleaning the data or functions for running blast.
    """
    virus_conf_path = f'viruses_conf/{virus}.py'
    if not Path(virus_conf_path).exists():
        spec = importlib.util.spec_from_file_location(
            'viruses_conf.default', 'viruses_conf/default.py')
    else:
        spec = importlib.util.spec_from_file_location(
            f'viruses_conf.{virus}', virus_conf_path)

    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)

    if not Path(virus_conf_path).exists():
        module.set_virus(virus)

    return module


def main():
    virus = select_virus()
    virus_obj = load_virus_obj(virus)
    run_blast = select_run_blast()

    virus_obj.build_blast_db()

    feature_list, reference_list, exclude_list = parse_genbank_records(
        virus_obj.genbank_file)

    df_excluded_seqs = pd.DataFrame(exclude_list)
    print("Number of excluded sequences:", len(df_excluded_seqs))
    print('Number of genbank records', len(feature_list))

    if run_blast == 1:
        feature_list = pooled_blast(feature_list, virus_obj)
        features_df = pd.DataFrame(feature_list)
        features_df.to_excel(str(virus_obj.genbank_feature_file), index=False)
        features_df = pd.read_excel(
            str(virus_obj.genbank_feature_file)).fillna('')

    # genbank_feature_file is created either using blast or without blast
    # depending on the user input
    # This statement loads the existing genbank feature file which may or may not
    # include blast data. This makes it possible to avoid re-running BLAST
    elif virus_obj.genbank_feature_file.exists():
        features_df = pd.read_excel(
            str(virus_obj.genbank_feature_file)).fillna('')

    # This statement creates features_df without running BLAST
    else:
        features_df = pd.DataFrame(feature_list)
        features_df.to_excel(str(virus_obj.genbank_feature_file), index=False)
        features_df = pd.read_excel(
            str(virus_obj.genbank_feature_file)).fillna('')

    # This uses data in the imported virus module to clean data in the feature table
    features_df = virus_obj.process_feature(features_df)

    features_df['RecordYear'] = features_df['record_date'].apply(extract_year_from_date_fields)
    features_df['IsolateYear'] = features_df['collection_date'].apply(extract_year_from_date_fields)

    features_df.to_excel(
        str(virus_obj.genbank_feature_check_file), index=False)

    # Aggregate by reference
    reference_df = pd.DataFrame(reference_list)
    reference_df['Year'] = reference_df['Journal'].apply(
        extract_year_from_journal)
    reference_df['Year'] = pd.to_numeric(reference_df['Year'], errors='coerce')
    reference_df['Year'] = reference_df['Year'].apply(
        lambda x: '' if pd.isna(x) else int(x))

    # Remove Submitted (date) from journal
    reference_df['Journal'] = reference_df['Journal'].str.replace(
        r"Submitted \(\d{2}-[A-Z]{3}-\d{4}\)", "", regex=True)

    reference_df['Journal'] = reference_df['Journal'].str.replace(
        r"(Patent).*", r"\1", regex=True)
    reference_df['Authors'] = reference_df['Authors'].apply(
        process_author_field)
    print("Number of original entries: ", len(reference_df))

    grouped_ref_df = reference_df.groupby(
        ['Authors', 'Title', 'Journal', 'PMID', 'Year'])[
            'accession'].apply(list).reset_index()
    grouped_ref_df['accession'] = grouped_ref_df['accession'].apply(
        lambda x: ', '.join(x))
    print("Number of entries following aggregation by exact matches: ",
          len(grouped_ref_df))

    merged_ref_df = process_authors_titles(grouped_ref_df)
    merged_ref_df['SetID'] = merged_ref_df.index + 1
    print("Number of entries following aggregation by similarity: ",
          len(merged_ref_df))

    # Combine references and features
    combined_df = combine_refs_and_features(merged_ref_df, features_df)

    combined_df.to_excel(str(virus_obj.combined_file), index=False)

    df_excluded_seqs.to_excel(str(virus_obj.exclude_seq_file), index=False)

    # Compare output file with saved file
    # saved_combined_df = pd.read_excel(str(virus_obj.comparison_file), na_values=[''])
    # compare_output_files(saved_combined_df, combined_df)

    # The virus_obj contains links to pubmed tables, genbank tables
    # compare_pubmed_genbank function will load these tables to pandas dataframes internally
    # the return values are: pubmed (the pubmed data file), pubmed_genbank (Pubmed and GenBank matches)
    pubmed, pubmed_genbank = compare_pubmed_genbank(virus_obj)

    # Create database using tables:
    #   GenBank Submission Set
    #   GenBank Features
    #   Pubmed literatures
    #   Pubmed GenBank Matches
    create_database(
        virus_obj, merged_ref_df, features_df, pubmed, pubmed_genbank)


def create_database(virus_obj, merged_ref_df, features_df, pubmed, pubmed_genbank):
    virus_obj.DB_FILE.unlink(missing_ok=True)

    tblSubmissionSet = merged_ref_df[[
        'SetID', 'Authors', 'Title', 'Journal', 'PMID', 'Year']]
    database.dump_table(virus_obj.DB_FILE, 'tblSubmissionSet', tblSubmissionSet)

    create_ref_link(virus_obj, merged_ref_df)

    features_df['Specimen'] = features_df['isolate_source']
    features_df['Virus'] = features_df['organism']

    tblIsolates = features_df[[
        'Accession', 'Country', 'Description', 'RecordYear',
        'IsolateYear', 'Host', 'Specimen', 'IsolateName', 'Virus']]
    database.dump_table(virus_obj.DB_FILE, 'tblIsolates', tblIsolates)

    features_df['PcntMatch'] = features_df['pcnt_id']
    features_df['HSPLength'] = features_df['align_len']
    tblSequences = features_df[[
        'Accession', 'Genes', 'AASeq', 'NASeq', 'NumAA', 'NumNA',
        'PcntMatch', 'HSPLength',
    ]]
    database.dump_table(virus_obj.DB_FILE, 'tblSequences', tblSequences)

    tblLitReferences = pubmed[[
        'RefID',
        'Authors',
        'Title',
        'Journal',
        'PMID',
        'Year'
    ]]
    database.dump_table(virus_obj.DB_FILE, 'tblLitReferences', tblLitReferences)

    tblLitTextReview = pubmed[[
        'RefID',
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
    database.dump_table(virus_obj.DB_FILE, 'tblLitTextReview', tblLitTextReview)

    tblLitSubmitLink = []
    for pubmed, genbank_list in pubmed_genbank:
        for g in genbank_list:
            tblLitSubmitLink.append((pubmed['RefID'], g['SetID']))

    tblLitSubmitLink = list(set(tblLitSubmitLink))
    tblLitSubmitLink = [
        {
            'RefID': i,
            'SetID': j
        }
        for i, j in tblLitSubmitLink
    ]

    database.dump_table(
        virus_obj.DB_FILE, 'tblLitSubmitLink', pd.DataFrame(tblLitSubmitLink))


def create_ref_link(virus_obj, ref):
    ref_link = []
    for i, row in ref.iterrows():
        accessions = row['accession']
        accessions = [i.strip() for i in accessions.split(',') if i.strip()]
        for acc in accessions:
            ref_link.append({
                'SetID': row['SetID'],
                'accession': acc
            })

    database.dump_table(virus_obj.DB_FILE, 'tblRefLink', pd.DataFrame(ref_link))


def parse_genbank_records(genbank_file):
    reference_list = []
    feature_list = []
    gene_list = []

    excluded_list = []
    total_record = 0
    with open(genbank_file, "r") as handle:
        for record in SeqIO.parse(handle, "genbank"):

            total_record += 1
            if total_record % 10000 == 0:
                print(total_record)

            taxonomy = record.annotations['taxonomy']
            if len(taxonomy) == 0 or taxonomy[0] != 'Viruses':
                excluded_seq_data = filter_by_taxonomy(record)
                excluded_list.append(excluded_seq_data)
                continue

            ref_data = extract_references(record.annotations, record.id)
            reference_list.extend(ref_data)

            features = {}
            features['Accession'] = record.id

            features['Description'] = record.description
            features['record_date'] = record.annotations['date']
            features['organism'] = record.annotations['organism']

            print(record.features)
            raise

            feature_data = extract_features(record.features, record.id)
            features['segment_source'] = feature_data.get('segment_source', '')
            features['cds'] = feature_data.get('product_CDS', '')
            features['Host'] = feature_data.get('host_source', '')
            features['isolate_source'] = feature_data.get(
                'isolation_source_source', '')
            features['IsolateName'] = feature_data.get('isolate_source', '')
            features['country_region'] = feature_data.get(
                'geo_loc_name_source', '')
            features['collection_date'] = feature_data.get(
                'collection_date_source', '')

            sample_seq = feature_data.get('translation_CDS', '')
            features['NumNA'] = len(record.seq)
            features['NumAA'] = len(sample_seq)
            features['AASeq'] = sample_seq
            features['NASeq'] = record.seq

            features['hit_name'] = ''
            features['e_value'] = 999
            features['pcnt_id'] = 0
            features['align_len'] = 0
            features['blast_name'] = ''
            features['AA_start'] = ''
            features['AA_stop'] = ''
            features['NA_start'] = ''
            features['NA_stop'] = ''

            feature_list.append(features)

    return feature_list, reference_list, excluded_list


def extract_references(annotations, accession):
    ref_list = annotations.get("references")
    ref_data = []
    for ref in ref_list:
        ref_items = {}
        ref_items['accession'] = accession
        ref_items['Authors'] = ref.authors
        ref_items['Title'] = ref.title
        ref_items['Journal'] = ref.journal
        ref_items['PMID'] = ref.pubmed_id
        ref_data.append(ref_items)
    return ref_data


def extract_features(features, accession):
    feature_items = {}
    feature_items["accession"] = accession
    for feature in features:
        for key, value in feature.qualifiers.items():
            key = key + '_' + feature.type
            feature_items[key] = value[0]
            # print(f"Key {key}: value {value}\n")
    return feature_items


if __name__ == '__main__':
    main()
