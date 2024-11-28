from Bio import SeqIO
from Bio import Entrez
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from datetime import datetime
import pandas as pd
import os
from pathlib import Path

from GenBankFunctions import (filter_by_taxonomy,
                              pooled_blast)

from DataFrameLogic import (process_authors_titles,
                            combine_refs_and_features,
                            translate_bio_term,
                            get_additional_host_data,
                            compare_output_files)

from Utilities import (extract_year_from_journal,
                       process_author_field)

# CONSTANTS
Entrez.email = "rshafer.stanford.edu"
timestamp = datetime.now().strftime('%m_%d')
pd.set_option('display.max_rows', 100)

VIRUS = "CCHF"
RUN_BLAST = 0
genbank_file = f"ReferenceData/{VIRUS}/{VIRUS}.gb"
reference_folder = Path(f"ReferenceData/{VIRUS}")
reference_aa_file = f"ReferenceData/{VIRUS}/{VIRUS}_RefAAs.fasta"
# comparison_file = f"OutputData/{VIRUS}/{VIRUS}_Combined_11_06a.xlsx"
output_dir = f"OutputData/{VIRUS}"
SEGMENTS = ['L', 'M', 'S']


def build_blast_db(segments):

    db_name = f"{VIRUS}_AA_db"
    os.system(
        f"makeblastdb -in {reference_aa_file} -dbtype prot -out {db_name}")

    aa_seqs = []
    na_seqs = []
    for s in segments:
        with open(reference_folder / f'{s}.gb', "r") as handle:
            for record in SeqIO.parse(handle, "genbank"):
                aa_seq = [
                    i
                    for i in record.features
                    if i.type == 'CDS'
                ][0].qualifiers['translation'][0]
                aa_seqs.append(SeqRecord(Seq(aa_seq), id=s, description=''))
                na_seqs.append(SeqRecord(Seq(record.seq), id=s, description=''))

    ref_aa_file = reference_folder / f"{VIRUS}_RefAAs.fasta"
    with open(ref_aa_file, "w") as output_handle:
        SeqIO.write(aa_seqs, output_handle, "fasta")

    ref_na_file = reference_folder / f"{VIRUS}_RefNAs.fasta"
    with open(ref_na_file, "w") as output_handle:
        SeqIO.write(na_seqs, output_handle, "fasta")

    db_name = f"{VIRUS}_NA_db"
    os.system(
        f"makeblastdb -in {ref_na_file} -dbtype nucl -out {db_name}")


def main():
    build_blast_db(SEGMENTS)

    feature_list, reference_list, exclude_list = parse_genbank_records(genbank_file)

    df_excluded_seqs = pd.DataFrame(exclude_list)
    print("NumExcludedSequences:", len(df_excluded_seqs))

    print('Process genbank records', len(feature_list))

    output_file = os.path.join(
        output_dir, f"{VIRUS}__GenBankFeatures_{timestamp}.xlsx")

    if RUN_BLAST == 1:
        feature_list = pooled_blast(feature_list, VIRUS)
        # Place sequence features in a data frame
        features_df = pd.DataFrame(feature_list)
        features_df.to_excel(output_file, index=False)
    else:
        features_df = pd.read_excel(output_file).fillna('')

    features_df = translate_bio_term(features_df)
    features_df = get_additional_host_data(features_df)
    features_df['host'] = features_df['host2']
    features_df['isolate_source'] = features_df['isolate_source2']
    features_df['country_region'] = features_df['country_region'].str.split(":").str[0]
    for i, row in features_df.iterrows():
        if str(row['segment_source']) not in SEGMENTS:
            features_df.at[i, 'segment_source'] = row['hit_name']

    output_check_file = os.path.join(
        output_dir, f"{VIRUS}__GenBankFeatures_{timestamp}_check.xlsx")
    features_df.to_excel(output_check_file, index=False)

    # Aggregate by reference
    reference_df = pd.DataFrame(reference_list)
    reference_df['year'] = reference_df['journal'].apply(
        extract_year_from_journal)
    reference_df['year'] = pd.to_numeric(reference_df['year'], errors='coerce')
    reference_df['year'] = reference_df['year'].apply(
        lambda x: '' if pd.isna(x) else int(x))

    # Remove Submitted (date) from journal
    reference_df['journal'] = reference_df['journal'].str.replace(
        r"Submitted \(\d{2}-[A-Z]{3}-\d{4}\)", "", regex=True)

    reference_df['journal'] = reference_df['journal'].str.replace(
        r"(Patent).*", r"\1", regex=True)
    reference_df['authors'] = reference_df['authors'].apply(
        process_author_field)
    print("Number of original entries: ", len(reference_df))

    grouped_ref_df = reference_df.groupby(
        ['authors', 'title', 'journal', 'pmid', 'year'])[
            'accession'].apply(list).reset_index()
    grouped_ref_df['accession'] = grouped_ref_df['accession'].apply(
        lambda x: ', '.join(x))
    print("Number of entries following aggregation by metadata: ", len(grouped_ref_df))

    merged_ref_df = process_authors_titles(grouped_ref_df)
    print("Number of entries following aggregation by metadata: ", len(merged_ref_df))

    # Combine references and features
    combined_df = combine_refs_and_features(merged_ref_df, features_df)

    # acc_set = set()
    # ref_acc_number = 0
    # for i, r in combined_df.iterrows():
    #     acc_set.update(set([j.strip() for j in r['accession'].split(',')]))
    #     ref_acc_number += len([j.strip() for j in r['accession'].split(',')])
    # print(len(acc_set), 'accession number')
    # print(ref_acc_number, 'Ref duplicated accession number')

    # Print output files

    output_file = os.path.join(
        output_dir, f"{VIRUS}_Combined_{timestamp}.xlsx")
    combined_df.to_excel(output_file, index=False)

    output_file = os.path.join(
        output_dir, f"{VIRUS}_Excluded_Seqs_{timestamp}.xlsx")
    df_excluded_seqs.to_excel(output_file, index=False)

    # Compare output file with saved file
    # saved_combined_df = pd.read_excel(comparison_file, na_values=[''])
    # compare_output_files(saved_combined_df, combined_df)


def parse_genbank_records(genbank_file):
    reference_list = []
    feature_list = []
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
            feature_data = extract_features(record.features, record.id)
            features['acc_num'] = record.id

            features['description'] = record.description
            features['record_date'] = record.annotations['date']
            features['organism'] = record.annotations['organism']
            features['segment_source'] = feature_data.get('segment_source', '')
            features['cds'] = feature_data.get('product_CDS', '')
            features['host'] = feature_data.get('host_source', '')
            features['isolate_source'] = feature_data.get(
                'isolation_source_source', '')
            features['isolate_name'] = feature_data.get('isolate_source', '')
            features['country_region'] = feature_data.get(
                'geo_loc_name_source', '')
            features['collection_date'] = feature_data.get(
                'collection_date_source', '')

            sample_seq = feature_data.get('translation_CDS', '')
            features['num_na'] = len(record.seq)
            features['num_aa'] = len(sample_seq)
            features['AASeq'] = sample_seq
            features['NASeq'] = record.seq

            features['hit_name'] = ''
            features['e_value'] = 999
            features['pcnt_id'] = 0
            features['align_len'] = 0
            features['blast_name'] = ''

            feature_list.append(features)

    return feature_list, reference_list, excluded_list


def extract_references(annotations, accession):
    ref_list = annotations.get("references")
    ref_data = []
    for ref in ref_list:
        ref_items = {}
        ref_items["accession"] = accession
        ref_items["authors"] = ref.authors
        ref_items["title"] = ref.title
        ref_items["journal"] = ref.journal
        ref_items["pmid"] = ref.pubmed_id
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
