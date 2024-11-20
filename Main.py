from Bio import SeqIO
from Bio import Entrez
from datetime import datetime
import pandas as pd
import os

from GenBankFunctions import (filter_by_taxonomy,
                              pooled_blast)

from DataFrameLogic import (process_authors_titles,
                            combine_refs_and_features,
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
reference_aa_file = f"ReferenceData/{VIRUS}/{VIRUS}_RefAAs.fasta"
comparison_file = f"OutputData/{VIRUS}/{VIRUS}_Combined_11_20.xlsx"
output_dir = f"OutputData/{VIRUS}"


def main():
    db_name = "ref_db"
    if (RUN_BLAST == 1):
        os.system(
            f"makeblastdb -in {reference_aa_file} -dbtype prot -out {db_name}")

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

    reference_list = []
    feature_list = []
    excluded_seqs = []
    count = 0
    with open(genbank_file, "r") as handle:

        for record in SeqIO.parse(handle, "genbank"):
            count += 1
            taxonomy = record.annotations['taxonomy']
            if len(taxonomy) == 0 or taxonomy[0] != 'Viruses':
                excluded_seq_data = filter_by_taxonomy(record)
                excluded_seqs.append(excluded_seq_data)
                continue

            ref_data = extract_references(record.annotations, record.id)
            reference_list.extend(ref_data)

            features = {}
            feature_data = extract_features(record.features, record.id)
            features['acc_num'] = record.id
            sample_seq = feature_data.get('translation_CDS', '')
            features['num_na'] = len(record.seq)
            features['num_aa'] = len(sample_seq)
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
            features['e_value'] = 999
            features['pcnt_id'] = 0
            features['align_len'] = 0
            features['_sample_seq'] = sample_seq
            feature_list.append(features)

    if RUN_BLAST == 1:
        feature_list = pooled_blast(feature_list, db_name)

    df_excluded_seqs = pd.DataFrame(excluded_seqs)
    print("NumExcludedSequences:", len(df_excluded_seqs))

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

    # output_file = os.path.join(
    #     output_dir, f"{VIRUS}__GenBankRefs_before_group_{timestamp}.xlsx")
    # reference_df.to_excel(output_file, index=False)

    grouped_ref_df = reference_df.groupby(
        ['authors', 'title', 'journal', 'pmid', 'year'])[
            'accession'].apply(list).reset_index()
    grouped_ref_df['accession'] = grouped_ref_df['accession'].apply(
        lambda x: ', '.join(x))
    print("Number of entries following aggregation by metadata: ", len(grouped_ref_df))

    # output_file = os.path.join(
    #     output_dir, f"{VIRUS}__GenBankRefs_before_merge_{timestamp}.xlsx")
    # grouped_ref_df.to_excel(output_file)

    # merged_ref_df = process_accession_lists(grouped_ref_df)
    # print("Number of entries following aggregation by accession numbers: ", len(merged_ref_df))
    merged_ref_df = process_authors_titles(grouped_ref_df)
    print("Number of entries following aggregation by metadata: ", len(merged_ref_df))

    # output_file = os.path.join(
    #     output_dir, f"{VIRUS}__GenBankRefs_after_merge_{timestamp}.xlsx")
    # merged_ref_df.to_excel(output_file)

    # Place sequence features in a data frame
    features_df = pd.DataFrame(feature_list)
    features_df = get_additional_host_data(features_df)
    features_df['country_region'] = features_df['country_region'].str.split(":").str[0]
    # Combine references and features
    combined_df = combine_refs_and_features(merged_ref_df, features_df)

    # Print output files
    output_file = os.path.join(
        output_dir, f"{VIRUS}__GenBankRefs_{timestamp}.xlsx")
    merged_ref_df.to_excel(output_file, index=False)

    output_file = os.path.join(
        output_dir, f"{VIRUS}__GenBankFeatures_{timestamp}.xlsx")
    features_df.to_excel(output_file, index=False)

    output_file = os.path.join(
        output_dir, f"{VIRUS}_Combined_{timestamp}.xlsx")
    combined_df.to_excel(output_file, index=False)

    output_file = os.path.join(
        output_dir, f"{VIRUS}_Excluded_Seqs_{timestamp}.xlsx")
    df_excluded_seqs.to_excel(output_file, index=False)

    # Compare output file with saved file
    saved_combined_df = pd.read_excel(comparison_file, na_values=[''])
    compare_output_files(saved_combined_df, combined_df)


if __name__ == '__main__':
    main()
