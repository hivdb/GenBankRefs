from Bio import Entrez
import pandas as pd

from viruses.load_virus import load_virus_obj
from viruses.load_virus import select_virus

from genbank_records import select_run_blast
from genbank_records import parse_genbank_records
from genbank_records import process_references
from genbank_records import process_features
from genbank_records import process_gene_list

from DataFrameLogic import aggregate_references
from DataFrameLogic import combine_refs_and_features

from summarize_genbank import summarize_genbank
from summarize_pubmed import summarize_pubmed
from match_pubmed_GB import match_pubmed_GB

from database import create_database

# Silences the warnings that occur when empty cells are replaced with 'NA'
pd.set_option('future.no_silent_downcasting', True)
Entrez.email = "rshafer.stanford.edu"


def main():
    virus = select_virus()
    # The virus_obj contains links to pubmed tables, genbank tables
    virus_obj = load_virus_obj(virus)
    run_blast = select_run_blast()

    references, features, genes, nonvirus, nonclinical = parse_genbank_records(
        virus_obj.genbank_file)

    excludes = pd.DataFrame(nonvirus)
    print("Number of excluded records:", len(excludes))
    print('Number of non clinical records:', len(nonclinical))
    # excludes.to_excel(str(virus_obj.exclude_seq_file), index=False)

    genes = process_gene_list(genes, run_blast, virus_obj)
    print("Number of Genes:", len(genes))

    print('Number of GenBank records:', len(features))
    features = process_features(features, genes, virus_obj)

    print("Number of GenBank References:", len(references))
    references = process_references(references)

    references = aggregate_references(references, virus_obj)

    # Combine references and features
    combined_df = combine_refs_and_features(references, features, genes)
    combined_df.to_excel(str(virus_obj.combined_file), index=False)

    # Compare output file with saved file
    # saved_combined_df = pd.read_excel(str(virus_obj.comparison_file), na_values=[''])
    # compare_output_files(saved_combined_df, combined_df)

    summarize_genbank(references, features, genes, virus_obj)

    pubmed = summarize_pubmed(virus_obj.pubmed_file, virus_obj)

    if pubmed.empty:
        return

    # The virus_obj contains links to pubmed tables, genbank tables
    # the return values are: pubmed (the pubmed data file), pubmed_genbank (Pubmed and GenBank matches)
    literature, lit_ref_match = match_pubmed_GB(pubmed, references, features, genes, virus_obj)

    if literature.empty or not lit_ref_match:
        return

    features = update_genbank_by_pubmed(features, lit_ref_match)
    features.to_excel(virus_obj.genbank_feature_filled_file)

    genes = update_genes_by_features(genes, features)
    genes.to_excel(virus_obj.genbank_gene_filled_file)

    # Create database using tables:
    #   GenBank Submission Set
    #   GenBank Features
    #   Pubmed literatures
    #   Pubmed GenBank Matches
    create_database(
        virus_obj, references, features, genes,
        literature, lit_ref_match)

    virus_obj.pick_phylo_sequence(genes)


def update_genbank_by_pubmed(features, matched):

    for pubmed, genbank_list, method in matched:
        acc_list = [
            i.strip()
            for g in genbank_list
            for i in g['accession'].split(',')
        ]

        for key in ['Country', 'Host', 'IsolateType', 'SampleYr']:
            if not pubmed[key].strip() or pubmed[key].upper() == 'NA':
                pubmed[key] = ''

        for i, row in features.iterrows():
            if row['Accession'] not in acc_list:
                continue

            if not row['Country'] and pubmed['Country']:
                features.at[i, 'Country'] = pubmed['Country'] + ' *'

            if not row['Host'] and pubmed['Host']:
                features.at[i, 'Host'] = pubmed['Host'] + ' *'

            if not row['isolate_source'] and pubmed['IsolateType']:
                features.at[i, 'isolate_source'] = pubmed['IsolateType'] + ' *'

            if not row['IsolateYear'] and pubmed['SampleYr']:
                features.at[i, 'IsolateYear'] = pubmed['SampleYr'] + ' *'

    return features


def update_genes_by_features(genes, features):
    for i, g in genes.iterrows():
        feature = features[features['Accession'] == g['Accession']]
        genes.at[i, 'Host'] = feature['Host'].tolist()[0]
        genes.at[i, 'IsolateYear'] = feature['IsolateYear'].tolist()[0]
        genes.at[i, 'RecordYear'] = feature['RecordYear'].tolist()[0]
        genes.at[i, 'NonClinical'] = feature['NonClinical'].tolist()[0]
        genes.at[i, 'isolate_source'] = feature['isolate_source'].tolist()[0]
        genes.at[i, 'Country'] = feature['Country'].tolist()[0]

    return genes


if __name__ == '__main__':
    main()
