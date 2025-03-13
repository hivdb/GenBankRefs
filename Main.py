from Bio import Entrez
import pandas as pd
from collections import defaultdict

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
    """
    Main function to process virus-related genomic and literature data.
    1. Loads virus object.
    2. Parses GenBank records and processes gene data.
    3. Processes features and references, filtering non-gene isolates.
    4. Aggregates references and integrates GenBank and PubMed data.
    5. Creates a database with processed genomic and literature information.
    """
    virus = select_virus()
    # The virus_obj contains links to pubmed tables, genbank tables
    virus_obj = load_virus_obj(virus)
    run_blast = select_run_blast()

    # Parse GenBank records into references, features, genes
    # Filtering out those that are not in the same virus category or are not clinical isolates
    total_references, references, features, genes, nonvirus, nonclinical = parse_genbank_records(
        virus_obj.genbank_file)

    excludes = pd.DataFrame(nonvirus)
    print("Number of excluded records:", len(excludes))
    print('Number of non clinical records:', len(nonclinical))
    # excludes.to_excel(str(virus_obj.exclude_seq_file), index=False)

    # Extract genes from all blast entries and additional detected via local alignment
    genes = process_gene_list(genes, run_blast, virus_obj)
    print("Number of Genes:", len(genes))

    # Extract features from all GenBank entries, filtering out isolates without detected genes in the features_df
    print('Number of GenBank records:', len(features))
    features = process_features(features, genes, virus_obj)

    print("Number of GenBank References:", len(references))
    acc_list = features['Accession'].tolist()

    # Extract reference (Author, Title, Journal, Year, Accessions) and combine
    # those that are from the same submission (title, author, pmid match)
    print("Number of GenBank References:", len(total_references))
    total_references = process_references(total_references)
    total_references = aggregate_references(total_references, virus_obj)

    print('-' * 80)

    # Filters the references list, keeping only those references that contains
    # at least one accession number found in features' acc_list (accessions with genes)
    references = [
        r for r in references
        if any([
            (a.strip() in acc_list) for a in r['accession'].split(',')
        ])
    ]

    print("Number of GenBank References after remove non clinical isolates:", len(references))
    references = process_references(references)
    references = aggregate_references(references, virus_obj, save_data=True)

    # Combine references and features
    combined_df = combine_refs_and_features(references, features, genes)
    combined_df.to_excel(str(virus_obj.combined_file), index=False)

    # Compare output file with saved file
    # saved_combined_df = pd.read_excel(str(virus_obj.comparison_file), na_values=[''])
    # compare_output_files(saved_combined_df, combined_df)

    # Summarize GenBank and PubMed data, see outut in datalog_genbank.txt and datalog_pubmed.txt
    summarize_genbank(references, features, genes, virus_obj)

    pubmed = summarize_pubmed(virus_obj.pubmed_file, virus_obj)

    if pubmed.empty:
        return

    # The virus_obj contains links to pubmed tables, genbank tables
    # the return values are: pubmed (the pubmed data file), pubmed_genbank (Pubmed and GenBank matches)
    literature, lit_ref_match, genbank2pubmed = match_pubmed_GB(pubmed, references, features, genes, virus_obj)

    if literature.empty or not lit_ref_match:
        return

    # Updates gene DataFrame with corresponding metadata from features DataFrame based on matching accession numbers
    genes = update_genes_by_features(genes, features)
    # Pick sequences for genes in each virus and generate phylogenetic tree - requirements vary for each
    virus_obj.pick_phylo_sequence(genes)

    # Updates features & genes DataFrame based on PubMed data on same accessions
    # features = update_genbank_by_pubmed(features, genbank2pubmed)
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


def update_genbank_by_pubmed(features, genbank2pubmed):
    # Pick best pubmed for update the accesion meta data

    feature_match_pub = defaultdict(list)
    for gen, publist, ref_id, method in genbank2pubmed:
        acc_list = [
            i.strip()
            for i in gen['accession'].split(',')
        ]
        for acc in acc_list:
            for _, pub in publist.iterrows():
                feature_match_pub[acc].append((pub, method))

    method_order = ['PMID', 'Hardlink', 'ACCESSION', 'Title']
    for acc, links in feature_match_pub.items():

        for order in method_order:
            link = [
                i
                for i in links
                if i[-1] == order
            ]
            if link:
                break

        (pubmed, method) = link[0]

        process_feature = features[features['Accession'] == acc]

        for key in ['Country', 'Host', 'IsolateType', 'SampleYr']:
            if not pubmed[key].strip() or pubmed[key].upper() == 'NA':
                pubmed[key] = ''

        for i, row in process_feature.iterrows():

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
