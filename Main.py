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
from DataFrameLogic import remove_no_pmid_ref_by_linked_accession
from DataFrameLogic import combine_refs_and_features

from summarize_genbank import summarize_genbank
from summarize_pubmed import summarize_pubmed
from match_pubmed_GB import match_pubmed_GB
from pubmed_search import search_by_pubmed_API
from AI_match_paper import using_ai_match

from database import create_database

from Utilities import yes_no

# Silences the warnings that occur when empty cells are replaced with 'NA'
pd.set_option('future.no_silent_downcasting', True)
Entrez.email = "rshafer.stanford.edu"


def main():
    """
    Main function to process virus-related genomic and literature data.
    - Loads virus object.
    - Parses GenBank records and processes gene data.
    - Processes features and references, filtering non-gene isolates.
    - Search publications using PubMed API, or GPT-4o API
    - Creates a database with processed genomic and literature information.
    - Aggregates references and integrates GenBank and PubMed data.
    """
    virus = select_virus()
    # The virus_obj contains links to pubmed tables, genbank tables
    virus_obj = load_virus_obj(virus)

    references, isolates, features, genes = extract_genbank_ref_feature_gene(virus_obj)

    for idx, row in references.iterrows():
        if row['PMID']:
            references.loc[idx, 'PMID_Source'] = 'GenBank'
        else:
            references.loc[idx, 'PMID_Source'] = ''
    references = find_publications_by_PubMed(virus_obj, references)

    references = find_publications_by_AI(virus_obj, references)

    literature, lit_ref_match, genbank2pubmed = find_publication_by_sys_review(
        virus_obj, references, features, genes)

    update_genbank_by_publication(virus_obj, features, genes, genbank2pubmed)

    # Create database using tables:
    #   GenBank Submission Set
    #   GenBank Features
    #   Pubmed literatures
    #   Pubmed GenBank Matches
    create_database(
        virus_obj, references, isolates, features, genes,
        literature, lit_ref_match)


def extract_genbank_ref_feature_gene(virus_obj):
    run_blast = select_run_blast()

    # Parse GenBank records into references, features, genes
    # Filtering out those that are not in the same virus category or are not clinical isolates
    (
        total_references, references,
        features, genes, exclude_acc_list
    ) = parse_genbank_records(
        virus_obj.genbank_file)

    # Extract genes from all blast entries and additional detected via local alignment
    genes = process_gene_list(genes, run_blast, virus_obj)

    # Extract features from all GenBank entries, filtering out isolates without detected genes in the features_df
    print('# GenBank Accessions:', len(features))
    features, exclude_features = process_features(features, genes, virus_obj)

    exclude_acc_list = pd.concat([exclude_acc_list, exclude_features])
    print('# Total exclude accessions:', len(exclude_acc_list))

    acc_list = features['Accession'].tolist()
    genes = genes[genes['Accession'].isin(acc_list)]
    print("# Genes:", len(genes))

    print('-' * 80)

    print("# Total GenBank References:", len(total_references))

    before_reference = references.copy()
    print("# Before excluding GenBank References:", len(before_reference))

    before_reference = process_references(before_reference)
    before_reference = aggregate_references(before_reference, virus_obj, save_data=False)
    before_reference = remove_no_pmid_ref_by_linked_accession(virus_obj, before_reference)

    # Extract reference (Author, Title, Journal, Year, Accessions) and combine
    # those that are from the same submission (title, author, pmid match)
    # total_references = process_references(total_references)
    # total_references = aggregate_references(total_references, virus_obj)

    print('-' * 80)

    # Filters the references list, keeping only those references that contains
    # at least one accession number found in features' acc_list (accessions with genes)
    references = [
        r for r in references
        if any([
            (a.strip() in acc_list) for a in r['accession'].split(',')
        ])
    ]

    print("# GenBank References after remove excluded accessions:", len(references))

    references = process_references(references)
    references = aggregate_references(references, virus_obj, save_data=True)
    references = remove_no_pmid_ref_by_linked_accession(virus_obj, references)

    # Combine references and features
    combined_df = combine_refs_and_features(references, features, genes)
    combined_df.to_excel(str(virus_obj.combined_file), index=False)

    # Compare output file with saved file
    # saved_combined_df = pd.read_excel(str(virus_obj.comparison_file), na_values=[''])
    # compare_output_files(saved_combined_df, combined_df)

    # Summarize GenBank and PubMed data, see outut in datalog_genbank.txt and datalog_pubmed.txt

    if yes_no('Summarize GenBank data?', True):
        summarize_genbank(references, features, genes, virus_obj)

    isolates = pd.read_excel(virus_obj.isolate_file)

    isolates = pd.DataFrame([
        {
            'IsolateID': idx + 1,
            'Accession': acc.strip(),
        }
        for idx, row in isolates.iterrows()
        for acc in row['Accession'].split(',')
    ])

    # Updates gene DataFrame with corresponding metadata from features DataFrame based on matching accession numbers
    genes = update_genes_by_features(genes, features)
    # Pick sequences for genes in each virus and generate phylogenetic tree - requirements vary for each

    if yes_no('Generate phylogenetics?', True):
        virus_obj.pick_phylo_sequence(genes)

    return references, isolates, features, genes


def find_publications_by_PubMed(virus, genbank):

    pubmed_result = search_by_pubmed_API(virus, genbank)

    counter = 0
    for idx, g in genbank.iterrows():
        if g['PMID']:
            continue
        search_r = pubmed_result[pubmed_result['RefID'] == g['RefID']]
        if search_r.empty:
            print('RefID', g['RefID'], 'is not found in file', virus.pubmed_search_result)
        pmid = search_r.iloc[0]['PMID']
        if not pd.isna(pmid):
            try:
                genbank.loc[idx, 'PMID'] = str(int(pmid))
            except ValueError:
                genbank.loc[idx, 'PMID'] = str(pmid)
            counter += 1
            genbank.loc[idx, 'PMID_Source'] = 'PubMed'

    print(counter, 'sets find publications by PubMed API.')
    return genbank


def find_publications_by_AI(virus, genbank):

    if not virus.AI_search_result:

        genbank_no_pmid_list = genbank[genbank['PMID'] == '']
        # print('# submission sets without PMID:', len(genbank_no_pmid_list))
        using_ai_match(virus, genbank_no_pmid_list, file_suffix='using_AI_find_paper')

        print('Please check AI Search result by hand.')
        exit()

    pubmed_result = pd.read_excel(virus.AI_search_result)

    counter = 0
    for idx, g in genbank.iterrows():
        if g['PMID']:
            continue
        search_r = pubmed_result[pubmed_result['RefID'] == g['RefID']]
        if search_r.empty:
            print('RefID', g['RefID'], 'is not found in file',
                  virus.AI_search_result)
        pmid = search_r.iloc[0]['PMID']
        if not pd.isna(pmid):
            try:
                genbank.loc[idx, 'PMID'] = str(int(pmid))
            except ValueError:
                genbank.loc[idx, 'PMID'] = str(pmid)
            counter += 1
            genbank.loc[idx, 'PMID_Source'] = 'AI'

    print(counter, 'sets find publications by AI.')
    return genbank


def find_publication_by_sys_review(virus_obj, references, features, genes):
    if not virus_obj.pubmed_file.exists():
        print('Please prepare a file with extract metadata from publications')
        exit()

    pubmed = pd.read_excel(virus_obj.pubmed_file, dtype=str).fillna('')
    pubmed['ref_source'] = 'PubMed search'
    if yes_no('Summarize pubmed systematic review?'):
        pubmed = summarize_pubmed(pubmed, virus_obj)

    # The virus_obj contains links to pubmed tables, genbank tables
    # the return values are: pubmed (the pubmed data file), pubmed_genbank (Pubmed and GenBank matches)
    literature, lit_ref_match, genbank2pubmed = match_pubmed_GB(
        pubmed, references, features, genes, virus_obj)

    return literature, lit_ref_match, genbank2pubmed


def update_genbank_by_publication(virus_obj, features, genes, genbank2pubmed):
    # Updates features & genes DataFrame based on PubMed data on same accessions
    features = update_genbank_by_pubmed(features, genbank2pubmed)
    features.to_excel(virus_obj.genbank_feature_filled_file)
    genes = update_genes_by_features(genes, features)
    genes.to_excel(virus_obj.genbank_gene_filled_file)


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

        for key in ['Country', 'Host', 'IsolateType', 'SampleYr', 'SeqMethod']:
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

            if not row['SeqMethod'] and pubmed['SeqMethod']:
                features.at[i, 'SeqMethod'] = pubmed['SeqMethod'] + ' *'

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
