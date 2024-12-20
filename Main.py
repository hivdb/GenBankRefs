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

    references, features, genes, excludes = parse_genbank_records(
        virus_obj.genbank_file)

    excludes = pd.DataFrame(excludes)
    print("Number of excluded records:", len(excludes))
    excludes.to_excel(str(virus_obj.exclude_seq_file), index=False)

    print("Number of GenBank References:", len(references))
    print('Number of GenBank records:', len(features))
    print("Number of Genes:", len(genes))

    genes = process_gene_list(genes, run_blast, virus_obj)

    features = process_features(features, genes, virus_obj)

    references = process_references(references)

    print("Number of original entries: ", len(references))

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

    # Create database using tables:
    #   GenBank Submission Set
    #   GenBank Features
    #   Pubmed literatures
    #   Pubmed GenBank Matches
    create_database(
        virus_obj, references, features, genes,
        literature, lit_ref_match)


if __name__ == '__main__':
    main()
