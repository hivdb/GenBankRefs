import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from pathlib import Path
from datetime import datetime

VIRUS = None
output_dir = None
reference_folder = None
genbank_file = None

BLAST_NA_DB_PATH = None
BLAST_AA_DB_PATH = None
genbank_feature_file = None
genbank_feature_check_file = None
combined_file = None
exclude_seq_file = None
comparison_file = None
timestamp = datetime.now().strftime('%m_%d')
DB_FILE = None

pubmed_folder = None


def set_virus(virus_name):
    global VIRUS
    global output_dir
    global genbank_file
    global genbank_feature_file
    global genbank_feature_check_file
    global combined_file
    global exclude_seq_file
    global comparison_file
    global DB_FILE
    global pubmed_folder
    global reference_folder
    global BLAST_NA_DB_PATH
    global BLAST_AA_DB_PATH
    VIRUS = virus_name

    output_dir = Path(f"OutputData/{VIRUS}")
    reference_folder = Path(f"ReferenceData/{VIRUS}")

    BLAST_NA_DB_PATH = reference_folder / f"blast/{VIRUS}_NA_db"
    BLAST_AA_DB_PATH = reference_folder / f"blast/{VIRUS}_AA_db"

    DB_FILE = output_dir / f"{VIRUS}.db"
    genbank_file = reference_folder / f"{VIRUS}.gb"
    genbank_feature_file = output_dir / f"{VIRUS}__GenBankFeatures_{timestamp}.xlsx"
    genbank_feature_check_file = output_dir / f"{VIRUS}__GenBankFeatures_{timestamp}_check.xlsx"
    combined_file = output_dir / f"{VIRUS}_Combined_{timestamp}.xlsx"
    exclude_seq_file = output_dir / f"{VIRUS}_Excluded_Seqs_{timestamp}.xlsx"
    comparison_file = output_dir / f"{VIRUS}_Combined_11_06a.xlsx"

    pubmed_folder = Path(f"Pubmed/{VIRUS}")


def build_blast_db():

    reference_aa_file = reference_folder f"/{VIRUS}_RefAAs.fasta"

    os.system(
        f"makeblastdb -in {reference_aa_file} -dbtype prot -out {BLAST_AA_DB_PATH}")


def process_feature(features_df):

    features_df['Country'] = features_df['country_region'].str.split(":").str[0]

    return features_df
