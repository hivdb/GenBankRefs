import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from pathlib import Path
from datetime import datetime

VIRUS = None
output_dir = None
genbank_file = None
genbank_feature_file = None
genbank_feature_check_file = None
combined_file = None
exclude_seq_file = None
comparison_file = None
timestamp = datetime.now().strftime('%m_%d')
DB_FILE = None


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
    VIRUS = virus_name

    output_dir = Path(f"OutputData/{VIRUS}")
    DB_File = output_dir / f"{VIRUS}.db"
    genbank_file = f"ReferenceData/{VIRUS}/{VIRUS}.gb"
    genbank_feature_file = output_dir / f"{VIRUS}__GenBankFeatures_{timestamp}.xlsx"
    genbank_feature_check_file = output_dir / f"{VIRUS}__GenBankFeatures_{timestamp}_check.xlsx"
    combined_file = output_dir / f"{VIRUS}_Combined_{timestamp}.xlsx"
    exclude_seq_file = output_dir / f"{VIRUS}_Excluded_Seqs_{timestamp}.xlsx"
    comparison_file = output_dir / f"{VIRUS}_Combined_11_06a.xlsx"


def build_blast_db():

    db_name = f"{VIRUS}_AA_db"

    reference_aa_file = f"ReferenceData/{VIRUS}/{VIRUS}_RefAAs.fasta"

    os.system(
        f"makeblastdb -in {reference_aa_file} -dbtype prot -out {db_name}")

    # reference_folder = Path(f"ReferenceData/{VIRUS}")

    # aa_seqs = []
    # na_seqs = []
    # for s in genes:
    #     with open(reference_folder / f'{s}.gb', "r") as handle:
    #         for record in SeqIO.parse(handle, "genbank"):
    #             aa_seq = [
    #                 i
    #                 for i in record.features
    #                 if i.type == 'CDS'
    #             ][0].qualifiers['translation'][0]
    #             aa_seqs.append(SeqRecord(Seq(aa_seq), id=s, description=''))
    #             na_seqs.append(SeqRecord(Seq(record.seq), id=s, description=''))

    # ref_aa_file = reference_folder / f"{VIRUS}_RefAAs.fasta"
    # with open(ref_aa_file, "w") as output_handle:
    #     SeqIO.write(aa_seqs, output_handle, "fasta")

    # ref_na_file = reference_folder / f"{VIRUS}_RefNAs.fasta"
    # with open(ref_na_file, "w") as output_handle:
    #     SeqIO.write(na_seqs, output_handle, "fasta")

    # db_name = f"{VIRUS}_NA_db"
    # os.system(
    #     f"makeblastdb -in {ref_na_file} -dbtype nucl -out {db_name}")


def process_feature(features_df):

    features_df['country_region'] = features_df['country_region'].str.split(":").str[0]

    return features_df
