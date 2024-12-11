import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from pathlib import Path
from datetime import datetime
from Utilities import get_logger


VIRUS = 'CCHF'
SEGMENTS = ['L', 'M', 'S']
GENES = SEGMENTS
timestamp = datetime.now().strftime('%m_%d')

output_dir = Path(f"OutputData/{VIRUS}")
reference_folder = Path(f"ReferenceData/{VIRUS}")

genbank_file = reference_folder / f"{VIRUS}.gb"

BLAST_NA_DB_PATH = reference_folder / f"blast/{VIRUS}_NA_db"
BLAST_AA_DB_PATH = reference_folder / f"blast/{VIRUS}_AA_db"

genbank_feature_file = output_dir / \
    f"{VIRUS}__GenBankFeatures_{timestamp}.xlsx"
genbank_feature_check_file = output_dir / \
    f"{VIRUS}__GenBankFeatures_{timestamp}_check.xlsx"
combined_file = output_dir / f"{VIRUS}_Combined_{timestamp}.xlsx"
exclude_seq_file = output_dir / f"{VIRUS}_Excluded_Seqs_{timestamp}.xlsx"
comparison_file = output_dir / f"{VIRUS}_Combined_11_06a.xlsx"
DB_FILE = output_dir / f"{VIRUS}.db"

pubmed_folder = Path(f"Pubmed/{VIRUS}")
pubmed_file = pubmed_folder / "ReferenceSummary_Dec4.xlsx"
pubmed_additional_from_gb = pubmed_folder / "ReferenceSummary_Genbank_Dec11.xlsx"
pubmed_genbank_combined = pubmed_folder / f"{VIRUS}_P_G_Combined_{timestamp}.xlsx"

logging_file = output_dir / f'{VIRUS}_summary.txt'
logger = get_logger(logging_file)


def build_blast_db():

    aa_seqs = []
    na_seqs = []
    for s in SEGMENTS:
        with open(reference_folder / f'{s}.gb', "r") as handle:
            for record in SeqIO.parse(handle, "genbank"):
                aa_seq = [
                    i
                    for i in record.features
                    if i.type == 'CDS'
                ][0].qualifiers['translation'][0]
                aa_seqs.append(SeqRecord(Seq(aa_seq), id=s, description=''))
                na_seqs.append(
                    SeqRecord(Seq(record.seq), id=s, description=''))

    ref_aa_file = reference_folder / f"{VIRUS}_RefAAs.fasta"
    with open(ref_aa_file, "w") as output_handle:
        SeqIO.write(aa_seqs, output_handle, "fasta")

    os.system(
        f"makeblastdb -in {ref_aa_file} -dbtype prot -out {BLAST_AA_DB_PATH}")

    ref_na_file = reference_folder / f"{VIRUS}_RefNAs.fasta"
    with open(ref_na_file, "w") as output_handle:
        SeqIO.write(na_seqs, output_handle, "fasta")

    os.system(
        f"makeblastdb -in {ref_na_file} -dbtype nucl -out {BLAST_NA_DB_PATH}")


# Provides directions for cleaning the information in the feature table
def process_feature(features_df):
    features_df = translate_bio_term(features_df)
    features_df = get_additional_host_data(features_df)
    features_df['Host'] = features_df['Host2']
    features_df['isolate_source'] = features_df['isolate_source2']

    features_df['Country'] = features_df['country_region'].str.split(
        ":").str[0]

    features_df['Genes'] = features_df['segment_source']
    for i, row in features_df.iterrows():
        if str(row['Genes']) not in SEGMENTS:
            features_df.at[i, 'Genes'] = row['hit_name']

    return features_df


def translate_bio_term(features_df):
    name_map = {
        r'\bRhipicephalus\s\w+\b': 'tick',
        r'\bHyalomma\s\w+\b': 'tick',
        r'\bDermacentor\s\w+\b': 'tick',
        r'Haemaphysalis\s\w+\b': 'tick',
        r'Alveonasus\s\w+\b': 'tick',
        r'Argas\s\w+\b': 'tick',
        r'Boophilus\s\w+\b': 'tick',
        r'Ixodes\s\w+\b': 'tick',
        r'Amblyomma\s\w+\b': 'tick',
        'nymph': 'tick',

        'Mice': 'Mouse',
        'Rattus rattus': 'Rat',
        'Mus musculus': 'Mouse',
        'Capricornis milneedwardsii': 'Serow',
        'Bos taurus': 'Cattle',
        'Camelus dromedarius': 'Camel',
        'Capra': 'Goat',
        'Euchoreutes naso': 'Jerboa',
        'Testudo graeca': 'Tortoise',
        'infected mouse brain': '',
        'Suckling mouse brain': '',
    }

    features_df['Host2'] = features_df['Host']
    features_df['isolate_source2'] = features_df['isolate_source']
    for k, v in name_map.items():
        features_df['Host2'] = features_df['Host2'].str.replace(
            k, v, regex=True)
        features_df['isolate_source2'] = features_df['isolate_source2'].str.replace(
            k, v, regex=True)

    features_df['organism'] = features_df['organism'].str.replace(
        'Orthonairovirus haemorrhagiae', 'CCHF', case=False)
    features_df['organism'] = features_df['organism'].str.replace(
        r'.*Crimean-Congo hemorrhagic fever.*', 'CCHF', case=False, regex=True)

    return features_df


def get_additional_host_data(features_df):
    blood_specimen = ['blood', 'serum', 'plasma', 'sera']
    other_speciman = ['nasopharyngeal swab']
    human_host = ['patient', 'human', 'homo sapiens']
    animal_host = [
        'mouse', 'rat', 'jerboa',
        'sheep', 'camel', 'cattle', 'goat', 'serow',
        'buffalo', 'calf',
        'animal',
        'tortoise']

    for index, row in features_df.iterrows():

        host = row['Host2'].lower()
        specimen = row['isolate_source2'].lower()

        if not host and not specimen:
            continue

        updated_host = []
        updated_specimen = []

        if any(key in specimen for key in human_host):
            updated_host.append("Homo sapiens")
        if any(key in host for key in human_host):
            updated_host.append("Homo sapiens")

        for a in animal_host:
            if a in specimen:
                updated_host.append(a.capitalize())
            if a in host:
                updated_host.append(a.capitalize())

        if 'tick' in specimen:
            updated_host.append("Ticks")
        if 'tick' in host:
            updated_host.append("Ticks")

        if any(key in specimen for key in blood_specimen):
            updated_specimen.append('blood')
        if any(key in host for key in blood_specimen):
            updated_specimen.append('blood')

        for a in other_speciman:
            if a in specimen:
                updated_specimen.append(a)
            if a in host:
                updated_specimen.append(a)

        if not updated_host and host:
            updated_host = ['Other']

        if not updated_specimen and specimen:
            # specieman other and NA are the same
            updated_specimen = ['NA']

        # features_df.at[index, 'Host'] = ",".join(sorted(list(set(updated_host))))
        # features_df.at[index, 'isolate_source'] = ",".join(sorted(list(set(updated_specimen))))
        features_df.at[index, 'Host2'] = ' and '.join(
            sorted(list(set(updated_host)))) if updated_host else ''
        features_df.at[index, 'isolate_source2'] = ' and '.join(
            sorted(list(set(updated_specimen)))) if updated_specimen else ''

    return features_df


def translate_gene(gene):
    return gene if gene in SEGMENTS else ('NA' if not gene or gene == 'NA' else 'Other')
