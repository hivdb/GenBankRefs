import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from pathlib import Path
from datetime import datetime


VIRUS = 'CCHF'
SEGMENTS = ['L', 'M', 'S']
timestamp = datetime.now().strftime('%m_%d')

output_dir = Path(f"OutputData/{VIRUS}")
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

    reference_folder = Path(f"ReferenceData/{VIRUS}")

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


def process_feature(features_df):
    features_df = translate_bio_term(features_df)
    features_df = get_additional_host_data(features_df)
    features_df['host'] = features_df['host2']
    features_df['isolate_source'] = features_df['isolate_source2']

    features_df['country_region'] = features_df['country_region'].str.split(":").str[0]
    for i, row in features_df.iterrows():
        if str(row['segment_source']) not in SEGMENTS:
            features_df.at[i, 'segment_source'] = row['hit_name']

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

    features_df['host2'] = features_df['host']
    features_df['isolate_source2'] = features_df['isolate_source']
    for k, v in name_map.items():
        features_df['host2'] = features_df['host2'].str.replace(k, v, regex=True)
        features_df['isolate_source2'] = features_df['isolate_source2'].str.replace(k, v, regex=True)

    features_df['organism'] = features_df['organism'].str.replace('Orthonairovirus haemorrhagiae', 'CCHF', case=False)
    features_df['organism'] = features_df['organism'].str.replace(r'.*Crimean-Congo hemorrhagic fever.*', 'CCHF', case=False, regex=True)

    return features_df


def get_additional_host_data(features_df):
    # This function is adapted for CCHF

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

        host = row['host2'].lower()
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

        # features_df.at[index, 'host'] = ",".join(sorted(list(set(updated_host))))
        # features_df.at[index, 'isolate_source'] = ",".join(sorted(list(set(updated_specimen))))
        features_df.at[index, 'host2'] = ' and '.join(
            sorted(list(set(updated_host)))) if updated_host else ''
        features_df.at[index, 'isolate_source2'] = ' and '.join(
            sorted(list(set(updated_specimen)))) if updated_specimen else ''

    return features_df
