import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from pathlib import Path
from datetime import datetime


VIRUS = 'Nipah'
GENES = ['N', 'P', 'M', 'F', 'G', 'L']
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
pubmed_file = pubmed_folder / "ReferenceSummary_Dec11.xlsx"
pubmed_additional_from_gb = None
pubmed_genbank_combined = pubmed_folder / f"{VIRUS}_P_G_Combined_{timestamp}.xlsx"


def build_blast_db():
    aa_seqs = []
    na_seqs = []
    with open(reference_folder / 'NC_002728.gb', "r") as handle:
        for record in SeqIO.parse(handle, "genbank"):
            na_seqs.append(
                SeqRecord(record.seq, id='genome', description=''))
            gene_seq = [
                i
                for i in record.features
                if i.type == 'CDS'
            ]
            for aa in gene_seq:
                gene = None
                if 'gene' in aa.qualifiers:
                    gene = aa.qualifiers['gene'][0].upper()
                else:
                    gene = aa.qualifiers['product'][0].split(' ')[0].upper()

                if gene not in GENES:
                    continue

                aa_seqs.append(
                    SeqRecord(Seq(aa.qualifiers['translation'][0]), id=gene, description=''))
                na_seqs.append(
                    SeqRecord(Seq(aa.location.extract(record.seq)), id=gene, description=''))

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


def process_feature(features_df):
    features_df = translate_bio_term(features_df)
    features_df = get_additional_host_data(features_df)
    features_df['Host'] = features_df['Host2']
    features_df['isolate_source'] = features_df['isolate_source2']

    features_df['Country'] = features_df['country_region'].str.split(
        ":").str[0]

    features_df['cds'] = features_df['cds'].apply(translate_cds_name)

    features_df['Genes'] = features_df['cds']
    for i, row in features_df.iterrows():
        if str(row['Genes']) not in GENES:
            features_df.at[i, 'Genes'] = row['hit_name']

        if int(row['NumNA']) > 17000:
            features_df.at[i, 'Genes'] = 'genome'

    return features_df


def translate_cds_name(cds):
    cds = cds.replace('protein', '').strip()
    name_map = {
        'RNA-directed RNA polymerase': 'L',
        'nucleocapsid': 'N',
        'polymerase': 'L',
        'fusion': 'F',
        'phosphoprotein': 'P',
        'large polymerase': 'L',
        'glycoprotein': 'G',
        'matrix': 'M',
        'nucleoprotein': 'N',
        'RNA polymerase': 'L'
    }

    for k, v in name_map.items():
        cds = cds.replace(k, v)

    return cds


def translate_bio_term(features_df):

    name_map = {
        r'Pteropus\s\w+\b': 'bat',
        'Homo sapiens; male': 'Homo sapiens',
        'Sus scrofa domesticus': 'Pig',
        'Canis lupus familiaris': 'Dog',
        'Sus scrofa (pig)': 'Pig',
    }

    features_df['Host2'] = features_df['Host']
    features_df['isolate_source2'] = features_df['isolate_source']
    for k, v in name_map.items():
        features_df['Host2'] = features_df['Host2'].str.replace(
            k, v, regex=True)
        features_df['isolate_source2'] = features_df['isolate_source2'].str.replace(
            k, v, regex=True)

    features_df['organism'] = features_df['organism'].str.replace(
        'Henipavirus nipahense', 'Nipah', case=False)

    return features_df


def get_additional_host_data(features_df):
    blood_specimen = ['blood', 'serum', 'plasma', 'sera']
    other_speciman = [
        'brain',
        'breast milk',
        'csf',
        'heart',
        'intestine',
        'kidney',
        'liver',
        'spleen',
        'lung',
        'oropharyngeal swab',
        'urine',
        'throat swab',
        ]
    human_host = ['patient', 'human', 'homo sapiens']
    animal_host = [
        'bat', 'pig', 'dog']

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
            sorted(list(set(updated_host)))) if updated_host else 'NA'
        features_df.at[index, 'isolate_source2'] = ' and '.join(
            sorted(list(set(updated_specimen)))) if updated_specimen else 'NA'

    return features_df
