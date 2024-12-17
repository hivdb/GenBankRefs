from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from .virus import Virus
import subprocess
from functools import partial


class Lassa(Virus):

    @property
    def SEGMENTS(self):
        return ['L segment', 'S segment']

    @property
    def GENES(self):
        return [
            'nucleoprotein',
            'glycoprotein',
            'Z',
            'L',
        ]

    @property
    def pubmed_file(self):
        return self.pubmed_folder / "Reference_Summary_Dec17.xlsx"

    @property
    def pubmed_additional_from_gb(self):
        return None

    def build_blast_db(self):
        build_blast_db(self)

    def process_feature(self, features_df):
        return process_feature(features_df)

    def process_gene_list(self, gene_df):
        return process_gene_list(self, gene_df)

    def translate_gene(self, gene):
        return translate_gene(self, gene)

    def process_pubmed(self, pubmed):
        pubmed['Gene'] = pubmed['Gene'].apply(partial(translate_gene, self))
        return categorize_host_specimen(self, pubmed)


Lassa("Lassa")


def build_blast_db(virus):

    aa_seqs = []
    na_seqs = []
    for s in virus.SEGMENTS:
        with open(virus.reference_folder / f'{s}.gb', "r") as handle:
            for record in SeqIO.parse(handle, "genbank"):
                # When
                na_seqs.append(
                    SeqRecord(record.seq, id=s, description=''))

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

                    if gene not in virus.GENES:
                        continue

                    aa_seqs.append(
                        SeqRecord(Seq(aa.qualifiers['translation'][0]), id=gene, description=''))
                    na_seqs.append(
                        SeqRecord(Seq(aa.location.extract(record.seq)), id=gene, description=''))

    ref_aa_file = virus.reference_folder / f"{virus.name}_RefAAs.fasta"
    with open(ref_aa_file, "w") as output_handle:
        SeqIO.write(aa_seqs, output_handle, "fasta")

    subprocess.run(
        f"makeblastdb -in {ref_aa_file} -dbtype "
        f"prot -out {virus.BLAST_AA_DB_PATH}",
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        shell=True
    )

    ref_na_file = virus.reference_folder / f"{virus.name}_RefNAs.fasta"
    with open(ref_na_file, "w") as output_handle:
        SeqIO.write(na_seqs, output_handle, "fasta")

    subprocess.run(
        f"makeblastdb -in {ref_na_file} -dbtype "
        f"nucl -out {virus.BLAST_NA_DB_PATH}",
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        shell=True
    )


# Provides directions for cleaning the information in the feature table
def process_feature(features_df):
    features_df = translate_bio_term(features_df)
    features_df = get_additional_host_data(features_df)
    features_df['Host'] = features_df['Host2']
    features_df['isolate_source'] = features_df['isolate_source2']

    features_df['Country'] = features_df['country_region'].str.split(
        ":").str[0]

    features_df['Genes'] = features_df['segment_source']

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


def process_gene_list(virus, gene_df):

    for i, row in gene_df.iterrows():
        if str(row['Gene']) not in virus.GENES:
            gene_df.at[i, 'Gene'] = row['hit_name']

    return gene_df


def translate_gene(virus, gene):
    return gene if gene in virus.SEGMENTS else ('NA' if not gene or gene == 'NA' else 'Other')


def categorize_host_specimen(self, pubmed):

    for index, row in pubmed.iterrows():
        host = row['Host'].lower()
        specimen = row['IsolateType'].lower()

        updated_host = []
        updated_specimen = []

        if 'homo sapiens' in host:
            updated_host.append('Homo sapiens')

        for a in ['animal', 'sheep', 'cattle', 'goat', 'mouse', 'boar', 'hare',
                  'livestock', 'cow', 'sheep', 'camel', 'monkey', 'deer', 'buffalo',
                  'rodent', 'serow']:
            if a in host:
                updated_host.append(a.capitalize())

        if 'tick' in host:
            updated_host.append('Ticks')

        if 'tick' in specimen:
            updated_host.append('Ticks')

        if not updated_host and host and host != 'NA'.lower():
            updated_host.append('Other')
            # updated_host.append(host)

        if not updated_host:
            updated_host.append('NA')

        for i in ['serum', 'blood', 'plasma', 'sera']:
            if i in specimen:
                updated_specimen.append('blood')

        for s in ['brain', 'spleen', 'nasal swab']:
            if s in specimen:
                updated_specimen.append(s)

        if not updated_specimen:
            updated_specimen.append('NA')

        pubmed.at[index, 'CleanedHost'] = ' and '.join(
            sorted(list(set(updated_host))))
        pubmed.at[index, 'CleanedSpecimen'] = ' and '.join(
            sorted(list(set(updated_specimen))))

        pubmed['Host'] = pubmed['CleanedHost']
        pubmed['Specimen'] = pubmed['CleanedSpecimen']

    return pubmed
