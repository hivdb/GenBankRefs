from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from .virus import Virus
import subprocess
import re

# Define a class for Lassa Virus that inherits from the Virus class
class Lassa(Virus):

    @property
    def SEGMENTS(self):
        return ['L segment', 'S segment']

    @property
    def GENES(self):
        return [
            'N',
            'G',
            'Z',
            'L',
        ]

    @property
    def pubmed_file(self):
        return self.pubmed_folder / "Reference_Summary_Mar25.xlsx"

    @property
    def pubmed_additional_from_gb(self):
        return self.pubmed_folder / "ReferenceSummary_Genbank_Feb19.xlsx"

    # @property
    # def pubmed_genbank_hardlink(self):
    #     return self.pubmed_folder / "Reference_Hardlink_Jan31.xlsx"

    @property
    def pubmed_search_missing(self):
        return self.pubmed_folder / "ReferenceSummary_PubMed_Missing_Feb19.xlsx"

    def build_blast_db(self):
        build_blast_db(self)

    def _process_features(self, features_df):
        return process_features(features_df)

    def process_pubmed(self, pubmed):
        # pubmed['Gene'] = pubmed['Gene'].apply(partial(translate_pubmed_genes, self))
        return categorize_host_specimen(self, pubmed)

    def translate_cds_name(self, cds):
        return translate_cds_name(cds)

    def pick_phylo_sequence(self, genes, picked_genes=['G', 'N', 'L']):
        return super().pick_phylo_sequence(genes, picked_genes)


Lassa("Lassa")


def build_blast_db(virus):

    aa_seqs = []
    na_seqs = []
    for s in virus.SEGMENTS:
        with open(virus.reference_folder / f'{s}.gb', "r") as handle:
            for record in SeqIO.parse(handle, "genbank"):

                # na_seqs.append(
                #     SeqRecord(record.seq, id=s, description=''))

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
                        gene = aa.qualifiers[
                            'product'][0].upper().replace(' PROTEIN', '').strip()

                    if gene == 'NUCLEOPROTEIN':
                        gene = 'N'
                    elif gene == 'GLYCOPROTEIN':
                        gene = 'G'

                    if gene not in virus.GENES:
                        continue

                    aa_seqs.append(
                        SeqRecord(Seq(aa.qualifiers['translation'][0]), id=gene, description=''))

                    na_seq = aa.location.extract(record.seq).upper()
                    if na_seq[-3:] in ['TAG', 'TGA', 'TAA']:
                        na_seq = na_seq[:-3]
                    na_seqs.append(
                        SeqRecord(Seq(na_seq), id=gene, description=''))

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
def process_features(features_df):
    features_df = translate_bio_term(features_df)
    features_df = get_additional_host_data(features_df)
    features_df['host_orig'] = features_df['Host']
    features_df['isolate_source_orig'] = features_df['isolate_source']

    features_df['Host'] = features_df['Host2']
    features_df['isolate_source'] = features_df['isolate_source2']

    features_df['isolate_source'] = features_df['isolate_source'].apply(
            lambda x: x.capitalize() if x.upper() != "NA" else x)

    for i, row in features_df.iterrows():
        if 'Kosovo' in row['IsolateName']:
            features_df.at[i, 'country_region'] = 'Kosovo'
        elif 'China' in row['IsolateName']:
            features_df.at[i, 'country_region'] = 'China'

    return features_df

# Function to standardize biological terms in feature tables
def translate_bio_term(features_df):
    name_map = {
        'Mastomys natalensis': 'rodent',
        'Mastomys erythroleucus': 'rodent',
        'Lophuromys sikapusi': 'rodent',
        'Cavia porcellus (guinea pig)': 'rodent',
        'Mus baoulei': 'rodent',
        'Mastomys': 'rodent',
        'Hylomyscus pamfi': 'rodent',
        'Rattus norvegicus': 'rodent',
        'Mastomys sp.': 'rodent',
        'Mastomys sp': 'rodent',
        'Rodents': 'rodent',
    }

    features_df['Host2'] = features_df['Host']
    features_df['isolate_source2'] = features_df['isolate_source']
    for k, v in name_map.items():
        features_df['Host2'] = features_df['Host2'].str.replace(
            k, v, regex=True)
        features_df['isolate_source2'] = features_df['isolate_source2'].str.replace(
            k, v, regex=True)

    features_df['organism'] = features_df['organism'].str.replace(
        r'.*Mammarenavirus lassaense.*', 'Lassa', case=False, regex=True)
    features_df['organism'] = features_df['organism'].str.replace(
        r'.*Lassa virus.*', 'Lassa', case=False, regex=True)

    return features_df


# Function to classify host and specimen types based on feature data
def get_additional_host_data(features_df):
    blood_specimen = ['blood', 'serum', 'plasma', 'sera']
    organs = [
        'brain',
        'spleen',
        'kidney',
        'liver',
        'lung',
        'tissue',
    ]
    other_specimen = [
        'pleural fluid',
        'urine',
        'csf',
        'breast milk',
        'rectal swab',
        'feces',
    ]
    human_host = ['patient', 'human', 'homo sapiens', 'homon sapiens']
    animal_host = [
        'rodent',
        'mouse',
        # 'rat',
        # 'goat', 'dog',
        # 'lizard'
    ]

    for index, row in features_df.iterrows():

        host = row['Host2'].lower().strip()
        specimen = row['isolate_source2'].lower().strip()

        if not host and not specimen:
            continue

        updated_host = []
        updated_specimen = []

        if any(key in specimen
               for key in human_host) or any(key in host
                                             for key in human_host):
            updated_host.append("Human")
            found_specimen = False
            for key in blood_specimen:
                if key in specimen:
                    updated_specimen.append("Blood")
                    found_specimen = True
            for key in other_specimen:
                if key in specimen:
                    updated_specimen.append(key)
                    found_specimen = True
            for key in organs:
                if key in specimen:
                    updated_specimen.append("Organs")
                    found_specimen = True
            if not found_specimen:
                updated_specimen.append("Human")

        for a in animal_host:
            if a in specimen:
                updated_host.append(a.capitalize())
                updated_specimen.append(a)
            if a in host:
                updated_host.append(a.capitalize())
                updated_specimen.append(a)

        if not updated_host and host:
            match = re.search(r"\(([^)]+)\)", host)
            if match:
                updated_host.append(match.group(1).capitalize())
            else:
                updated_host.append(host.capitalize())  # ['Other']

        if not updated_specimen and specimen:
            # specieman other and NA are the same
            updated_specimen = ['']

        # features_df.at[index, 'Host'] = ",".join(sorted(list(set(updated_host))))
        # features_df.at[index, 'isolate_source'] = ",".join(sorted(list(set(updated_specimen))))
        features_df.at[index, 'Host2'] = ', '.join(
            sorted(list(set(updated_host)))) if updated_host else ''
        features_df.at[index, 'isolate_source2'] = ', '.join(
            sorted(list(set(updated_specimen)))) if updated_specimen else ''

    return features_df


def translate_cds_name(cds):
    name_map = {
        'GPC': 'G',
        'GP': 'G',
        'GLYCOPROTEIN': 'G',
        'GLYCOPROTEIN PRECURSOR': 'G',
        'UNNAMED PRODUCT; GLYCOPROTEIN PRECURSOR (AA 1-490)': 'G',

        'N': 'N',
        'NP': 'N',
        'NUCLEOCAPSID (N)': 'N',
        'UNNAMED PRODUCT; NUCLEOCAPSID (AA 1-570)': 'N',
        'NUCLEOPROTEIN': 'N',
        'NUCLEOCAPSID': 'N',

        'POL': 'L',
        'L POLYMERASE': 'L',
        'RING FINGER': 'L',
        'RING-FINGER': 'L',
        'RNA-DEPENDENT RNA POLYMERASE': 'L',
        'POLYMERASE': 'L',
        'L': 'L',

        'Z': 'Z',

        'UNNAMED PRODUCT': '',
        'HYPOTHETICAL': '',
    }

    for v in name_map.values():
        assert (not v or (v in Virus.get_virus('Lassa').GENES))

    if cds in name_map:
        return name_map[cds]
    elif cds in ('isolate', 'isolate_complete'):
        return ''
    elif not cds:
        return ''
    else:
        print('Missing CDS translation', cds)
        return ''


def translate_pubmed_genes(virus, gene):
    return gene if gene in virus.GENES else ('NA' if not gene or gene == 'NA' else 'Other')


def categorize_host_specimen(self, pubmed):
    """
    This function categorizes host and specimen types from a PubMed dataset.
    It cleans and standardizes host and specimen information based on known categories.
    """

    tissue = ['tissue', 'brain', 'lung', 'spleen', 'kidney', 'liver']
    other_source = [
        'pleural fluid', 'urine', 'csf', 'breast milk', 'rectal swab', 'feces'
    ]
    for index, row in pubmed.iterrows():
        host = row['Host'].lower()
        specimen = row['IsolateType'].lower()

        updated_host = []
        updated_specimen = []

        if 'homo sapiens' in host:
            updated_host.append('Human')
            updated_specimen.append('Human')
            #only update for human
            for s in tissue:
                if s in specimen:
                    updated_specimen.append("Tissue & Organ")
            for s in other_source:
                if s in specimen:
                    updated_specimen.append(s)
            for i in ['serum', 'blood', 'plasma', 'sera']:
                if i in specimen:
                    updated_specimen.append('Blood')

        if 'mice' in host:
            updated_host.append("Mouse")

        for a in ['rodent', 'mouse', 'rat', 'goat', 'dog', 'lizard', 'pig']:
            if a in host:
                updated_host.append(a.capitalize())

        for a in ['cell', 'vero', 'biopsy', 'lab', 'culture', 'recombinant']:
            if a in specimen:
                updated_host.append("Lab Sample")
            if a in host:
                updated_host.append("Lab Sample")

        if not updated_host:
            if host and host != 'NA'.lower():
                updated_host.append(host)
            else:
                updated_host.append('NA')

        if not updated_specimen:
            updated_specimen.append('NA')

        pubmed.at[index, 'CleanedHost'] = ', '.join(
            sorted(list(set(updated_host))))
        pubmed.at[index, 'CleanedSpecimen'] = ', '.join(
            sorted(list(set(updated_specimen))))

    pubmed['Host'] = pubmed['CleanedHost']
    pubmed['Specimen'] = pubmed['CleanedSpecimen']
    pubmed['Specimen'] = pubmed['Specimen'].apply(lambda x: x.capitalize()
                                                  if x.upper() != "NA" else x)

    return pubmed
