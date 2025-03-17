from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from .virus import Virus
import subprocess
from functools import partial


# Define a class for CCHF Virus that inherits from the Virus class
class CCHF(Virus):

    @property
    def SEGMENTS(self):
        return ['L', 'M', 'S']

    @property
    def GENES(self):
        return self.SEGMENTS

    @property
    def comparison_file(self):
        return self.output_dir / f"{self.name}_Combined_11_06a.xlsx"

    @property
    def pubmed_file(self):
        return self.pubmed_folder / "ReferenceSummary_Feb25.xlsx"

    @property
    def pubmed_additional_from_gb(self):
        return self.pubmed_folder / "ReferenceSummary_Genbank_Feb25.xlsx"

    @property
    def pubmed_search_missing(self):
        return self.pubmed_folder / "ReferenceSummary_PubMed_Missing_Mar17.xlsx"

    @property
    def pubmed_genbank_hardlink(self):
        return self.pubmed_folder / "Reference_Hardlink_Jan31.xlsx"

    def build_blast_db(self):
        build_blast_db(self)

    def _process_features(self, features_df):
        return process_features(features_df)

    def process_pubmed(self, pubmed):
        pubmed['Gene'] = pubmed['Gene'].apply(partial(translate_pubmed_genes, self))
        return categorize_host_specimen(self, pubmed)

    def pick_phylo_sequence(self, genes, picked_genes=['S', 'M', 'L']):
        return super().pick_phylo_sequence(genes, picked_genes, coverage_pcnt=0.85)


CCHF("CCHF", full_name="Crimean–Congo hemorrhagic fever virus")


def build_blast_db(virus):

    aa_seqs = []
    na_seqs = []
    for s in virus.SEGMENTS:
        with open(virus.reference_folder / f'{s}.gb', "r") as handle:
            for record in SeqIO.parse(handle, "genbank"):

                # TODO, support all genes, or support segments

                aa_seq = [
                    i
                    for i in record.features
                    if i.type == 'CDS'
                ][0].qualifiers['translation'][0]
                aa_seqs.append(SeqRecord(Seq(aa_seq), id=s, description=''))
                na_seqs.append(
                    SeqRecord(Seq(record.seq), id=s, description=''))

    ref_aa_file = virus.reference_folder / f"{virus.name}_RefAAs.fasta"
    with open(ref_aa_file, "w") as output_handle:
        SeqIO.write(aa_seqs, output_handle, "fasta")

    # Run BLAST database creation for protein sequences
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

    # Run BLAST database creation for nucleotide sequences
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

    return features_df


# Function to standardize biological terms in feature tables
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

# Function to classify host and specimen types based on feature data
def get_additional_host_data(features_df):
    blood_specimen = ['blood', 'serum', 'plasma', 'sera']
    organs = [
        'tissue',
        'brain',
        'spleen',
        'kidney',
        'liver',
        'lung',
    ]
    other_specimen = [
        'nasopharyngeal swab',
        'pleural fluid',
        'urine',
        'csf',
        'breast milk',
        'rectal swab',
        'feces',
    ]
    human_host = ['patient', 'human', 'homo sapiens']
    animal_host = [
        'tick', 'mouse', 'rat', 'jerboa',
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

        if any(key in specimen for key in human_host) or any(key in host for key in human_host):
            updated_host.append("Human")
            found_specimen = False
            # Only get specimen for human hosts
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

        # if animal, simply append host name to specimen, ticks special case
        for a in animal_host:
            if a in specimen:
                if a == 'tick':
                    updated_host.append("Tick")
                    updated_specimen.append("Tick")
                # just put animal for other animals
                else:
                    updated_host.append(a.capitalize())
                    updated_specimen.append("Animal")
            elif a in host:
                if a == 'tick':
                    updated_host.append("Tick")
                    updated_specimen.append("Tick")
                else:
                    updated_host.append(a.capitalize())
                    updated_specimen.append("Animal")

        if not updated_host and host:
            updated_host.append(host)  # ['Other']

        if not updated_specimen and specimen:
            # specimen other and NA are the same
            updated_specimen = ['']

        # features_df.at[index, 'Host'] = ",".join(sorted(list(set(updated_host))))
        # features_df.at[index, 'isolate_source'] = ",".join(sorted(list(set(updated_specimen))))
        features_df.at[index, 'Host2'] = ', '.join(
            sorted(list(set(updated_host)))) if updated_host else ''
        features_df.at[index, 'isolate_source2'] = ', '.join(
            sorted(list(set(updated_specimen)))) if updated_specimen else ''

    return features_df


def translate_pubmed_genes(virus, gene):
    return gene if gene in virus.SEGMENTS else ('NA' if not gene or gene == 'NA' else 'Other')


def categorize_host_specimen(self, pubmed):
    """
    This function categorizes host and specimen types from a PubMed dataset.
    It cleans and standardizes host and specimen information based on known categories.
    """

    for index, row in pubmed.iterrows():
        host = row['Host'].lower()
        specimen = row['IsolateType'].lower()

        updated_host = []
        updated_specimen = []

        # Human, specimen matters
        if 'homo sapiens' in host:
            updated_host.append('Human')
            # only check for human what's the source
            for i in ['serum', 'blood', 'plasma', 'sera']:
                if i in specimen:
                    updated_specimen.append('blood')

            for s in ['brain', 'spleen', 'nasal swab']:
                if s in specimen:
                    updated_specimen.append(s)

        # For ticks specimen & hosts are both
        if 'tick' in host or 'tick' in specimen:
            updated_host.append('Tick')
            updated_specimen.append('Tick')

        # Rest of the animals
        for a in ['animal', 'sheep', 'cattle', 'goat', 'mouse', 'boar', 'hare',
                  'livestock', 'cow', 'sheep', 'camel', 'monkey', 'deer', 'buffalo',
                  'rodent', 'serow']:
            if a in host:
                updated_host.append(a.capitalize())
                updated_specimen.append('Animal')

        for a in ['cell', 'vero', 'biopsy', 'lab', 'culture', 'recombinant']:
            if a in specimen:
                updated_host.append("Lab Sample")
            if a in host:
                updated_host.append("Lab Sample")

        if not updated_host and host and host != 'NA'.lower():
            # updated_host.append('Other')
            updated_host.append(host)

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

    pubmed['Specimen'] = pubmed['Specimen'].apply(
        lambda x: x.capitalize() if x.upper() != "NA" else x)

    return pubmed
