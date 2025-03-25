from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from .virus import Virus
import subprocess

# Define a class for Lassa Virus that inherits from the Virus class
class Nipah(Virus):

    @property
    def GENES(self):
        return ['N', 'P', 'M', 'F', 'G', 'L']

    @property
    def pubmed_file(self):
        return self.pubmed_folder / "ReferenceSummary_Mar25.xlsx"

    @property
    def pubmed_additional_from_gb(self):
        return self.pubmed_folder / "ReferenceSummary_Genbank_Jan13.xlsx"

    # @property
    # def pubmed_genbank_hardlink(self):
    #     return self.pubmed_folder / "Reference_Hardlink_Jan31.xlsx"

    @property
    def pubmed_search_missing(self):
        return self.pubmed_folder / "ReferenceSummary_PubMed_Missing_Mar17.xlsx"

    def build_blast_db(self):
        build_blast_db(self)

    def _process_features(self, features_df):
        return process_features(features_df)

    def process_pubmed(self, pubmed):
        # pubmed['Gene'] = pubmed['Gene'].apply(partial(translate_pubmed_genes, self))
        return categorize_host_specimen(self, pubmed)

    def pick_phylo_sequence(self, genes, picked_genes=['N', 'L']):
        return super().pick_phylo_sequence(genes, picked_genes)

    def translate_cds_name(self, cds):
        return translate_cds_name(cds)


Nipah("Nipah")


def build_blast_db(virus):
    aa_seqs = []
    na_seqs = []
    with open(virus.reference_folder / 'NC_002728.gb', "r") as handle:
        for record in SeqIO.parse(handle, "genbank"):
            # na_seqs.append(
            #     SeqRecord(record.seq, id='genome', description=''))

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

                if gene == 'P PHOSPHOPROTEIN':
                    gene = 'P'

                if gene not in virus.GENES:
                    # print(gene)
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
        if 'NV/MY' in row['IsolateName']:
            features_df.at[i, 'country_region'] = 'Malaysia'
        elif 'MALAYSIA' in row['IsolateName']:
            features_df.at[i, 'country_region'] = 'Malaysia'

    # for i, row in features_df.iterrows():
    #     print(row)
    #     if int(row['SeqLength']) > 17000:
    #         features_df.at[i, 'Genes'] = 'genome'

    return features_df


def translate_bio_term(features_df):

    name_map = {
        r'Pteropus\s\w+\b': 'bat',
        'Homo sapiens; male': 'Homo sapiens',
        'Sus scrofa domesticus': 'Pig',
        'Canis lupus familiaris': 'Dog',
        'Sus scrofa (pig)': 'Pig',
        'oropharyngeal swab': 'throat swab',
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
    organs = [
        'kidney',
        'liver',
        'spleen',
        'lung',
        'heart',
        'brain',
    ]
    other_specimen = [
        'breast milk',
        'csf',
        'intestine',
        'urine',
        'throat swab',
    ]
    human_host = ['patient', 'human', 'homo sapiens']
    animal_host = ['bat', 'pig', 'dog', 'flying foxes']

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
                if a == 'flying foxes':
                    updated_host.append("Bat")
                else:
                    updated_host.append(a.capitalize())



        if not updated_host and host:
            updated_host = ['Other']

        if not updated_specimen and specimen:
            # specieman other and NA are the same
            updated_specimen = ['']

        # features_df.at[index, 'Host'] = ",".join(sorted(list(set(updated_host))))
        # features_df.at[index, 'isolate_source'] = ",".join(sorted(list(set(updated_specimen))))
        features_df.at[index, 'Host2'] = ' and '.join(
            sorted(list(set(updated_host)))) if updated_host else ''
        features_df.at[index, 'isolate_source2'] = ' and '.join(
            sorted(list(set(updated_specimen)))) if updated_specimen else ''

    return features_df


def translate_cds_name(cds):
    name_map = {
        'RNA-directed RNA polymerase': 'L',
        'POLYMERASE': 'L',
        'large polymerase': 'L',
        'RNA polymerase': 'L',
        'L': 'L',

        'N': 'N',
        'NUCLEOCAPSID': 'N',
        'NUCLEOPROTEIN': 'N',

        'FUSION': 'F',
        'F': 'F',

        'GLYCOPROTEIN': 'G',
        'G': 'G',
        'ATTACHMENT GLYCOPROTEIN': 'G',

        'MATRIX': 'M',
        'M': 'M',

        'P': 'P',
        'P/V/M/C': 'P',
        'P/V/C': 'P',
        'P/V/W/C': 'P',
        'V': 'P',
        'W': 'P',
        'C': 'P',
        'P PHOSPHOPROTEIN': 'P',
        'PHOSPHOPROTEIN': 'P',

    }

    for v in name_map.values():
        assert (v in Virus.get_virus('Nipah').GENES)

    if cds in name_map:
        return name_map[cds]
    elif cds in ('isolate', 'isolate_complete'):
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

    for index, row in pubmed.iterrows():
        host = row['Host'].lower().strip()
        specimen = row['IsolateType'].lower().strip()

        updated_host = []
        updated_specimen = []

        # Only collect those specimen for human
        if 'homo sapiens' in host:
            updated_host.append('Human')
            for s in [
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
            ]:
                if s in specimen:
                    updated_specimen.append(s)

        for a in ['bat', 'pig', 'dog']:
            if a in host:
                updated_host.append(a.capitalize())

        if 'flying fox' in host:
            updated_host.append('Bat')

        for a in ['cell', 'vero', 'biopsy', 'lab', 'culture', 'recombinant']:
            if a in specimen:
                updated_host.append("Lab Sample")
            if a in host:
                updated_host.append("Lab Sample")

        # if not updated_host and host and host != 'NA'.lower():
        #     # updated_host.append('Other')
        #     updated_host.append(host)

        if not updated_host:
            if host and host != 'NA'.lower():
                updated_host.append(host.capitalize())
            else:
                updated_host.append('NA')

        for i in ['serum', 'blood', 'plasma', 'sera']:
            if i in specimen:
                updated_specimen.append('blood')

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
