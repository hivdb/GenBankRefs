import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from .virus import Virus
import subprocess
from functools import partial


class Nipah(Virus):

    @property
    def GENES(self):
        return ['N', 'P', 'M', 'F', 'G', 'L']

    @property
    def pubmed_file(self):
        return self.pubmed_folder / "ReferenceSummary_Dec18.xlsx"

    def build_blast_db(self):
        build_blast_db(self)

    def _process_features(self, features_df):
        return process_features(features_df)

    def process_gene_list(self, gene_df):
        return process_gene_list(self, gene_df)

    def process_pubmed(self, pubmed):
        # pubmed['Gene'] = pubmed['Gene'].apply(partial(translate_pubmed_genes, self))
        return categorize_host_specimen(self, pubmed)


Nipah("Nipah")


def build_blast_db(virus):
    aa_seqs = []
    na_seqs = []
    with open(virus.reference_folder / 'NC_002728.gb', "r") as handle:
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
                    gene = aa.qualifiers[
                        'product'][0].upper().replace(' PROTEIN', '').strip()

                if gene == 'P PHOSPHOPROTEIN':
                    gene = 'P'

                if gene not in virus.GENES:
                    # print(gene)
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


def process_features(features_df):
    features_df = translate_bio_term(features_df)
    features_df = get_additional_host_data(features_df)
    features_df['Host'] = features_df['Host2']
    features_df['isolate_source'] = features_df['isolate_source2']

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


def process_gene_list(virus, gene_df):

    gene_df['Gene'] = gene_df['CDS_NAME'].apply(translate_cds_name)

    for i, row in gene_df.iterrows():
        if str(row['Gene']) not in virus.GENES:
            gene_df.at[i, 'Gene'] = row['hit_name']

    return gene_df


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
    elif cds == 'isolate':
        return ''
    else:
        print(cds)
        return ''


def translate_pubmed_genes(virus, gene):
    return gene if gene in virus.GENES else ('NA' if not gene or gene == 'NA' else 'Other')


def categorize_host_specimen(self, pubmed):

    for index, row in pubmed.iterrows():
        host = row['Host'].lower()
        specimen = row['IsolateType'].lower()

        updated_host = []
        updated_specimen = []

        if 'homo sapiens' in host:
            updated_host.append('Homo sapiens')

        for a in ['bat', 'pig', 'dog']:
            if a in host:
                updated_host.append(a.capitalize())

        if not updated_host and host and host != 'NA'.lower():
            updated_host.append('Other')
            # updated_host.append(host)

        if not updated_host:
            updated_host.append('NA')

        for i in ['serum', 'blood', 'plasma', 'sera']:
            if i in specimen:
                updated_specimen.append('blood')

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

        if not updated_specimen:
            updated_specimen.append('NA')

        pubmed.at[index, 'CleanedHost'] = ' and '.join(
            sorted(list(set(updated_host))))
        pubmed.at[index, 'CleanedSpecimen'] = ' and '.join(
            sorted(list(set(updated_specimen))))

        pubmed['Host'] = pubmed['CleanedHost']
        pubmed['Specimen'] = pubmed['CleanedSpecimen']

    return pubmed
