from Bio import SeqIO
from Bio import Entrez


Entrez.email = "rshafer.stanford.edu"


def fetch_genbank_by_accession(accession):
    # Use Entrez to fetch the GenBank file using the accession number
    handle = Entrez.efetch(db="nucleotide", id=accession,
                           rettype="gb", retmode="text")
    return SeqIO.read(handle, "genbank")


def create_ref_aa_seq(accession_list):
    print("AccessionList: ", accession_list)
    # ref_files = []
    combined_ref_aa_seq = ''
    for acc in accession_list:
        record = fetch_genbank_by_accession(acc)
        # ref_files.append(record)
        for feature in record.features:
            if feature.type == "CDS":
                if 'translation' in feature.qualifiers:
                    protein_seq = feature.qualifiers['translation'][0]
                    combined_ref_aa_seq = combined_ref_aa_seq + protein_seq
                else:
                    print("No translation available for this CDS feature.")
    return combined_ref_aa_seq


def filter_by_taxonomy(record):
    excluded_seq = {}
    taxonomy = record.annotations['taxonomy']
    # print(taxonomy)
    excluded_seq['Accession'] = record.id
    excluded_seq['Taxonomy'] = ', '.join(taxonomy)
    excluded_seq['SeqLen'] = len(record.seq)
    excluded_seq['Organism'] = record.annotations['organism']
    excluded_seq['Description'] = record.description
    return excluded_seq


def is_reference_genome(acc):
    prefix_list = ['NC', 'NG', 'NM', 'NR']
    for item in prefix_list:
        if acc.startswith(item):
            return True
    return False
