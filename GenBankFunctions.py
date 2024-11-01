from Bio import SeqIO
from Bio import Entrez
import pandas as pd
import numpy as np
import re
import os
from collections import Counter
from Bio.Blast.Applications import NcbiblastpCommandline
from Bio.Blast import NCBIXML
import Levenshtein
from multiprocessing import Pool

Entrez.email = "rshafer.stanford.edu"


def fetch_genbank_by_accession(accession):
    # Use Entrez to fetch the GenBank file using the accession number
    handle = Entrez.efetch(db="nucleotide", id=accession, rettype="gb", retmode="text")
    return SeqIO.read(handle, "genbank")


def create_ref_aa_seq(accession_list):
    print("AccessionList: ", accession_list)
    #ref_files = []
    combined_ref_aa_seq = ''
    for acc in accession_list:
        record = fetch_genbank_by_accession(acc)
        #ref_files.append(record)
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
    excluded_seq['Accession'] = record.id
    excluded_seq['Taxonomy'] = ', '.join(taxonomy)
    excluded_seq['SeqLen'] = len(record.seq)
    excluded_seq['Organism'] = record.annotations['organism']
    excluded_seq['Description'] = record.description
    return excluded_seq 



def perform_blastp(idx, sample_seq, db_name):
    """
    Input: sample_seq (str): sequence to compare/blast
        db_name (str), prebuilt db by calling makeblastdb
        output_file (str): tmp file name to store results, cleared each iteration
    Output:
        a dictionary of keys including {e_value, percent_identity, alignment_length, overlap}
            of each alignment

    """
    with open(f"sample{idx}.fasta", "w") as sample_file:
        sample_file.write(f">sample_seq\n{sample_seq}\n")

    output_file = f"sample{idx}.xml"

    # Run BLASTP with the sample sequence against the reference database
    blastp_cline = NcbiblastpCommandline(query=f"sample{idx}.fasta", db=db_name, outfmt=5, out=output_file)
    stdout, stderr = blastp_cline()

    # Parse the BLAST results
    with open(output_file, "r") as result_handle:
        blast_records = NCBIXML.read(result_handle)

    # Extract statistics (assuming a single hit, adjust as needed)
    # Try blast_records.alignments[0].hsps[0]
    result = {}
    for alignment in blast_records.alignments:
        for hsp in alignment.hsps:  # High-scoring segment pairs
            e_value = hsp.expect
            alignment_length = hsp.align_length
            identity = hsp.identities
            percent_identity = (identity / alignment_length) * 100
            overlap = hsp.align_length

            result = {
                "e_value": e_value,
                "pcnt_id": percent_identity,
                "align_len": alignment_length,
                "overlap": overlap
            }
            break

    os.remove(f"sample{idx}.fasta")
    os.remove(output_file)

    return result


def is_reference_genome(acc):
    prefix_list = ['NC', 'NG', 'NM', 'NR']
    for item in prefix_list:
        if acc.startswith(item):
            return True
    return False


def blast_sequence(idx, features, db_name):

    if len(features['_sample_seq']) <= 30:
        return features

    blast_data = perform_blastp(idx, features['_sample_seq'], db_name)
        # print(blast_data)
    features['e_value'] = blast_data['e_value']
    features['pcnt_id'] = blast_data['pcnt_id']
    features['align_len'] = blast_data['align_len']

    return features


def pooled_blast(features_list, db_name, poolsize=20):

    with Pool(poolsize) as pool:
        parameters = [
            (idx, f, db_name)
            for idx, f in enumerate(features_list)
        ]
        feature_list = []
        for count, i in enumerate(pool.starmap(blast_sequence, parameters)):
            feature_list.append(i)

    return feature_list
