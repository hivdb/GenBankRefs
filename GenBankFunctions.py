from Bio import SeqIO
from Bio import Entrez
import pandas as pd
import numpy as np
import re
import os
from collections import Counter
from Bio.Blast.Applications import NcbiblastpCommandline
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast.Applications import NcbiblastxCommandline
from Bio.Blast import NCBIXML
import Levenshtein
from multiprocessing import Pool

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
    excluded_seq['Accession'] = record.id
    excluded_seq['Taxonomy'] = ', '.join(taxonomy)
    excluded_seq['SeqLen'] = len(record.seq)
    excluded_seq['Organism'] = record.annotations['organism']
    excluded_seq['Description'] = record.description
    return excluded_seq


def perform_blast(idx, query_seq, db_name, func, blast_name):
    """
    Input: sample_seq (str): sequence to compare/blast
        db_name (str), prebuilt db by calling makeblastdb
        output_file (str): tmp file name to store results, cleared each iteration
    Output:
        a dictionary of keys including {e_value, percent_identity, alignment_length, overlap}
            of each alignment

    """
    input_file = f"/tmp/query_{idx}.fasta"
    output_file = f"/tmp/query_{idx}.xml"

    with open(input_file, "w") as fd:
        fd.write(f">query\n{query_seq}\n")

    blastp_cline = func(
        query=input_file,
        db=db_name, outfmt=5, out=output_file)
    stdout, stderr = blastp_cline()

    # Parse the BLAST results
    with open(output_file, "r") as result_handle:
        blast_records = NCBIXML.read(result_handle)

    # Extract statistics (assuming a single hit, adjust as needed)
    # Try blast_records.alignments[0].hsps[0]
    blast_result = []
    for alignment in blast_records.alignments:
        # max_alignment = sorted([
        #     (alignment.length, alignment)
        #     for alignment in blast_records.alignments
        # ], key=lambda x: x[0])[-1][-1]

        # hsp = sorted([
        #     (hsp.score, hsp)
        #     for hsp in max_alignment.hsps
        # ], key=lambda x: x[0])[-1][-1]

        # print(dir(max_alignment), dir(hsp))
        # print(idx, hsp.query, len(hsp.query), hsp.query_start, hsp.query_end)

        for hsp in alignment.hsps:

            e_value = hsp.expect
            alignment_length = hsp.align_length
            identity = hsp.identities
            percent_identity = (identity / alignment_length) * 100
            overlap = hsp.align_length

            blast_result.append({
                'hit_name': alignment.hit_def,
                "e_value": e_value,
                "pcnt_id": percent_identity,
                "align_len": alignment_length * 3 if blast_name == 'blastp' else alignment_length,
                "overlap": overlap,
                'blast_name': blast_name,
            })

    os.remove(input_file)
    os.remove(output_file)

    return blast_result


def is_reference_genome(acc):
    prefix_list = ['NC', 'NG', 'NM', 'NR']
    for item in prefix_list:
        if acc.startswith(item):
            return True
    return False


def blast_sequence(idx, features, virus):

    blast_result = perform_blast(
        features['acc_num'], features['NASeq'], f"{virus}_NA_db",
        func=NcbiblastnCommandline, blast_name='blastn')

    if len(features['AASeq']) > 30:
        blast_result.extend(perform_blast(
            features['acc_num'], features['AASeq'], f"{virus}_AA_db",
            func=NcbiblastpCommandline, blast_name='blastp'))

    blast_result.extend(perform_blast(
        features['acc_num'], features['NASeq'], f"{virus}_AA_db",
        func=NcbiblastxCommandline, blast_name='blastx'))

    # if 'e_value' not in blast_data:
    #     print(blast_data, idx, len(features['AASeq']) <= 30)

    blast_result = sorted([
            b
            for b in blast_result
        ], key=lambda x: int(x['align_len']), reverse=True)

    blast_data = blast_result[0]

    features['hit_name'] = blast_data.get('hit_name', '')
    features['e_value'] = blast_data.get('e_value', '')
    features['pcnt_id'] = blast_data.get('pcnt_id', '')
    features['align_len'] = blast_data.get('align_len', '')
    features['blast_name'] = blast_data.get('blast_name', '')

    features['hit_name_list'] = ', '.join([
        str(i['hit_name'])
        for i in blast_result])
    features['e_value_list'] = ', '.join([
        str(i['e_value'])
        for i in blast_result])
    features['pcnt_id_list'] = ', '.join([
        str(i['pcnt_id'])
        for i in blast_result])
    features['align_len_list'] = ', '.join([
        str(i['align_len'])
        for i in blast_result])
    features['blast_name_list'] = ', '.join([
        str(i['blast_name'])
        for i in blast_result])

    return features


def pooled_blast(features_list, db_name, poolsize=20):

    # features_list = [
    #     i
    #     for i in features_list
    #     if i['acc_num'] in ['FV537249.1', 'FV537248.1']
    # ]

    with Pool(poolsize) as pool:
        parameters = [
            (idx, f, db_name)
            for idx, f in enumerate(features_list)
        ]
        alignment_result = []
        for count, i in enumerate(pool.starmap(blast_sequence, parameters)):
            alignment_result.append(i)

    return alignment_result
