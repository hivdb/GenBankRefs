from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Data.CodonTable import TranslationError

from Bio.Blast.Applications import NcbiblastpCommandline
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast.Applications import NcbiblastxCommandline
from Bio.Blast import NCBIXML
from multiprocessing import Pool
from pathlib import Path
from xml.parsers.expat import ExpatError
from collections import defaultdict
from skbio.alignment import StripedSmithWaterman
from bioinfo import dump_fasta
from bioinfo import load_fasta
import subprocess


def check_genbank_coding_seq(gene_list, virus_obj, poolsize=20):
    """
        For speeding up blast running, using multiprocessing method
        Input:
            gene_list: a list of items, for each one contains NA or AA sequences
            db_name: blast database name
            poolsize: depending on how many cores the CPU has, you can adjust this number for running multiple blast program.
    """

    gene_list = [
        i
        for i in gene_list
        if i['CDS_NAME'] not in ('isolate', 'isolate_complete')
    ]

    with Pool(poolsize) as pool:
        parameters = [
            (gene,
             virus_obj)
            for gene in gene_list
        ]
        alignment_result = []
        for count, i in enumerate(
                pool.imap_unordered(blast_gene, parameters)):
            if count % 10 == 0:
                print(f'\rProgress: {count}/{len(parameters)}', end='', flush=True)
            alignment_result.append(i)

    return alignment_result


def blast_gene(args):
    """
        Try blastn, blastp, blastx for detecting the genes or segments of an isolate
        This function will decide the best blast result by alignment length of nucleic acid,
        if it's using blastp, the alignment length will be multiplied by 3.
    """
    gene, virus_obj = args

    blast_na = perform_blast(
        gene['Accession'],
        gene['Order'],
        gene['NA_raw_seq'],
        virus_obj.BLAST_NA_DB_PATH,
        func=NcbiblastnCommandline,
        blast_name='blastn')
    blast_na = [i for i in blast_na if i]
    blast_na = get_best_blast(blast_na, by='align_len')

    if (not blast_na) or (blast_na['sense'] != 'positive'):
        blast_na = perform_blast(
            gene['Accession'],
            gene['Order'],
            Seq(gene['NA_raw_seq']).reverse_complement(),
            virus_obj.BLAST_NA_DB_PATH,
            func=NcbiblastnCommandline,
            blast_name='blastn')
        blast_na = [i for i in blast_na if i]
        blast_na = get_best_blast(blast_na, by='align_len')

    blast_aa = {}
    if len(gene['AA_raw_seq']) > 30:
        blast_aa = perform_blast(
            gene['Accession'],
            gene['Order'],
            gene['AA_raw_seq'],
            virus_obj.BLAST_AA_DB_PATH,
            func=NcbiblastpCommandline,
            blast_name='blastp')

        blast_aa = [i for i in blast_aa if i]
        blast_aa = get_best_blast(blast_aa, by='align_len')

    # blastx deprecated
    # perform_blast(
    #     gene['Accession'], gene['Order'],
    #     gene['NA_raw_seq'], blast_aa_db_path,
    #     func=NcbiblastxCommandline, blast_name='blastx')

    blast_result = get_blast_result(gene, blast_na, blast_aa)

    # Update Gene, hit_method, NA_raw_seq

    gene_name = virus_obj.translate_cds_name(blast_result['CDS_NAME'])
    if gene_name and gene_name in virus_obj.GENES:
        blast_result['Gene'] = gene_name
    elif blast_result.get('hit_name'):
        blast_result['Gene'] = blast_result['hit_name']

    blast_result['hit_method'] = 'coding_seq'

    # print(new_gene['Accession'])

    return blast_result


def perform_blast(acc, order, query_seq, db_name, func, blast_name):
    """
    Input: sample_seq (str): sequence to compare/blast
        db_name (str), prebuilt db by calling makeblastdb
        output_file (str): tmp file name to store results, cleared each iteration
    Output:
        a dictionary of keys including {e_value, percent_identity, alignment_length}
            of each alignment

    """
    input_file = f"/tmp/query_{acc}_{order}.fasta"
    output_file = f"/tmp/query_{acc}_{order}.xml"

    find_db = False
    for i in Path(db_name.parent).resolve().iterdir():
        if i.stem == db_name.stem:
            find_db = True
            break

    if not find_db:
        return []

    if not query_seq.strip():
        return []

    with open(input_file, "w") as fd:
        fd.write(f">query\n{query_seq}\n")

    blastp_cline = func(
        query=input_file,
        db=db_name, outfmt=5, out=output_file)
    stdout, stderr = blastp_cline()

    if not Path(output_file).exists():
        return []

    # Parse the BLAST results
    with open(output_file, "r") as result_handle:
        try:
            blast_records = NCBIXML.read(result_handle)
        except (ExpatError, ValueError):
            return []

    # Extract statistics (assuming a single hit, adjust as needed)
    blast_result = []
    for alignment in blast_records.alignments:

        # sbjct is the reference
        # Query is the query sequence

        for hsp in alignment.hsps:
            alignment_length = hsp.align_length
            # identity = hsp.identities
            # percent_identity = (identity / alignment_length) * 100

            # seq_cut = query_seq[hsp.query_start - 1: hsp.query_end]
            # seq = hsp.query

            # _ins = 0 if len(seq) <= len(seq) else len(seq) - len(seq_cut)
            # _del = 0 if len(seq) >= len(seq) else len(seq_cut) - len(seq)

            sense = (
                'positive'
                if hsp.sbjct_end >= hsp.sbjct_start
                else 'negative'
            )

            # trim_length = len(query_seq) - len(seq)

            blast_result.append({
                'blast_name': blast_name,
                'hit_name': alignment.hit_def,
                'e_value': hsp.expect,
                'score': hsp.score,
                'sense': sense,

                "align_len": alignment_length,

                # "pcnt_id": int(percent_identity),
                # 'query_seq': query_seq,

                # 'start': hsp.sbjct_start,
                # 'stop': hsp.sbjct_end,
                # 'seq': hsp.query,
                # 'length': len(hsp.query),
                # 'raw_start': hsp.query_start,
                # 'raw_stop': hsp.query_end,
                # 'num_ins': _ins,
                # 'num_del': _del,
                # 'num_N': hsp.query.count('N') if blast_name == 'blastn' else hsp.query.count('X') * 3,
                # 'trim_length': trim_length,
            })

    Path(input_file).unlink(missing_ok=True)
    Path(output_file).unlink(missing_ok=True)

    return blast_result


def get_best_blast(blasts, by='align_len'):
    blasts = sorted(blasts, key=lambda x: int(x[by]), reverse=True)

    if not blasts:
        blast = {}
    else:
        blast = blasts[0]

    return blast


def detect_non_annot_gene_by_blast(
        gene_list, checked_gene_list, virus_obj, poolsize=20):

    exist_genes = defaultdict(list)
    for i in checked_gene_list:
        gene_name = i['Gene']
        if not gene_name:
            continue
        exist_genes[i['Accession']].append(gene_name)

    isolates = [
        i for i in gene_list
        if i['CDS_NAME'] in ('isolate', 'isolate_complete')
    ]

    additional_genes = []

    # Step 1: Run BLAST First
    with Pool(poolsize) as pool:
        parameters = [(
            isolate,
            exist_genes[isolate['Accession']],
            virus_obj)
            for isolate in isolates]
        for count, blast_results in enumerate(
                pool.imap_unordered(detect_gene_by_blast, parameters)):
            additional_genes.extend(blast_results)

    return additional_genes


def detect_gene_by_blast(args):

    isolate, exist_genes, virus_obj = args

    blast_result_list = perform_blast(
        isolate['Accession'],
        isolate['Order'],
        isolate['NA_raw_seq'],
        virus_obj.BLAST_NA_DB_PATH,
        func=NcbiblastnCommandline, blast_name='blastn')

    additional_genes = []

    for blast in blast_result_list:
        if blast['hit_name'] not in virus_obj.GENES:
            continue

        if (blast['sense'] != 'positive'):
            blast = perform_blast(
                isolate['Accession'],
                isolate['Order'],
                Seq(isolate['NA_raw_seq']).reverse_complement(),
                virus_obj.BLAST_NA_DB_PATH,
                func=NcbiblastnCommandline,
                blast_name='blastn')

            blast = [i for i in blast if i]
            blast = get_best_blast(blast, by='align_len')

        if blast['hit_name'] in exist_genes:
            continue

        blast_result = get_blast_result(isolate, blast, {})

        blast_result['CDS_NAME'] = ''
        blast_result['Gene'] = blast_result['hit_name']

        blast_result['hit_method'] = 'blast'

        additional_genes.append(blast_result)

    return additional_genes


def detect_non_annot_gen_by_biopython(isolates, virus_obj):
    additional_genes = []

    for isolate in isolates:
        seq = isolate['NA_raw_seq']
        # Attempt gene detection using local alignment
        matched_genes = local_align_genes(seq, virus_obj)

        for gene in matched_genes:
            # print(isolate['Accession'], 'new')
            new_gene = {
                'Accession': isolate['Accession'],
                'Gene': gene["Gene"],
                'CDS_NAME': '',
                'Order': isolate['Order'],
                'hit_method': 'biopython',
                'AA_raw_seq': '',
                'NA_raw_seq': seq,
            }

            additional_genes.append(new_gene)

    return additional_genes


def local_align_genes(seq, virus_obj, sim_threshold=0.80):

    matched_genes = []
    for gene, ref_seq in virus_obj.ref_na_gene_map.items():
        query = StripedSmithWaterman(ref_seq)
        alignment = query(seq)

        if alignment.optimal_alignment_score > len(ref_seq) * sim_threshold:

            start = min(alignment.query_begin + 1, alignment.query_end + 1)
            stop = max(alignment.query_begin + 1, alignment.query_end + 1)

            if start > stop:
                seq = seq[::-1]
            matched_genes.append({
                'Gene': gene,
                'NA_raw_seq': seq,
                # 'Alignment Length': len(alignment.aligned_query_sequence),  # remove gaps?
                # 'Percent Identity': len(alignment.aligned_query_sequence) / len(ref_seq),
                # 'NA_seq': alignment.aligned_target_sequence,
                # 'NA_len': len(alignment.aligned_target_sequence),
                # 'AA_len': len(alignment.aligned_target_sequence) // 3,
                # 'NA_start': start,
                # 'NA_stop': stop,
            })

    return matched_genes


def get_blast_result(gene, blast_na, blast_aa):

    columns_order = [
        'Accession',
        'Gene',
        'CDS_NAME',
        'Order',
        'AA_raw_seq',
        'AA_raw_length',
        'NA_raw_seq',
        'NA_raw_length',
    ]

    result = {}
    for c in columns_order:
        result[c] = gene.get(c, '')

    blast_columns = [
        ('hit_name', ''),

        ('blast_name', ''),
        # ('e_value', 999),
        # ('score', 0),
        ('sense', ''),
        # ('pcnt_id', 0),
        # ('align_len', 0),
    ]

    for key, default in blast_columns:
        if not blast_na and blast_aa:
            result[key] = blast_aa.get(key, default)
        elif blast_na and not blast_aa:
            result[key] = blast_na.get(key, default)
        else:
            result[key] = blast_na.get(key, blast_aa.get(key, default))

    if blast_na:
        result['sense'] = blast_na['sense']

    # missing_columns = [
    #     i
    #     for i in gene.keys()
    #     if i not in columns_order
    # ]

    # if missing_columns:
    #     print('Missing Gene data columns:', missing_columns)

    # fix reverse complement
    na_raw_seq = result['NA_raw_seq']
    if result['sense'] == 'negative':
        na_raw_seq = Seq(result['NA_raw_seq']).reverse_complement()

    if na_raw_seq[-3:] in ['TAG', 'TGA', 'TAA']:
        na_raw_seq = na_raw_seq[:-3]

    result['NA_raw_seq'] = na_raw_seq
    result['NA_raw_length'] = len(na_raw_seq)

    return result


def align_genes(virus, genes_df, poolsize=20):
    with Pool(poolsize) as pool:
        parameters = [
            (row, virus)
            for idx, row in genes_df.iterrows()
        ]
        alignment_result = []
        for count, i in enumerate(
                pool.imap_unordered(align_gene_seq, parameters)):
            if count % 10 == 0:
                print(f'\rProgress: {count}/{len(parameters)}', end='', flush=True)
            alignment_result.append(i)

        alignment_result.sort(key=lambda x: x['SeqID'])

        return alignment_result


def align_gene_seq(args):
    row, virus = args

    gene_name = row['Gene']
    ref_na = virus.ref_na_gene_map[gene_name]
    ref_aa = virus.ref_aa_gene_map[gene_name]

    # query is the references

    query = StripedSmithWaterman(ref_na)

    alignment = query(row['NA_raw_seq'])

    row['NA_seq'] = alignment.aligned_target_sequence
    row['NA_length'] = len(alignment.aligned_target_sequence)
    row['NA_start'] = alignment.query_begin + 1
    row['NA_stop'] = alignment.query_end + 1
    row['NA_num_ins'] = alignment.aligned_query_sequence.count('-')
    row['NA_num_del'] = alignment.aligned_target_sequence.count('-')
    row['num_N'] = alignment.aligned_target_sequence.count('N')
    row['align_len'] = len(alignment.aligned_target_sequence)
    row['pcnt_id'] = count_pcnt_identity(
        alignment.aligned_query_sequence,
        alignment.aligned_target_sequence,
        len(ref_na)
    )

    row['AA_seq'] = row['AA_raw_seq']
    row['AA_length'] = len(row['AA_raw_seq'])

    if not row['AA_seq']:
        try:
            trans_seq = Seq(alignment.aligned_target_sequence)
            trans_seq = trans_seq[:len(trans_seq) - (len(trans_seq) % 3)]
            aa_seq = Seq(trans_seq).translate()
            row['AA_seq'] = str(aa_seq)
            row['AA_length'] = len(aa_seq)
        except TranslationError:
            row['AA_seq'] = ''
            row['AA_length'] = 0
            row['translation_error'] = 1

    if row['AA_seq']:
        query = StripedSmithWaterman(ref_aa)
        alignment = query(row['AA_seq'])
        row['AA_start'] = alignment.query_begin + 1
        row['AA_stop'] = alignment.query_end + 1
        row['AA_num_ins'] = alignment.aligned_query_sequence.count('-')
        row['AA_num_del'] = alignment.aligned_target_sequence.count('-')
        row['AA_num_stop'] = alignment.aligned_target_sequence.count('*')
    return row


def count_pcnt_identity(aligned_ref, aligne_seq, ref_length):
    same = [
        1
        for i, j in zip(aligned_ref, aligne_seq)
        if i == j
    ]
    return int(len(same) * 100 / ref_length)
