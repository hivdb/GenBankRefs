from Bio import SeqIO, pairwise2
from Bio.Seq import Seq
from Bio import Entrez
from Bio.Blast.Applications import NcbiblastpCommandline
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast.Applications import NcbiblastxCommandline
from Bio.Blast import NCBIXML
from multiprocessing import Pool
from pathlib import Path
from xml.parsers.expat import ExpatError
from collections import defaultdict
from skbio.alignment import StripedSmithWaterman

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
            identity = hsp.identities
            percent_identity = (identity / alignment_length) * 100

            seq_cut = query_seq[hsp.query_start - 1: hsp.query_end]
            seq = hsp.query

            _ins = 0 if len(seq) <= len(seq) else len(seq) - len(seq_cut)
            _del = 0 if len(seq) >= len(seq) else len(seq_cut) - len(seq)

            sense = (
                'positive'
                if hsp.sbjct_end >= hsp.sbjct_start
                else 'negative'
            )

            trim_length = len(query_seq) - len(seq)

            blast_result.append({
                'hit_name': alignment.hit_def,
                "e_value": hsp.expect,
                'score': hsp.score,
                "pcnt_id": percent_identity,
                "align_len": alignment_length,
                'blast_name': blast_name,
                'query_seq': query_seq,

                'start': hsp.sbjct_start,
                'stop': hsp.sbjct_end,
                'seq': hsp.query,
                'length': len(hsp.query),
                'raw_start': hsp.query_start,
                'raw_stop': hsp.query_end,
                'num_ins': _ins,
                'num_del': _del,
                'num_N': hsp.query.count('N') if blast_name == 'blastn' else hsp.query.count('X') * 3,
                'trim_length': trim_length,
                'sense': sense,
            })

    Path(input_file).unlink(missing_ok=True)
    Path(output_file).unlink(missing_ok=True)

    return blast_result


def blast_gene(gene, blast_aa_db_path, blast_na_db_path):
    """
        Try blastn, blastp, blastx for detecting the genes or segments of an isolate
        This function will decide the best blast result by alignment length of nucleic acid,
        if it's using blastp, the alignment length will be multiplied by 3.
    """

    blast_na = perform_blast(
        gene['Accession'], gene['Order'],
        gene['NA_raw_seq'], blast_na_db_path,
        func=NcbiblastnCommandline, blast_name='blastn')
    blast_na = [i for i in blast_na if i]
    blast_na = get_best_blast(blast_na, by='align_len')

    if (not blast_na) or (blast_na['sense'] != 'positive'):
        blast_na = perform_blast(
            gene['Accession'], gene['Order'],
            Seq(gene['NA_raw_seq']).reverse_complement(), blast_na_db_path,
            func=NcbiblastnCommandline, blast_name='blastn')
        blast_na = [i for i in blast_na if i]
        blast_na = get_best_blast(blast_na, by='align_len')

    blast_aa = []
    if len(gene['AA_raw_seq']) > 30:
        blast_aa = perform_blast(
            gene['Accession'], gene['Order'],
            gene['AA_raw_seq'], blast_aa_db_path,
            func=NcbiblastpCommandline, blast_name='blastp')

        blast_aa = [i for i in blast_aa if i]

    blast_aa = get_best_blast(blast_aa, by='align_len')

    # blastx deprecated
    # perform_blast(
    #     gene['Accession'], gene['Order'],
    #     gene['NA_raw_seq'], blast_aa_db_path,
    #     func=NcbiblastxCommandline, blast_name='blastx')

    # if 'e_value' not in blast_data:
    #     print(blast_data, idx, len(gene['AASeq']) <= 30)

    new_gene = get_new_gene(gene, blast_na, blast_aa)

    return new_gene


def get_best_blast(blasts, by='align_len'):
    blasts = sorted(blasts, key=lambda x: int(x[by]), reverse=True)

    if not blasts:
        blast = {}
    else:
        blast = blasts[0]

    return blast


def pooled_blast_genes(gene_list, virus_obj, poolsize=20):
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
            (gene, virus_obj.BLAST_AA_DB_PATH, virus_obj.BLAST_NA_DB_PATH)
            for gene in gene_list
        ]
        alignment_result = []
        for count, i in enumerate(pool.starmap(blast_gene, parameters)):
            alignment_result.append(i)

    for i in alignment_result:
        gene_name = virus_obj.translate_cds_name(i['CDS_NAME'])
        if gene_name:
            i['Gene'] = gene_name
            continue

        if i.get('hit_name'):
            i['Gene'] = i['hit_name']

    return alignment_result


def local_align_genes(seq, virus_obj):

    gene_dict = {}  # Dictionary to store gene names and sequences

    with open(virus_obj.ref_na_path) as fasta_file:
        for record in SeqIO.parse(fasta_file, "fasta"):
            gene_dict[record.id] = str(record.seq)

    matched_genes = []
    for gene, ref_seq in gene_dict.items():
        # alignments = pairwise2.align.localms(seq, ref_seq, 2, -3, -5, -2)

        # best_alignment = alignments[0]  # Take the best alignment
        # aligned_seq1, aligned_seq2, align_score, start, end = best_alignment
        # seq1_length = len(aligned_seq1.replace("-", ""))  # Exclude gaps

        # matches = sum(1 for a, b in zip(aligned_seq1, aligned_seq2) if a == b and a != '-' and b != '-')
        # percent_identity = (matches / len(ref_seq)) * 100 if seq1_length > 0 else 0

        query = StripedSmithWaterman(ref_seq)
        alignment = query(seq)

        if alignment.optimal_alignment_score > len(ref_seq) * 0.80:  # 80% similarity threshold
            matched_genes.append({
                'Gene': gene,
                'Alignment Length': len(alignment.aligned_query_sequence),  # remove gaps?
                'Percent Identity': len(alignment.aligned_query_sequence) / len(ref_seq),
                'NA_len': len(alignment.aligned_target_sequence),
                'AA_len': len(alignment.aligned_target_sequence) // 3
            })

    return matched_genes


def detect_additional_genes(
        gene_list, gene_list2, virus_obj, poolsize=20):

    isolate_genes = defaultdict(list)
    for i in gene_list2:
        gene_name = i['Gene']
        if not gene_name:
            continue
        isolate_genes[i['Accession']].append(i['Gene'])

    isolates = [
        i for i in gene_list
        if i['CDS_NAME'] in ('isolate', 'isolate_complete')
    ]

    additional_genes = []

    # Step 1: Run BLAST First
    with Pool(poolsize) as pool:
        parameters = [(isolate, isolate_genes[isolate['Accession']],
                       virus_obj.BLAST_NA_DB_PATH) for isolate in isolates]
        for count, blast_results in enumerate(
                pool.starmap(blast_addi_isolates, parameters)):
            additional_genes.extend(blast_results)

    for i in additional_genes:
        if not i['Gene'] and i.get('hit_name'):
            i['Gene'] = i['hit_name']

    return additional_genes


# Step 2: Run Local Alignment in Parallel (for Isolates Without BLAST Hits)
def detect_gene_by_biopython(isolates, virus_obj):
    additional_genes = []

    for isolate in isolates:
        seq = isolate.get('NA_raw_seq', '')
        # Attempt gene detection using local alignment
        aligned_genes = local_align_genes(seq, virus_obj)

        if aligned_genes:
            for gene in aligned_genes:
                # print(isolate['Accession'], 'new')
                new_gene = {
                    'Accession': isolate['Accession'],
                    'Gene': gene["Gene"],
                    'CDS_NAME': '',
                    'Order': isolate['Order'],
                    'detected_gene': 1,
                    'NA_length': gene["NA_len"],
                    'AA_length': gene["AA_len"],
                    'NA_raw_seq': seq,
                    'align_len': gene["Alignment Length"],
                    'pcnt_id': gene["Percent Identity"]
                }
                additional_genes.append(new_gene)

    # print(len(additional_genes), 'add genes')
    return additional_genes


def blast_addi_isolates(isolate, isolate_genes, blast_na_db_path):

    blast_result = perform_blast(
        isolate['Accession'], isolate['Order'],
        isolate['NA_raw_seq'], blast_na_db_path,
        func=NcbiblastnCommandline, blast_name='blastn')

    additional_genes = []

    for blast in blast_result:
        hit_name = blast['hit_name']
        if hit_name == 'genome':
            continue

        if (blast['sense'] != 'positive'):
            blast = perform_blast(
                isolate['Accession'], isolate['Order'],
                Seq(isolate['NA_raw_seq']).reverse_complement(), blast_na_db_path,
                func=NcbiblastnCommandline, blast_name='blastn')
            blast = [i for i in blast if i]
            blast = get_best_blast(blast, by='align_len')

        if hit_name in isolate_genes:
            continue

        new_gene = get_new_gene(isolate, blast, {})
        new_gene['detected_gene'] = 1

        new_gene['CDS_NAME'] = ''
        additional_genes.append(new_gene)

    return additional_genes


def get_new_gene(gene, blast_na, blast_aa):

    blast_columns = [
        ('hit_name', ''),
        ('e_value', 999),
        ('pcnt_id', 0),
        ('align_len', 0),
    ]

    for key, default in blast_columns:
        if not blast_na and blast_aa:
            gene[key] = blast_aa.get(key, default)
        elif blast_na and not blast_aa:
            gene[key] = blast_na.get(key, default)
        else:
            gene[key] = blast_na.get(key, blast_aa.get(key, default))

    if (
        blast_na.get('hit_name') and
        blast_aa.get('hit_name') and
        blast_na['hit_name'] != blast_aa['hit_name']
            ):
        gene['diff_hit_name'] = 1

    align_columns = [
        ('raw_start', 0),
        ('raw_stop', 0),
        ('seq', ''),
        ('start', 0),
        ('stop', 0),
        ('length', 0),
        ('num_ins', 0),
        ('num_del', 0),
        ('trim_length', 0),  # trim length calculation is hard
        ('sense', ''),
    ]

    for key, default in align_columns:
        gene[f"NA_{key}"] = blast_na.get(key, default)

    for key, default in align_columns:
        gene[f"AA_{key}"] = blast_aa.get(key, default)

    gene['num_N'] = max(blast_aa.get('num_N', 0), blast_na.get('num_N', 0))

    # dont use NA to support AA or vice versa, the alignment need more work
    # if not gene['AA_seq'] and gene['NA_seq']:
    #     gene['AA_start'] = ((gene['NA_start'] - 1) // 3) + 1
    #     gene['AA_stop'] = (gene['NA_stop'] // 3)
    #     # gene['AA_seq'] = Seq(gene['NA_seq']).translate(to_stop=True)
    #     gene['AA_seq'] = Seq(
    #         gene['NA_raw_seq'][gene['NA_raw_start'] - 1: gene['NA_raw_stop']]
    #         ).translate(to_stop=True)
    #     gene['AA_length'] = len(gene['AA_seq'])

    if not gene['NA_seq']:
        gene['NA_blast_failed'] = 1

    if not gene['AA_seq']:
        gene['AA_blast_failed'] = 1

    # if not gene['NA_seq'] and gene['AA_seq']:
    #     gene['NA_start'] = (gene['AA_start'] - 1) * 3 + 1
    #     gene['NA_stop'] = gene['AA_stop'] * 3
    #     gene['NA_seq'] = gene['NA_raw_seq'][gene['NA_start'] - 1: gene['NA_stop']]
    #     gene['NA_length'] = len(gene['NA_seq'])

    if (
            gene['NA_seq'] and
            gene['AA_seq'] and
            gene['NA_length'] != gene['AA_length'] * 3
            ):
        gene['translation_issue'] = 1

    # use na for aa seq
    # gene['AASeq'] = Seq(gene['NASeq']).translate(to_stop=False)

    columns_order = [
        'Accession',
        'Gene',
        'CDS_NAME',
        'Order',
    ]

    raw_columns = [
        ('raw_seq', ''),
        ('raw_length', 0),
    ]

    for key, default in (raw_columns + align_columns):
        columns_order.append(f"AA_{key}")

    for key, default in (raw_columns + align_columns):
        columns_order.append(f"NA_{key}")

    columns_order.append('num_N')

    columns_order.extend([
        'detected_gene',
        "AA_blast_failed",
        'NA_blast_failed',
        'diff_hit_name',
        'translation_issue'
    ])

    for key, default in blast_columns:
        columns_order.append(key)

    missing_columns = [
        i
        for i in gene.keys()
        if i not in columns_order
    ]

    if missing_columns:
        print('Missing Gene data columns:', missing_columns)

    new_gene = {}
    for c in columns_order:
        new_gene[c] = gene.get(c, '')

    return new_gene
