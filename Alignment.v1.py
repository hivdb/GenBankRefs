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

        exist_genes.append(blast['hit_name'])

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

    if (len(ref_na) % 3) != 0:
        print(gene_name, 'has frameshift in reference!')

    # query is the references

    query = StripedSmithWaterman(ref_na)

    alignment = query(row['NA_raw_seq'])
    row['cigar'] = alignment.cigar

    aligned_ref = alignment.aligned_query_sequence
    aligned_seq = alignment.aligned_target_sequence

    # if row['Accession'] == 'MK575070' and row['Gene'] == 'L':
    #     print(len([1 for i, j in zip(aligned_seq, aligned_ref) if i != j]), 'www', alignment.cigar)

    # TODO, pre adjust will change AA alignment
    # aligned_ref, aligned_seq = adjust_alignment(aligned_ref, aligned_seq)

    # TODO add prefix or suffix of seq, will make the alignment a bit better.
    (
        aligned_ref_codon, aligned_seq_codon,
        na_start, na_stop) = get_codon_alignment(
            ref_na, ref_aa,
            row['NA_raw_seq'],
            aligned_ref, aligned_seq,
            alignment, adjacent_window=10)

    aligned_ref = ''.join(aligned_ref_codon)
    aligned_seq = ''.join(aligned_seq_codon)

    # aligned_ref, aligned_seq, na_start, na_stop = get_new_aligned_ref_seq(
    #     ref_na, aligned_ref_codon, aligned_seq_codon
    # )

    row['NA_seq'] = aligned_seq
    row['NA_length'] = len(aligned_seq)
    row['NA_start'] = na_start
    row['NA_stop'] = na_stop

    row['NA_num_ins'] = aligned_ref.count('-')
    row['NA_num_del'] = aligned_seq.count('-')
    row['NA_ins_del_diff_3'] = abs(row['NA_num_ins'] - row['NA_num_del']) % 3
    row['NA_num_N'] = aligned_seq.count('N')

    row['align_len'] = len(aligned_ref.replace('-', ''))
    row['pcnt_id'] = count_pcnt_identity(
        aligned_ref,
        aligned_seq,
        len(ref_na)
    )

    aa_start = (na_start - 1) // 3 + 1
    row = translate_aligned_codon(
        row, ref_aa, aa_start, aligned_ref_codon, aligned_seq_codon)

    return row


# def get_new_aligned_ref_seq(ref_na, aligned_ref_codon, aligned_seq_codon):

#     aligned_ref = ''.join(aligned_ref_codon)
#     aligned_seq = ''.join(aligned_seq_codon)

#     p1 = get_pos_pre(aligned_seq, '-')
#     p2 = get_pos_post(aligned_seq, '-')

#     assert (len(aligned_ref) == len(aligned_seq))

#     aligned_ref = aligned_ref[p1: p2]
#     aligned_seq = aligned_seq[p1: p2]
#     start = p1 + 1
#     stop = len(ref_na) - (len(aligned_seq) - p2)
#     return aligned_ref, aligned_seq, start, stop


def count_pcnt_identity(aligned_ref, aligne_seq, ref_length):
    same = [
        1
        for i, j in zip(aligned_ref, aligne_seq)
        if i == j
    ]
    return int(len(same) * 100 / ref_length)


def adjust_alignment(ref, seq, window_size=30):
    # TODO auto guess window size, by diff '-' % 3 == 0
    # This method is not an ideal one because
    # 1. it blindly adjust the window
    # 2. it didn't consider the codon, amino acid alignment.

    new_ref = []
    mid_ref = []
    new_seq = []
    mid_seq = []
    num_ins = 0

    for r, s in zip(ref, seq):
        if r != '-' and s != '-':
            if num_ins == 0:
                new_ref.append(r)
                new_seq.append(s)
            else:
                mid_ref.append(r)
                mid_seq.append(s)
        elif r == '-' and s == '-':
            continue
        elif r == '-':
            num_ins += 1
            mid_ref.append(r)
            mid_seq.append(s)
            if (num_ins % 3) == 0:
                mid_ref = ''.join(mid_ref).replace('-', '')
                mid_seq = ''.join(mid_seq).replace('-', '')
                if num_ins == 0:
                    new_ref.append(mid_ref)
                    new_seq.append(mid_seq)
                elif num_ins > 0:
                    new_ref.append(mid_ref + '-' * num_ins)
                    new_seq.append(mid_seq)
                elif num_ins < 0:
                    new_ref.append(mid_ref)
                    new_seq.append(mid_seq + '-' * (-num_ins))
                mid_ref = []
                mid_seq = []
                num_ins = 0
                continue
        elif s == '-':
            num_ins -= 1
            mid_ref.append(r)
            mid_seq.append(s)
            if (num_ins % 3) == 0:
                mid_ref = ''.join(mid_ref).replace('-', '')
                mid_seq = ''.join(mid_seq).replace('-', '')
                if num_ins == 0:
                    new_ref.append(mid_ref)
                    new_seq.append(mid_seq)
                elif num_ins > 0:
                    new_ref.append(mid_ref + '-' * num_ins)
                    new_seq.append(mid_seq)
                elif num_ins < 0:
                    new_ref.append(mid_ref)
                    new_seq.append(mid_seq + '-' * (-num_ins))
                mid_ref = []
                mid_seq = []
                num_ins = 0
                continue

        # Not to extend forever
        # TODO: any better way?
        if len(mid_ref) >= window_size:
            mid_ref = ''.join(mid_ref)
            mid_seq = ''.join(mid_seq)
            new_ref.append(mid_ref)
            new_seq.append(mid_seq)
            mid_ref = []
            mid_seq = []
            num_ins = 0

    if mid_ref:
        new_ref.append(''.join(mid_ref))
        new_seq.append(''.join(mid_seq))

    new_ref = ''.join(new_ref)
    new_seq = ''.join(new_seq)

    # TODO
    # if new_ref.count('-') == new_seq.count('-'):
    #     new_ref = new_ref.replace('-', '')
    #     new_seq = new_seq.replace('-', '')

    assert (new_ref.replace('-', '') == ref.replace('-', ''))
    assert (new_seq.replace('-', '') == seq.replace('-', ''))

    assert (len(new_ref) == len(new_seq))

    return new_ref, new_seq


def get_codon_alignment(
        ref_na, ref_aa,
        seq_na,
        aligned_ref, aligned_seq,
        alignment,
        adjacent_window=10):

    na_start = alignment.query_begin
    na_stop = alignment.query_end

    pre_ref = ref_na[:alignment.query_begin]

    chop_length = len(pre_ref) % 3
    if chop_length:
        chop_length = 3 - chop_length
        aligned_ref = aligned_ref[chop_length:]
        aligned_seq = aligned_seq[chop_length:]
        na_start += chop_length

    # remove tail and head that are not codon, in case the query seq missing
    # tail and head.
    # although build the full sequence is an option, the tail and head codon
    # in cds won't be changed

    ref_codon_list = get_ref_codon_list(aligned_ref)

    aligned_ref = ''.join(ref_codon_list)
    aligned_seq = aligned_seq[: len(aligned_ref)]
    na_stop = na_start + len(ref_codon_list) * 3 - 1
    assert (len(aligned_ref) == len(aligned_seq))

    seq_codon_list = []
    cursor = 0
    for codon in ref_codon_list:
        seq_codon = aligned_seq[cursor: cursor + len(codon)]
        cursor += len(codon)
        seq_codon_list.append(seq_codon)

    if adjacent_window:
        seq_codon_list = adj_adjacent_seq_codon_list(
            ref_codon_list, seq_codon_list, adjacent_window)
        ref_codon_list, seq_codon_list = adjust_codon_tail_del(
            ref_codon_list, seq_codon_list)

    # assert (len(''.join(ref_codon_list)) == len(aligned_ref))
    assert (len(ref_codon_list) == len(seq_codon_list))
    # TODO: issue, the frameshift of references translation is not considered,
    # should be reflexted in codon alignment
    assert (len(''.join(ref_codon_list)) == len(''.join(seq_codon_list)))
    # assert (len(ref_codon_list) == len(ref_aa))

    return ref_codon_list, seq_codon_list, na_start + 1, na_stop + 1


def get_ref_codon_list(aligned_ref):
    ref_codon_list = []
    ref_codon = []

    for idx, r in enumerate(aligned_ref):
        ref_codon.append(r)
        codon_core = ''.join(ref_codon).replace('-', '')
        if len(codon_core) == 3:
            ref_codon = ''.join(ref_codon)

            codon_pre = ref_codon[:get_pos_pre(ref_codon, '-')]
            this_codon = ref_codon[get_pos_pre(ref_codon, '-'):]
            # codon_post = ref_codon[get_pos_post(ref_codon, '-'):]
            if codon_pre:
                if ref_codon_list:
                    ref_codon_list[-1] += codon_pre
                else:
                    this_codon = ref_codon

            ref_codon_list.append(this_codon)

            ref_codon = []

    return ref_codon_list


def adj_adjacent_seq_codon_list(ref_codon_list, seq_codon_list, adjacent_window):
    adj_adjacent = []
    cursor = 0
    while cursor < len(seq_codon_list):
        codon = seq_codon_list[cursor]
        # Normal codon
        if '-' not in codon and len(codon) == 3:
            adj_adjacent.append(codon)
            cursor += 1
            continue

        forword_window = adjust_look_forward(
            ref_codon_list[cursor: cursor + adjacent_window],
            seq_codon_list[cursor: cursor + adjacent_window])

        # Get codons in a window
        seq_list = []
        # depend on AA to align codon
        ref_list = []

        for i in range(forword_window):
            next_cursor = cursor + i
            if next_cursor >= len(seq_codon_list):
                break
            next_seq = seq_codon_list[next_cursor]
            next_ref = ref_codon_list[next_cursor]

            # Normal codon
            # if '-' not in next_seq and len(next_seq) == 3:
            #     break

            seq_list.append(next_seq)
            ref_list.append(next_ref)

        # print(ref_list)
        # print(seq_list)

        join_codon = ''.join(seq_list)
        codon_core = join_codon.replace('-', '')

        # Only codon can be translated will be adjusted.
        if len(codon_core) % 3 != 0:
            adj_adjacent.extend(seq_list)
            cursor += len(seq_list)
            continue

        # Adjust by original codon length
        codon_core_list = [
            codon_core[i:i + 3]
            for i in range(0, len(codon_core), 3)]

        new_codon_list = []
        for idx, c in enumerate(seq_list):
            if idx < len(codon_core_list):
                new_c = codon_core_list[idx]
                new_c += '-' * (len(c) - len(new_c))
                new_codon_list.append(new_c)
            else:
                new_codon_list.append('-' * len(c))

        adj_adjacent.extend(new_codon_list)
        cursor += len(new_codon_list)

    return adj_adjacent


def adjust_look_forward(ref_codons, seq_codons):
    # Match long deletion in ref and seq and remove them in this window
    ref_del = 0
    seq_del = 0
    for idx, (r, s) in enumerate(zip(ref_codons, seq_codons)):
        ref_del += r.count('-')
        seq_del += s.count('-')

        if ref_del and seq_del and (ref_del == seq_del):
            break

    return (idx + 1)


def adjust_codon_tail_del(ref_codon_list, seq_codon_list):
    new_ref = []
    new_seq = []
    for i, j in zip(ref_codon_list, seq_codon_list):
        if ('-' in i) and ('-' in j):
            if i.count('-') == j.count('-'):
                new_ref.append(i.replace('-', ''))
                new_seq.append(j.replace('-', ''))
                continue

        new_ref.append(i)
        new_seq.append(j)

    return new_ref, new_seq


def translate_aligned_codon(
        row, ref_aa, na_start, ref_codon, seq_codon):
    AA_list = []
    issue_list = []
    ins_list = []
    del_list = []
    stop_list = []

    for idx, (r, s) in enumerate(zip(ref_codon, seq_codon)):
        if (len(s) % 3) != 0:
            AA_list.append((idx + na_start, 'X'))
            issue_list.append((idx + na_start, r, s))
            continue
        try:
            trans_seq = Seq(s)
            aa = str(Seq(trans_seq).translate())

            if aa == '*':
                stop_list.append((idx + na_start, aa))
            elif aa == '-':
                del_list.append((idx + na_start, aa))
            elif len(aa) > 1:
                ins_list.append((idx + na_start, aa))

            AA_list.append((idx + na_start, aa))
        except TranslationError:
            AA_list.append((idx + na_start, 'X'))
            issue_list.append((idx + na_start, r, s))

    row['AA_codon_issue'] = ', '.join([
        f"{pos} ({r}, {s})"
        for pos, r, s in issue_list
    ])
    row['AA_num_codon_issue'] = len(issue_list)

    aa_seq = ''.join([aa for (pos, aa) in AA_list])
    row['AA_seq'] = aa_seq
    row['AA_length'] = len(aa_seq)

    row['AA_start'] = min([pos for (pos, aa) in AA_list])
    row['AA_stop'] = max([pos for (pos, aa) in AA_list])

    row['AA_ins_pos'] = ', '.join([
        str(pos)
        for (pos, aa) in ins_list
    ])
    row['AA_del_pos'] = ', '.join([
        str(pos)
        for (pos, aa) in del_list
    ])
    row['AA_stop_pos'] = ', '.join([
        str(pos)
        for (pos, aa) in stop_list
    ])

    row['AA_num_ins'] = len(ins_list)
    row['AA_num_del'] = len(del_list)
    row['AA_num_stop'] = len(stop_list)

    ref_aa = ref_aa[row['AA_start'] - 1: row['AA_stop']]
    mutations = [
        f"{a}{pos}{b}"
        for (a, (pos, b)) in zip(ref_aa, AA_list)
        if (a != b)
    ]

    row['Mutations'] = ', '.join(mutations)
    row['Num_mutations'] = len(mutations)
    # if row['AA_seq']:
    #     query = StripedSmithWaterman(ref_aa)
    #     alignment = query(row['AA_seq'])
    #     row['AA_start'] = alignment.query_begin + 1
    #     row['AA_stop'] = alignment.query_end + 1
    #     row['AA_num_ins'] = alignment.aligned_query_sequence.count('-')
    #     row['AA_num_del'] = alignment.aligned_target_sequence.count('-')
    #     row['AA_num_stop'] = alignment.aligned_target_sequence.count('*')

    return row


def get_pos_pre(str, char):
    pos = 0
    for c in str:
        if c == char:
            pos += 1
        else:
            break
    return pos


def get_pos_post(str, char):
    pos = len(str)
    for c in str[::-1]:
        if c == char:
            pos -= 1
        else:
            break
    return pos


# TODO auto detect AA alignment error
# consecutive positions mutations, or X
