from Bio import SeqIO
from Bio import Entrez
from Bio.Blast.Applications import NcbiblastpCommandline
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast.Applications import NcbiblastxCommandline
from Bio.Blast import NCBIXML
from multiprocessing import Pool
from pathlib import Path
from xml.parsers.expat import ExpatError

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


def perform_blast(acc, order, query_seq, db_name, func, blast_name):
    """
    Input: sample_seq (str): sequence to compare/blast
        db_name (str), prebuilt db by calling makeblastdb
        output_file (str): tmp file name to store results, cleared each iteration
    Output:
        a dictionary of keys including {e_value, percent_identity, alignment_length, overlap}
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

            if blast_name == 'blastn':
                NA_start = hsp.sbjct_start
                NA_stop = hsp.sbjct_end
                AA_start = ((NA_start - 1) // 3) + 1
                AA_stop = (NA_stop // 3)
            else:
                AA_start = hsp.sbjct_start
                AA_stop = hsp.sbjct_end
                NA_start = (AA_start - 1) * 3 + 1
                NA_stop = AA_stop * 3

            blast_result.append({
                'hit_name': alignment.hit_def,
                "e_value": e_value,
                "pcnt_id": percent_identity,
                "align_len": alignment_length * 3 if blast_name == 'blastp' else alignment_length,
                "overlap": overlap,
                'blast_name': blast_name,
                'AA_start': AA_start,
                'AA_stop': AA_stop,
                'NA_start': NA_start,
                'NA_stop': NA_stop,
            })

    Path(input_file).unlink(missing_ok=True)
    Path(output_file).unlink(missing_ok=True)

    return blast_result


def is_reference_genome(acc):
    prefix_list = ['NC', 'NG', 'NM', 'NR']
    for item in prefix_list:
        if acc.startswith(item):
            return True
    return False


def blast_sequence(idx, features, blast_aa_db_path, blast_na_db_path):
    """
        Try blastn, blastp, blastx for detecting the genes or segments of an isolate
        This function will decide the best blast result by alignment length of nucleic acid,
        if it's using blastp, the alignment length will be multiplied by 3.
    """

    blast_result = perform_blast(
        features['Accession'], features['Order'],
        features['NASeq'], blast_na_db_path,
        func=NcbiblastnCommandline, blast_name='blastn')

    if len(features['AASeq']) > 30:
        blast_result.extend(perform_blast(
            features['Accession'], features['Order'],
            features['AASeq'], blast_aa_db_path,
            func=NcbiblastpCommandline, blast_name='blastp'))

    blast_result.extend(perform_blast(
        features['Accession'], features['Order'],
        features['NASeq'], blast_aa_db_path,
        func=NcbiblastxCommandline, blast_name='blastx'))

    # if 'e_value' not in blast_data:
    #     print(blast_data, idx, len(features['AASeq']) <= 30)

    blast_result = sorted([
            b
            for b in blast_result
        ], key=lambda x: int(x['align_len']), reverse=True)

    if not blast_result:
        blast_data = {}
    else:
        blast_data = blast_result[0]

    features['hit_name'] = blast_data.get('hit_name', '')
    features['e_value'] = blast_data.get('e_value', '')
    features['pcnt_id'] = blast_data.get('pcnt_id', 0)
    features['align_len'] = blast_data.get('align_len', 0)
    features['blast_name'] = blast_data.get('blast_name', '')
    features['NA_start'] = blast_data.get('NA_start', '')
    features['NA_stop'] = blast_data.get('NA_stop', '')
    features['AA_start'] = blast_data.get('AA_start', '')
    features['AA_stop'] = blast_data.get('AA_stop', '')

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


def pooled_blast(features_list, virus_obj, poolsize=20):
    """
        For speeding up blast running, using multiprocessing method
        Input:
            features_list: a list of items, for each one contains NA or AA sequences
            db_name: blast database name
            poolsize: depending on how many cores the CPU has, you can adjust this number for running multiple blast program.
    """

    with Pool(poolsize) as pool:
        parameters = [
            (idx, f, virus_obj.BLAST_AA_DB_PATH, virus_obj.BLAST_NA_DB_PATH)
            for idx, f in enumerate(features_list)
        ]
        alignment_result = []
        for count, i in enumerate(pool.starmap(blast_sequence, parameters)):
            alignment_result.append(i)

    return alignment_result
