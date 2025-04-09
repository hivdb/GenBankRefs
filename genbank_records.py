from Bio import SeqIO
import pandas as pd
from operator import itemgetter

from GenBankFunctions import filter_by_taxonomy
from GenBankFunctions import pooled_blast_genes
from GenBankFunctions import detect_additional_genes
from GenBankFunctions import detect_gene_by_biopython

from Utilities import extract_year_from_date_fields
from Utilities import extract_year_from_journal
from Utilities import process_author_field


# This is called by the main function to allow users to determine whether
# BLAST is run
def select_run_blast(message='Run blast?', show_option=True):

    option = input(f'{message} [y/n]: ')
    option = option.lower()

    option = option == 'y'
    if show_option:
        if option:
            print('=> Choose Yes')
        else:
            print('=> Choose No')
    return option


def parse_genbank_records(genbank_file):
    reference_list = []
    feature_list = []
    gene_list = []

    nonvirus_list = []
    nonclinical_list = []
    lab_host_list = []

    total_ref_list = []

    exclusion_keywords = [
        'patent',
        'FDA',
        'Modified Microbial Nucleic Acid',
        'CONSTRUCT',
        'COMPOSITIONS',
        'monoclonal antibody',
        'MICROARRAY',
        'conformation',
        # 'Chain',
    ]

    seq_methods = parse_seq_method(genbank_file)

    with open(genbank_file, "r") as handle:
        for record in SeqIO.parse(handle, "genbank"):
            # total_record += 1
            # if total_record % 10000 == 0:
            #     print(total_record)

            taxonomy = record.annotations['taxonomy']
            if len(taxonomy) == 0 or taxonomy[0] != 'Viruses':
                excluded_seq_data = filter_by_taxonomy(record)
                nonvirus_list.append(excluded_seq_data)
                continue

            # Check if the description contains any exclusion keyword
            description = record.description.lower()
            should_exclude = any(
                keyword.lower() in description
                for keyword in exclusion_keywords
            )

            refs, features, genes = process_one_record(record)

            acc = record.id.split('.')[0]
            features['SeqMethod'] = seq_methods.get(acc, '')

            total_ref_list.extend(refs)

            # if 'MAG' in record.annotations['keywords']:
            #     should_exclude = True
            if not should_exclude:
                if all(
                    'patent' in r['Journal'].lower()
                    for r in refs
                ):
                    should_exclude = True

            if should_exclude:
                nonclinical_list.append(features)
                continue

            if features['lab_host']:
                lab_host_list.append(features)
                nonclinical_list.append(features)
                continue

            if 'Biochemical characterization'.lower() in refs[0]['Title'].lower():
                lab_host_list.append(features)
                nonclinical_list.append(features)
                continue

            if 'mouse-adapted' in features['note']:
                lab_host_list.append(features)
                nonclinical_list.append(features)
                continue

            reference_list.extend(refs)
            feature_list.append(features)
            gene_list.extend(genes)

    print('# Total accessions:', len(feature_list) + len(nonvirus_list) + len(nonclinical_list))
    print('# Lab host accessions:', len(lab_host_list))

    print("# excluded non virus records:", len(nonvirus_list))
    print('# non clinical records (including lab host):', len(nonclinical_list))

    exclude_list = nonvirus_list + nonclinical_list

    return (
        total_ref_list,
        reference_list,
        feature_list,
        gene_list,
        pd.DataFrame(exclude_list)
        )


def process_one_record(record):

    accession = record.id.split('.', -1)[0]

    refs = extract_references(record.annotations, accession)

    features = {}
    features['Accession'] = accession

    features['Description'] = record.description
    features['record_date'] = record.annotations['date']
    features['organism'] = record.annotations['organism']

    feature_data = extract_features(record.features, accession)
    features['segment_source'] = feature_data.get('segment_source', '')
    features['Host'] = feature_data.get('host_source', '')
    features['isolate_source'] = feature_data.get(
        'isolation_source_source', '')
    features['IsolateName'] = feature_data.get('isolate_source') or feature_data.get('strain_source', '')
    features['country_region'] = feature_data.get(
        'geo_loc_name_source', '')
    features['collection_date'] = feature_data.get(
        'collection_date_source', '')

    features['lab_host'] = feature_data.get('lab_host_source', '')

    features['Seq'] = str(record.seq)
    features['SeqLength'] = len(record.seq)

    features['note'] = feature_data.get('note_source', '')

    gene_seq = [
        i
        for i in record.features
        if i.type == 'CDS'
    ]

    cds_names = []
    genes = []
    for idx, aa in enumerate(gene_seq):
        gene_name = None
        if 'gene' in aa.qualifiers:
            gene_name = aa.qualifiers['gene'][0].upper()
        elif 'product' in aa.qualifiers:
            gene_name = aa.qualifiers[
                'product'][0].upper().replace(' PROTEIN', '').strip()
        elif 'note' in aa.qualifiers:
            gene_name = aa.qualifiers[
                'note'][0].upper().replace(' PROTEIN', '').strip()

        # TODO, this many need check translation exists.
        cds_names.append(gene_name)

        if 'translation' in aa.qualifiers:
            aa_seq = str(aa.qualifiers['translation'][0])
            na_seq = str(aa.location.extract(record.seq))
            genes.append({
                'Accession': accession,
                'Gene': '',
                'CDS_NAME': gene_name,
                'Order': idx + 1,

                'AA_raw_seq': aa_seq,
                'AA_raw_length': len(aa_seq) if aa_seq else '',

                'NA_raw_seq': na_seq,
                'NA_raw_length': len(na_seq) if na_seq else '',

            })

    aa_seq = ''
    na_seq = str(record.seq)

    genes.append({
        'Accession': accession,
        'Gene': '',
        'CDS_NAME': 'isolate',
        'Order': 0,

        'AA_raw_seq': aa_seq,
        'AA_raw_length': len(aa_seq),

        'NA_raw_seq': na_seq,
        'NA_raw_length': len(na_seq),
    })

    features['cds'] = ', '.join(cds_names)
    features['NumSubSeqs'] = len(cds_names)

    return refs, features, genes


def extract_references(annotations, accession):
    ref_list = annotations.get("references")
    ref_data = []
    for ref in ref_list:
        ref_items = {}
        ref_items['accession'] = accession
        ref_items['Authors'] = ref.authors
        ref_items['Title'] = ref.title
        ref_items['Journal'] = ref.journal
        ref_items['PMID'] = ref.pubmed_id
        ref_data.append(ref_items)
    return ref_data


def extract_features(features, accession):
    feature_items = {}
    feature_items["accession"] = accession
    for feature in features:
        for key, value in feature.qualifiers.items():
            key = key + '_' + feature.type
            feature_items[key] = value[0]
            # print(f"Key {key}: value {value}\n")
    return feature_items


def process_references(references):
    """
    Processes reference data from GenBank

    Steps:
    1. Converts reference list into a DataFrame.
    2. Extracts the publication year from the 'Journal' field.
    3. Cleans and standardizes the 'Journal' field by removing submission dates and simplifying patent references.
    4. Processes the 'Authors' field for standard formatting.

    Returns:
    DataFrame: Processed reference dataset with cleaned metadata.
    """
    references = pd.DataFrame(references)
    references['Year'] = references['Journal'].apply(
        extract_year_from_journal)

    references['Year'] = pd.to_numeric(
        references['Year'], errors='coerce')
    references['Year'] = references['Year'].apply(
        lambda x: '' if pd.isna(x) else int(x))

    # Remove Submitted (date) from journal
    references['Journal'] = references['Journal'].str.replace(
        r"Submitted \(\d{2}-[A-Z]{3}-\d{4}\)", "", regex=True)

    references['Journal'] = references['Journal'].str.replace(
        r"(Patent).*", r"\1", regex=True)
    references['Authors'] = references['Authors'].apply(
        process_author_field)

    return references


def process_features(feature_list, genes, virus_obj):
    """
    Processes genomic feature data for a given virus.

    Steps:
    1. Converts feature list into a DataFrame.
    2. Cleans and processes the feature data using the virus-specific module.
    3. Extracts isolate and record years from date fields.
    4. Saves and reloads the processed feature data.
    5. Filters out isolates that lack detected genes.

    Returns:
    DataFrame: Cleaned and filtered genomic features.
    """
    features_df = pd.DataFrame(feature_list)
    # This uses data in the imported virus module to clean data in the feature table
    features_df = virus_obj.process_features(features_df, genes)

    features_df['RecordYear'] = features_df['record_date'].apply(
        extract_year_from_date_fields)
    features_df['IsolateYear'] = features_df['collection_date'].apply(
        extract_year_from_date_fields)

    features_df.to_excel(
        str(virus_obj.genbank_feature_check_file), index=False)
    features_df = pd.read_excel(
        str(virus_obj.genbank_feature_check_file)).fillna('')

    excluded_features = features_df[(features_df["Genes"].isna() | (features_df["Genes"] == ""))]
    print('# Excluded isolates without Gene', len(excluded_features))

    # Drop sequences with no detected genes
    features_df = features_df[~(features_df["Genes"].isna() | (features_df["Genes"] == ""))]
    print('# Accessions included', len(features_df))

    return features_df, excluded_features


def process_gene_list(gene_list, run_blast, virus_obj):
    """
    Processes the list of genes associated with a virus, optionally running BLAST to detect additional genes.

    Steps:
    1. If BLAST is enabled (run_blast == True), build a BLAST database and detect additional genes.
    2. Identify genes that were missed initially and attempt detection using Biopython.
    3. Sort and store the gene data into an Excel file.
    4. If BLAST isn't enabled, check if a previously processed gene excel file exists and load it.
    5. If no pre-existing file is found, create a new gene dataset and save it.
    6. If BLAST isn't enabled & no gene excel file, load directly from gene_list (check)

    Returns:
    DataFrame: Processed gene dataset.
    """
    if run_blast:
        virus_obj.build_blast_db()

        gene_list2 = pooled_blast_genes(gene_list, virus_obj)
        additional_gene_list = detect_additional_genes(
            gene_list, gene_list2, virus_obj)

        with_gene_acc = [
            i['Accession']
            for i in gene_list2 + additional_gene_list
            if i['Gene']
        ]
        without_gene = [
            i
            for i in gene_list
            if i['Accession'] not in with_gene_acc
        ]
        print("Accessions no gene after blast", len(without_gene))
        print('CDS_NAME', set([i['CDS_NAME'] for i in without_gene]))
        missing_genes = detect_gene_by_biopython(without_gene, virus_obj)
        print('Detected missing genes', len(missing_genes))
        print('Genes', [i['Gene'] for i in missing_genes])
        gene_list = gene_list2 + additional_gene_list + missing_genes
        gene_list.sort(key=itemgetter('Accession', 'Gene'))

        gene_df = pd.DataFrame(gene_list)
        gene_df['SeqID'] = gene_df.index + 1
        gene_df.to_excel(str(virus_obj.genbank_gene_file), index=False)
        gene_df = pd.read_excel(
            str(virus_obj.genbank_gene_file)).fillna('')

    # genbank_gene_file is created either using blast or without blast
    # depending on the user input
    # This statement loads the existing genbank feature file which may or may not
    # include blast data. This makes it possible to avoid re-running BLAST
    elif virus_obj.genbank_gene_file.exists():
        gene_df = pd.read_excel(
            str(virus_obj.genbank_gene_file)).fillna('')

    # This statement creates gene df without running BLAST
    else:
        gene_df = pd.DataFrame(gene_list)
        gene_df.to_excel(str(virus_obj.genbank_gene_file), index=False)
        gene_df = pd.read_excel(
            str(virus_obj.genbank_gene_file)).fillna('')
    gene_df = virus_obj.process_gene_list(gene_df)
    gene_df.to_excel(str(virus_obj.genbank_gene_file), index=False)

    return gene_df


def parse_seq_method(genbank_file):

    matches = {}
    acc = None
    with open(genbank_file, "r") as fd:
        for line in fd.readlines():
            if line.startswith('ACCESSION'):
                acc = line.strip().split()[-1]
            elif 'Sequencing Technology ::' in line:
                matches[acc] = line.split('::')[-1].strip()

    matches = {
        k: translate_seq_method(v)
        for k, v in matches.items()
    }

    return matches


def translate_seq_method(seq_method):
    seq_method_list = set()
    if 'Sanger dideoxy sequencing' in seq_method:
        seq_method_list.add('Sanger')

    NGS_methods = [
        'Illumina',
        'ONT',
        'Oxford Nanopore',
        'IonTorrent',
        'Oxford Nanopore Technologies',
        '454',
        'MGISEQ-2000',
        'MinIon',
        'Illumina MiSeq and MiSeq FGx',
        'NexSeq500',
        'Illumina HiSeq 2500',
        'Illumina MiSeq',
        'Oxford Nanopore MinION',
        'Illumina NovaSeq, HiSeq 2500 and MiSeq',
        'Illumina HiSeq1500',
        'HotShot'
    ]
    for i in NGS_methods:
        if i in seq_method:
            seq_method_list.add('NGS')

    if len(seq_method_list) == 2:
        return 'Sanger and NGS'
    else:
        return ''.join(list(seq_method_list))
