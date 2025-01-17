from Utilities import median_year
from Utilities import with_country
from Utilities import count_number
from Utilities import int_sorter
from Utilities import create_binned_pcnts
from Utilities import create_binned_seq_lens
from Utilities import create_binnned_year
from collections import defaultdict


def summarize_genbank(genbank_ref, genbank_feature, genbank_genes, virus_obj):
    genbank_ref = genbank_ref.fillna('')

    genbank_feature = genbank_feature.fillna('')

    genbank_genes = genbank_genes.fillna('')

    logger = virus_obj.get_logger('genbank')
    logger.report(summarize_genbank_by_ref(genbank_ref))

    logger.report(summarize_genbank_by_seq(
        genbank_feature, genbank_genes))

    logger.report(summarize_genbank_full_genome(
        genbank_ref, genbank_feature,
        full_gene_set=virus_obj.GENES))


def summarize_genbank_by_ref(df):
    summarize_report = []

    section = ['Summarize Genbank By Ref']
    summarize_report.append(section)

    df['MedianPublishYear'] = df['Year'].apply(median_year)
    publish_year = count_number([
        v for i, v in df.iterrows()], 'MedianPublishYear', sorter=int_sorter)
    section = ['Publish Year']
    section.append(publish_year)

    publish_year = [
        int(v['MedianPublishYear']) for i, v in df.iterrows()
        if v['MedianPublishYear'] and v['MedianPublishYear'] != 'NA']
    section.append(create_binnned_year(publish_year))

    summarize_report.append(section)

    # Journal information not included
    # journal_values = [row['Journal'].split(',')[0].strip()
    #                   for _, row in df.iterrows() if 'Journal' in row and pd.notnull(row['Journal'])]
    # cleaned_entries = [remove_parenthesis(entry) for entry in journal_values]
    # journals = count_number([{'Journal': value}
    #                         for value in cleaned_entries], 'Journal')
    # section = ['Journals']
    # section.append(journals)

    section = ['References Num Isolates (Sequences)']
    df['NumSeq (GB)'] = df['accession'].apply(lambda x: len(x.split(',')))
    num_seqs = count_number(
        [v for i, v in df.iterrows()], 'NumSeq (GB)', sorter=int_sorter)
    section.append(num_seqs)

    section.append((
        'Total', len(set([
            j.strip()
            for i, v in df.iterrows() if v['NumSeq (GB)']
            for j in v['accession'].split(',')
            if j.strip()
        ]))))
    summarize_report.append(section)

    section = ["End of Report"]
    summarize_report.append(section)

    return summarize_report


def summarize_genbank_by_seq(df, genes_df):
    summarize_report = []

    section = ['Summarize Genbank By Seq']
    summarize_report.append(section)

    section = ['Host']
    hosts = count_number([v for i, v in df.iterrows()], 'Host')
    section.append(hosts)
    summarize_report.append(section)

    section = ['Specimens']
    specimen = count_number([v for i, v in df.iterrows()], 'isolate_source')
    section.append(specimen)
    summarize_report.append(section)

    section = ['RecordYears']
    year = count_number(
        [v for i, v in df.iterrows()], 'RecordYear', sorter=int_sorter)
    section.append(year)

    year = [int(v['RecordYear']) for i, v in df.iterrows() if v['RecordYear']]
    section.append(create_binnned_year(year))
    summarize_report.append(section)

    section = ['Sample Years']
    year = count_number([
        int(v['IsolateYear'])
        if v['IsolateYear'] else '' for i, v in df.iterrows()],
        sorter=int_sorter)
    section.append(year)

    year = [int(v['IsolateYear']) for i, v in df.iterrows() if v['IsolateYear'] and v['IsolateYear'] != 'NA']
    section.append(create_binnned_year(year))
    summarize_report.append(section)

    section = ['Countries']
    country = count_number(
        [v for i, v in df.iterrows()], 'Country')
    section.append(country)
    summarize_report.append(section)

    section = ['Countries W/WO']
    country = count_number(
        [v for i, v in df.iterrows()], 'Country',
        translater=with_country)
    section.append(country)
    section.append('=' * 40)

    section = ['Genes']
    genes = count_number(
        [v for i, v in df.iterrows()], 'Genes')
    section.append(genes)
    summarize_report.append(section)

    section = ['AlignLens']
    aligns = [int(v['align_len']) for i, v in genes_df.iterrows()]
    section.append(create_binned_seq_lens(aligns))
    summarize_report.append(section)

    section = ['NA length']
    num_na = [int(v['NA_length']) for i, v in genes_df.iterrows()]
    section.append(create_binned_seq_lens(num_na))
    summarize_report.append(section)

    section = ['AA length']
    num_aa = [int(v['AA_length']) for i, v in genes_df.iterrows()]
    section.append(create_binned_seq_lens(num_aa))
    summarize_report.append(section)

    section = ['PcntIDs']
    pcnt_ident = [float(v['pcnt_id']) for i, v in genes_df.iterrows()]
    section.append(create_binned_pcnts(pcnt_ident))

    summarize_report.append(section)
    section = ['End of report']
    summarize_report.append(section)

    return summarize_report


def summarize_genbank_full_genome(
        ref_df, features_df, full_gene_set):

    summarize_report = []

    section = ['Summarize references with full genomes']
    summarize_report.append(section)

    num_ref = 0
    # num_seq = 0
    for i, row in ref_df.iterrows():
        accessions = set([
            a.strip()
            for a in row['accession'].split(',')
        ])

        features = features_df[
            features_df['Accession'].isin(accessions)]

        gene_count = defaultdict(int)
        for gene_list in features['Genes']:
            for gene in gene_list.split(','):
                gene = gene.strip()
                gene_count[gene] += 1

        genome = 0
        if 'genome' in gene_count:
            genome += gene_count['genome']
            del gene_count['genome']

        if set(gene_count.keys()) == set(full_gene_set):
            genome += min(gene_count.values())

        if genome:
            # num_seq += genome
            num_ref += 1

    section = ['Number of References with Full genome:']
    section.append(num_ref)
    # section.append('Number of full genome seq', num_seq)

    summarize_report.append(section)

    section = ['End of report']
    summarize_report.append(section)

    return summarize_report
