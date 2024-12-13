from .translate_value import median_year
from .translate_value import translate_country
from .utils import count_number
from .utils import int_sorter
from Utilities import create_binned_pcnts
from Utilities import create_binned_seq_lens
from Utilities import create_binnned_year
from collections import defaultdict


def summarize_genbank_by_ref(df, logger):
    logger.info('Summarize Genbank By Ref')

    df['MedianPublishYear'] = df['Year'].apply(median_year)
    publish_year = count_number([
        v for i, v in df.iterrows()], 'MedianPublishYear', sorter=int_sorter)
    logger.info('Publish Year')
    logger.info(publish_year)
    publish_year = [
        int(v['MedianPublishYear']) for i, v in df.iterrows()
        if v['MedianPublishYear'] and v['MedianPublishYear'] != 'NA']
    logger.info(create_binnned_year(publish_year))
    logger.info('=' * 40)

    # Journal information not included
    # journal_values = [row['Journal'].split(',')[0].strip()
    #                   for _, row in df.iterrows() if 'Journal' in row and pd.notnull(row['Journal'])]
    # cleaned_entries = [remove_parenthesis(entry) for entry in journal_values]
    # journals = count_number([{'Journal': value}
    #                         for value in cleaned_entries], 'Journal')
    # logger.info('Journals')
    # logger.info(journals)
    # logger.info('=' * 40)

    df['NumSeq (GB)'] = df['accession'].apply(lambda x: len(x.split(',')))
    num_seqs = count_number(
        [v for i, v in df.iterrows()], 'NumSeq (GB)', sorter=int_sorter)
    logger.info('NumSeq')
    logger.info(num_seqs)

    logger.info('Total', len(set([
        j.strip()
        for i, v in df.iterrows() if v['NumSeq (GB)']
        for j in v['accession'].split(',')
        if j.strip()
        ])))
    logger.info('=' * 40)


def summarize_genbank_full_genome(
        ref_df, features_df, logger, full_gene_set, col_name='Gene'):

    num_ref = 0
    num_seq = 0
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

        if set(gene_count.keys()) == full_gene_set:
            genome += min(gene_count.values())

        if genome:
            num_seq += genome
            num_ref += 1

    logger.info('Full genome Ref')
    logger.info(num_ref)
    logger.info('Full genome seq')
    logger.info(num_seq)


def summarize_genbank_by_seq(df, genes_df, logger):
    logger.info('Summarize Genbank By Seq')

    hosts = count_number([v for i, v in df.iterrows()], 'Host')
    logger.info('Host')
    logger.info(hosts)
    logger.info('=' * 40)

    specimen = count_number([v for i, v in df.iterrows()], 'isolate_source')
    logger.info('Specimens')
    logger.info(specimen)
    logger.info('=' * 40)

    year = count_number(
        [v for i, v in df.iterrows()], 'RecordYear', sorter=int_sorter)
    logger.info('RecordYears')
    logger.info(year)
    year = [int(v['RecordYear']) for i, v in df.iterrows() if v['RecordYear']]
    logger.info(create_binnned_year(year))
    logger.info('=' * 40)

    year = count_number(
        [v for i, v in df.iterrows()], 'IsolateYear', sorter=int_sorter)
    logger.info('Sample Years')
    logger.info(year)
    year = [int(v['IsolateYear']) for i, v in df.iterrows() if v['IsolateYear'] and v['IsolateYear'] != 'NA']
    logger.info(create_binnned_year(year))
    logger.info('=' * 40)

    country = count_number(
        [v for i, v in df.iterrows()], 'Country')
    logger.info('Countries')
    logger.info(country)
    logger.info('=' * 40)

    country = count_number(
        [v for i, v in df.iterrows()], 'Country',
        translater=translate_country)
    logger.info('Countries W/WO')
    logger.info(country)
    logger.info('=' * 40)

    genes = count_number(
        [v for i, v in df.iterrows()], 'Genes')
    logger.info('Genes')
    logger.info(genes)
    logger.info('=' * 40)

    aligns = [int(v['align_len']) for i, v in genes_df.iterrows()]
    logger.info('AlignLens')
    logger.info(create_binned_seq_lens(aligns))
    logger.info('=' * 40)

    num_na = [int(v['NumNA']) for i, v in genes_df.iterrows()]
    logger.info('NA length')
    logger.info(create_binned_seq_lens(num_na))
    logger.info('=' * 40)

    num_aa = [int(v['NumAA']) for i, v in genes_df.iterrows()]
    logger.info('AA length')
    logger.info(create_binned_seq_lens(num_aa))
    logger.info('=' * 40)

    pcnt_ident = [float(v['pcnt_id']) for i, v in genes_df.iterrows()]
    logger.info('PcntIDs')
    logger.info(create_binned_pcnts(pcnt_ident))
    logger.info('=' * 40)

    logger.info('\n\n', '*' * 40, '\n\n')
