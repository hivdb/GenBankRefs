from .utils import count_number
from .utils import int_sorter
from .sum_genbank import summarize_genbank_by_seq
from .utils import get_values_of_value_count_list
from Utilities import create_binnned_year
import re


def summarize_combined_data(combined, features, genes, logger):
    logger.info('Summarize PubMed GenBank Match')
    logger.info('=' * 80)

    matches = combined[(combined['match'] == 'Yes')]

    publish_year = count_number([
        v for i, v in matches.iterrows()], 'Year', sorter=int_sorter)
    logger.info('Publish Year')
    logger.info(publish_year)

    logger.info('*' * 80)
    publish_year = [
        int(v['Year']) for i, v in matches.iterrows()
        if v['Year'] and v['Year'] != 'NA']
    logger.info(create_binnned_year(publish_year))
    logger.info('=' * 80)

    # journals = count_number([v for i, v in matches.iterrows()], 'Journal')
    # logger.info('Journals')
    # logger.info(journals)
    # logger.info('=' * 80)

    methods = count_number(
        [v for i, v in matches.iterrows()], 'SeqMethod (PM)')
    logger.info('Seq method')
    logger.info(methods)
    logger.info('=' * 80)

    accessions = set([
        j.strip()
        for i, v in matches.iterrows()
        for j in v['GenBank (GB)'].split(',')
        ])
    num_seq = len(accessions)
    logger.info('After matching, NumSeq from GenBank:', num_seq)
    logger.info('=' * 80)

    features = features[features['Accession'].isin(list(accessions))]

    genes = genes[genes['Accession'].isin(list(accessions))]

    logger.info('Matched pubmed genbank:')
    summarize_genbank_by_seq(features, genes, logger)

    logger.info('Pubmed Supplement GenBank')
    for name in ['Hosts', 'Specimen', 'SampleYr', 'Countries', 'Genes', 'SeqMethod']:

        count = 0
        for idx, row in combined[combined['match'] == 'Yes'].iterrows():
            p_value = row[f"{name} (PM)"]
            g_value = row[f"{name} (GB)"]
            g_value = re.sub(r"\s\(\d+\)", "", g_value)

            p_value = [p for p in p_value.split(',') if p.strip() and p.strip().upper() != 'NA']
            g_value = [g for g in g_value.split(',') if g.strip() and g.strip().upper() != 'NA']

            if p_value and not g_value:
                count += 1

        logger.info(name, count)

    # logger.info('Similar virus')
    # logger.info(summarize_similarity(combined, 'Viruses'))

    # logger.info('Similar hosts')
    # logger.info(summarize_similarity(combined, 'Hosts'))

    # logger.info('Similar Specimens')
    # logger.info(summarize_similarity(combined, 'Specimen'))

    # logger.info('Similar countries')
    # logger.info(summarize_similarity(combined, 'Countries'))

    # logger.info('Similar Genes')
    # logger.info(summarize_similarity(combined, 'Genes'))


def summarize_similarity(df, col_name):

    count = 0
    for i, row in df.iterrows():
        pubmed = row[f'{col_name} (PM)']
        pubmed = [i.strip().lower() for i in pubmed.split(',') if i.strip().lower() != 'NA']

        genbank = row[f'{col_name} (GB)']
        genbank = get_values_of_value_count_list(genbank) if genbank else set()
        genbank = [i.lower() for i in genbank if i.lower() != 'NA']

        if (set(pubmed) & set(genbank)):
            count += 1

    return count
