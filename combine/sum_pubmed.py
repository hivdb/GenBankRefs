from .utils import count_number
from .utils import split_value_by_comma
from .translate_value import median_year
from .utils import int_sorter
from Utilities import create_binnned_year
import numpy as np


def summarize_pubmed(df, logger):
    logger.info("Summarize PubMed")
    logger.info('=' * 80)

    publish_year = count_number([v for i, v in df.iterrows()], 'Year', sorter=int_sorter)
    logger.info('Publish Year')
    logger.info(publish_year)

    logger.info('*' * 80)

    publish_year = [
        int(v['Year']) for i, v in df.iterrows()
        if v['Year'] and v['Year'] != 'NA']
    logger.info(create_binnned_year(publish_year))
    logger.info('=' * 80)

    # journals = count_number([v for i, v in df.iterrows()], 'Journal')
    # logger.info('Journals')
    # logger.info(journals)
    # logger.info('=' * 80)

    # num_seq = count_number([v for i, v in df.iterrows()], 'NumSeqs', sorter=int_sorter)
    # logger.info('NumSeq By Journal')
    # logger.info(num_seq)
    # logger.info('Total', sum([int(v['NumSeqs']) for i, v in df.iterrows() if v['NumSeqs']]))
    # logger.info('=' * 80)

    hosts = count_number([v for i, v in df.iterrows()], 'Host')
    logger.info('Host')
    logger.info(hosts)
    logger.info('=' * 80)

    specimen = count_number([v for i, v in df.iterrows()], 'Specimen')
    logger.info('Specimen')
    logger.info(specimen)
    logger.info('=' * 80)

    df['MedianYear'] = df['SampleYr'].apply(median_year)
    year = count_number([v for i, v in df.iterrows()], 'MedianYear', sorter=int_sorter)
    logger.info('Median of Sample Year')
    logger.info(year)
    year = [int(v['MedianYear']) for i, v in df.iterrows() if v['MedianYear'] and v['MedianYear'] != 'NA']

    logger.info('*' * 80)
    logger.info('Median IQR', np.percentile(year, 25), np.percentile(year, 50), np.percentile(year, 75))

    logger.info('*' * 80)
    logger.info(create_binnned_year(year))
    logger.info('=' * 80)

    # country_list = split_value_by_comma(df, 'Country')
    country_list = []
    for i, v in df.iterrows():
        countries = v['Country'].split(',')
        countries = [i.strip().capitalize() for i in countries]
        countries = sorted(list(set(countries)))
        country_list.append(', '.join(countries) if len(countries) < 4 else 'Multinational')
    country = count_number(country_list)
    logger.info('Country')
    logger.info(country)
    logger.info('=' * 80)

    # country = count_number(
    #     [v for i, v in df.iterrows()], 'Country', translater=translate_country)
    # logger.info('Country W/WO')
    # logger.info(country)
    # logger.info('=' * 80)

    gene_list = split_value_by_comma(df, 'Gene')
    gene = count_number(gene_list)
    logger.info('Gene')
    logger.info(gene)
    logger.info('=' * 80)

    methods = count_number([v for i, v in df.iterrows()], 'SeqMethod')
    logger.info('Seq method')
    logger.info(methods)
    logger.info('=' * 80)
    logger.info('# End of Section')
    logger.info('=' * 80)
