from .utils import count_number
from .utils import split_value_by_comma
from .translate_value import categorize_host_specimen
from .translate_value import median_year
from .utils import sum_value_count
from .utils import int_sorter
from .translate_value import translate_country
from Utilities import create_binnned_year
import numpy as np


def summarize_pubmed_data(df, logger):
    logger.info("PubMed")

    publish_year = count_number([v for i, v in df.iterrows()], 'Year', sorter=int_sorter)
    logger.info('Publish Year')
    logger.info(publish_year)
    publish_year = [
        int(v['Year']) for i, v in df.iterrows()
        if v['Year'] and v['Year'] != 'NA']
    logger.info(create_binnned_year(publish_year))
    logger.info('=' * 40)

    journals = count_number([v for i, v in df.iterrows()], 'Journal')
    logger.info('Journals')
    logger.info(journals)
    logger.info('=' * 40)

    num_seq = count_number([v for i, v in df.iterrows()], 'NumSeqs', sorter=int_sorter)
    logger.info('NumSeq By Journal')
    logger.info(num_seq)
    logger.info('Total', sum([int(v['NumSeqs']) for i, v in df.iterrows() if v['NumSeqs']]))
    logger.info('=' * 40)

    categorize_host_specimen(df, 'Host', 'IsolateType')

    hosts = count_number([v for i, v in df.iterrows()], 'CleanedHost')
    logger.info('Host')
    logger.info(hosts)
    logger.info('=' * 40)

    specimen = count_number([v for i, v in df.iterrows()], 'CleanedSpecimen')
    logger.info('Specimen')
    logger.info(specimen)
    logger.info('=' * 40)

    df['MedianYear'] = df['SampleYr'].apply(median_year)
    year = count_number([v for i, v in df.iterrows()], 'MedianYear', sorter=int_sorter)
    logger.info('Median of Sample Year')
    logger.info(year)
    year = [int(v['MedianYear']) for i, v in df.iterrows() if v['MedianYear'] and v['MedianYear'] != 'NA']
    logger.info('median', np.percentile(year, 25), np.percentile(year, 50), np.percentile(year, 75))
    logger.info(create_binnned_year(year))
    logger.info('=' * 40)

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
    logger.info('=' * 40)

    country = count_number(
        [v for i, v in df.iterrows()], 'Country', translater=translate_country)
    logger.info('Country W/WO')
    logger.info(country)
    logger.info('=' * 40)

    gene_list = split_value_by_comma(df, 'Gene')
    gene = count_number(gene_list)
    logger.info('Gene')
    logger.info(gene)
    logger.info('=' * 40)

    methods = count_number([v for i, v in df.iterrows()], 'SeqMethod')
    logger.info('Seq method')
    logger.info(methods)
    logger.info('=' * 40)

    logger.info('\n\n', '*' * 40, '\n\n')
