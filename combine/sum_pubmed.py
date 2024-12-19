from .utils import count_number
from .utils import split_value_by_comma
from .translate_value import median_year
from .utils import int_sorter
from Utilities import create_binnned_year
import numpy as np


def summarize_pubmed(df, logger):
    summarize_report = []

    section = ["Summarize PubMed"]
    summarize_report.append(section)

    # logger.info('Human Likely, AI Likely:')
    # logger.info('Human Likely, AI Unlikely:')
    # logger.info('Reviwer2 Likely:')
    # logger.info('Reviwer2 Unlikely:')
    # logger.info('Human Likely, AI Unlikely:')

    # logger.info('Human Unlikely, AI Unlikely (not full):')

    section = ["Publish Year"]
    section.append(
        count_number([v for i, v in df.iterrows()], 'Year', sorter=int_sorter)
    )
    publish_year = [
        int(v['Year']) for i, v in df.iterrows()
        if v['Year'] and v['Year'] != 'NA']
    section.append(
        create_binnned_year(publish_year)
    )
    summarize_report.append(section)

    # section = ['Journals']
    # journals = count_number([v for i, v in df.iterrows()], 'Journal')
    # section.append(journals)
    # summarize_report.append(section)

    # section = ["NumSeq By Journal"]
    # num_seq = count_number([v for i, v in df.iterrows()], 'NumSeqs', sorter=int_sorter)
    # section.append(num_seq)
    # section.append('Total')
    # section.append(sum([int(v['NumSeqs']) for i, v in df.iterrows() if v['NumSeqs']]))
    # summarize_report.append(section)

    section = ["Host"]
    hosts = count_number([v for i, v in df.iterrows()], 'Host')
    section.append(hosts)
    summarize_report.append(section)

    section = ["Specimen"]
    specimen = count_number([v for i, v in df.iterrows()], 'Specimen')
    section.append(specimen)
    summarize_report.append(section)

    section = ["Median of Sample Year"]
    df['MedianYear'] = df['SampleYr'].apply(median_year)
    year = count_number([v for i, v in df.iterrows()], 'MedianYear', sorter=int_sorter)
    section.append(year)

    year = [int(v['MedianYear']) for i, v in df.iterrows() if v['MedianYear'] and v['MedianYear'] != 'NA']
    section.append(('Median IQR', np.percentile(year, 25), np.percentile(year, 50), np.percentile(year, 75)))
    section.append(create_binnned_year(year))
    summarize_report.append(section)

    section = ["Country"]
    # country_list = split_value_by_comma(df, 'Country')
    country_list = []
    for i, v in df.iterrows():
        countries = v['Country'].split(',')
        countries = [i.strip().capitalize() for i in countries]
        countries = sorted(list(set(countries)))
        country_list.append(', '.join(countries) if len(countries) < 4 else 'Multinational')
    country = count_number(country_list)
    section.append(country)
    summarize_report.append(section)

    # section = ["Country W/WO"]
    # country = count_number(
    #     [v for i, v in df.iterrows()], 'Country', translater=translate_country)
    # section.append(country)
    # summarize_report.append(section)

    section = ["Gene"]
    gene_list = split_value_by_comma(df, 'Gene')
    gene = count_number(gene_list)
    logger.info('Gene')
    section.append(gene)
    summarize_report.append(section)

    section = ["Seq method"]
    methods = count_number([v for i, v in df.iterrows()], 'SeqMethod')
    section.append(methods)
    summarize_report.append(section)

    section = ["End of Report"]
    summarize_report.append(section)


    for section in summarize_report:
        for pid, part in enumerate(section):
            if isinstance(part, tuple):
                logger.info(*part)
            elif isinstance(part, list):
                logger.info(*part)
            else:
                logger.info(part)
            if pid < len(section) - 1:
                logger.info('-' * 80)
        logger.info('=' * 80)
