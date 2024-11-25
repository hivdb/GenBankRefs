from .utils import count_number
from .utils import split_value_by_comma
from .translate_value import categorize_host_specimen
from .translate_value import median_year
from .utils import sum_value_count
from .utils import int_sorter
from .translate_value import translate_gene
from .translate_value import translate_country
from Utilities import create_binnned_year


def summarize_pubmed_data(df):
    print("PubMed")

    publish_year = count_number([v for i, v in df.iterrows()], 'Year', sorter=int_sorter)
    print('Publish Year')
    print(publish_year)
    # sum_value_count(publish_year)
    print('=' * 40)

    journals = count_number([v for i, v in df.iterrows()], 'Journal')
    print('Journals')
    print(journals)
    print('=' * 40)

    num_seq = count_number([v for i, v in df.iterrows()], 'NumSeqs', sorter=int_sorter)
    print('NumSeq')
    print(num_seq)
    print('=' * 40)

    categorize_host_specimen(df, 'Host', 'IsolateType')

    hosts = count_number([v for i, v in df.iterrows()], 'CleanedHost')
    print('Host')
    print(hosts)
    print('=' * 40)

    specimen = count_number([v for i, v in df.iterrows()], 'CleanedSpecimen')
    print('Specimen')
    print(specimen)
    print('=' * 40)

    df['MedianYear'] = df['SampleYr'].apply(median_year)
    year = count_number([v for i, v in df.iterrows()], 'MedianYear', sorter=int_sorter)
    print('Median of Sample Year')
    print(year)
    year = [int(v['MedianYear']) for i, v in df.iterrows() if v['MedianYear'] and v['MedianYear'] != 'NA']
    print(create_binnned_year(year))
    print('=' * 40)

    country_list = split_value_by_comma(df, 'Country')
    country = count_number(country_list)
    print('Country')
    print(country)
    print('=' * 40)

    country = count_number(
        [v for i, v in df.iterrows()], 'Country', translater=translate_country)
    print('Country W/WO')
    print(country)
    print('=' * 40)

    gene_list = split_value_by_comma(df, 'Gene')
    gene = count_number(gene_list, translater=translate_gene)
    print('Gene')
    print(gene)
    print('=' * 40)

    methods = count_number([v for i, v in df.iterrows()], 'SeqMethod')
    print('Seq method')
    print(methods)
    print('=' * 40)

    print('\n\n', '*' * 40, '\n\n')
