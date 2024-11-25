from .translate_value import median_year
from .translate_value import translate_country
from .translate_value import translate_gene
from .translate_value import translate_specimen
from .translate_value import translate_hosts
from .utils import count_number
from .utils import merge_genbank_list_columns


def summarize_genbank_by_seq(df):
    df['MedianPublishYear'] = df['year'].apply(median_year)
    publish_year = count_number([v for i, v in df.iterrows()], 'MedianPublishYear')
    print('Publish Year')
    print(publish_year)
    print('=' * 40)

    # journal_values = [row['journal'].split(',')[0].strip()
    #                   for _, row in df.iterrows() if 'journal' in row and pd.notnull(row['journal'])]
    # cleaned_entries = [remove_parenthesis(entry) for entry in journal_values]
    # journals = count_number([{'journal': value}
    #                         for value in cleaned_entries], 'journal')
    # print('Journals')
    # print(journals)
    # print('=' * 40)

    df['NumSeq (GB)'] = df['accession'].apply(lambda x: len(x.split(',')))
    num_seqs = count_number([v for i, v in df.iterrows()], 'NumSeq (GB)')
    print('NumSeq')
    print(num_seqs)
    print('=' * 40)

    host = merge_genbank_list_columns(
        [v for i, v in df.iterrows()], 'Hosts',
        translater=translate_hosts)
    print('Host')
    print(host)
    print('=' * 40)

    specimen = merge_genbank_list_columns(
        [v for i, v in df.iterrows()], 'Specimens',
        translater=translate_specimen)
    print('Specimens')
    print(specimen)
    print('=' * 40)

    year = merge_genbank_list_columns(
        [v for i, v in df.iterrows()], 'RecordYears')
    print('RecordYears')
    print(year)
    print('=' * 40)

    year = merge_genbank_list_columns(
        [v for i, v in df.iterrows()], 'IsolateYears')
    print('Sample Years')
    print(year)
    print('=' * 40)

    country = year = merge_genbank_list_columns(
        [v for i, v in df.iterrows()], 'IsolateYears',
        translater=translate_country)
    print('Countries')
    print(country)
    print('=' * 40)

    genes = year = merge_genbank_list_columns(
        [v for i, v in df.iterrows()], 'Gene', translater=translate_gene)
    print('Gene')
    print(genes)
    print('=' * 40)

    # Works?
    # df['AverageAlignLength'] = df['AlignLens'].apply(calculate_average_length)
    # host = count_number([v for i, v in df.iterrows()], 'AverageAlignLength')
    # print('AlignLens')
    # print(host)
    # print('=' * 40)

    aligns = merge_genbank_list_columns(
        [v for i, v in df.iterrows()], 'AlignLens')
    print('AlignLens')
    print(aligns)
    print('=' * 40)

    # df['MostFrequentRange'] = df['PcntIDs'].apply(most_frequent_range)
    # host = count_number([v for i, v in df.iterrows()], 'MostFrequentRange')
    # print('MostFrequentRange')
    # print(host)
    # print('=' * 40)

    pcnt_ident = merge_genbank_list_columns(
        [v for i, v in df.iterrows()], 'PcntIDs')
    print('PcntIDs')
    print(pcnt_ident)
    print('=' * 40)

    print('\n\n', '*' * 40, '\n\n')
