from .translate_value import median_year
from .translate_value import translate_country
from .translate_value import translate_gene
from .translate_value import translate_specimen
from .translate_value import translate_hosts
from .utils import count_number
from .utils import merge_genbank_list_columns
from .utils import sum_value_count
from .utils import int_sorter
from .translate_value import categorize_host_specimen
from Utilities import extract_year_from_date_fields
from Utilities import create_binned_pcnts
from Utilities import create_binned_seq_lens
from Utilities import create_binnned_year


def summarize_genbank_by_ref(df):
    print('Summarize Genbank By Ref')

    df['MedianPublishYear'] = df['year'].apply(median_year)
    publish_year = count_number([v for i, v in df.iterrows()], 'MedianPublishYear', sorter=int_sorter)
    print('Publish Year')
    print(publish_year)
    publish_year = [int(v['MedianPublishYear']) for i, v in df.iterrows() if v['MedianPublishYear'] and v['MedianPublishYear'] != 'NA']
    print(create_binnned_year(publish_year))
    print('=' * 40)

    # Journal information not included
    # journal_values = [row['journal'].split(',')[0].strip()
    #                   for _, row in df.iterrows() if 'journal' in row and pd.notnull(row['journal'])]
    # cleaned_entries = [remove_parenthesis(entry) for entry in journal_values]
    # journals = count_number([{'journal': value}
    #                         for value in cleaned_entries], 'journal')
    # print('Journals')
    # print(journals)
    # print('=' * 40)

    df['NumSeq (GB)'] = df['accession'].apply(lambda x: len(x.split(',')))
    num_seqs = count_number(
        [v for i, v in df.iterrows()], 'NumSeq (GB)', sorter=int_sorter)
    print('NumSeq')
    print(num_seqs)
    print('=' * 40)


def summarize_genbank_by_seq(df):
    print('Summarize Genbank By Seq')

    categorize_host_specimen(df, 'host', 'isolate_source')
    hosts = count_number([v for i, v in df.iterrows()], 'CleanedHost')
    print('Host')
    print(hosts)
    print('=' * 40)

    specimen = count_number([v for i, v in df.iterrows()], 'CleanedSpecimen')
    print('Specimens')
    print(specimen)
    print('=' * 40)

    df['record_year'] = df['record_date'].apply(extract_year_from_date_fields)
    year = count_number(
        [v for i, v in df.iterrows()], 'record_year', sorter=int_sorter)
    print('RecordYears')
    print(year)
    year = [int(v['record_year']) for i, v in df.iterrows() if v['record_year']]
    print(create_binnned_year(year))
    print('=' * 40)

    df['isolate_year'] = df['collection_date'].apply(
        extract_year_from_date_fields)
    year = count_number(
        [v for i, v in df.iterrows()], 'isolate_year', sorter=int_sorter)
    print('Sample Years')
    print(year)
    year = [int(v['isolate_year']) for i, v in df.iterrows() if v['isolate_year'] and v['isolate_year'] != 'NA']
    print(create_binnned_year(year))
    print('=' * 40)

    country = count_number(
        [v for i, v in df.iterrows()], 'country_region')
    print('Countries')
    print(country)
    print('=' * 40)

    country = count_number(
        [v for i, v in df.iterrows()], 'country_region',
        translater=translate_country)
    print('Countries W/WO')
    print(country)
    print('=' * 40)

    genes = count_number(
        [v for i, v in df.iterrows()], 'segment_source',
        translater=translate_gene)
    print('Genes')
    print(genes)
    print('=' * 40)

    aligns = [int(v['align_len']) for i, v in df.iterrows()]
    print('AlignLens')
    print(create_binned_seq_lens(aligns))
    print('=' * 40)

    pcnt_ident = [float(v['pcnt_id']) for i, v in df.iterrows()]
    print('PcntIDs')
    print(create_binned_pcnts(pcnt_ident))
    print('=' * 40)

    print('\n\n', '*' * 40, '\n\n')
