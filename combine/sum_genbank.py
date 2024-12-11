from .translate_value import median_year
from .translate_value import translate_country
from .translate_value import translate_gene
from .utils import count_number
from .utils import int_sorter
from .utils import split_value_count
from Utilities import create_binned_pcnts
from Utilities import create_binned_seq_lens
from Utilities import create_binnned_year


def summarize_genbank_by_ref(df):
    print('Summarize Genbank By Ref')

    df['MedianPublishYear'] = df['Year'].apply(median_year)
    publish_year = count_number([
        v for i, v in df.iterrows()], 'MedianPublishYear', sorter=int_sorter)
    print('Publish Year')
    print(publish_year)
    publish_year = [
        int(v['MedianPublishYear']) for i, v in df.iterrows()
        if v['MedianPublishYear'] and v['MedianPublishYear'] != 'NA']
    print(create_binnned_year(publish_year))
    print('=' * 40)

    # Journal information not included
    # journal_values = [row['Journal'].split(',')[0].strip()
    #                   for _, row in df.iterrows() if 'Journal' in row and pd.notnull(row['Journal'])]
    # cleaned_entries = [remove_parenthesis(entry) for entry in journal_values]
    # journals = count_number([{'Journal': value}
    #                         for value in cleaned_entries], 'Journal')
    # print('Journals')
    # print(journals)
    # print('=' * 40)

    df['NumSeq (GB)'] = df['accession'].apply(lambda x: len(x.split(',')))
    num_seqs = count_number(
        [v for i, v in df.iterrows()], 'NumSeq (GB)', sorter=int_sorter)
    print('NumSeq')
    print(num_seqs)

    print('Total', len(set([
        j.strip()
        for i, v in df.iterrows() if v['NumSeq (GB)']
        for j in v['accession'].split(',')
        if j.strip()
        ])))
    print('=' * 40)


def summarize_genbank_full_genome(
        df, col_name='Gene', full_gene_set={'L', 'S', 'M'}):

    potential = []
    total = 0
    for i, row in df.iterrows():
        count_list = []
        value_list = []
        for i in row[col_name].split(','):
            value, count = split_value_count(i)
            count_list.append(int(count))
            value_list.extend([value] * int(count))

        if set(value_list) == full_gene_set and len(set(count_list)) == 1:
            potential.append(row)
            total += count_list[0]

    print('Full genome Ref')
    print(len(potential))
    print('Full genome seq')
    print(total)


def summarize_genbank_by_seq(df):
    print('Summarize Genbank By Seq')

    hosts = count_number([v for i, v in df.iterrows()], 'Host')
    print('Host')
    print(hosts)
    print('=' * 40)

    specimen = count_number([v for i, v in df.iterrows()], 'isolate_source')
    print('Specimens')
    print(specimen)
    print('=' * 40)

    year = count_number(
        [v for i, v in df.iterrows()], 'RecordYear', sorter=int_sorter)
    print('RecordYears')
    print(year)
    year = [int(v['RecordYear']) for i, v in df.iterrows() if v['RecordYear']]
    print(create_binnned_year(year))
    print('=' * 40)

    year = count_number(
        [v for i, v in df.iterrows()], 'IsolateYear', sorter=int_sorter)
    print('Sample Years')
    print(year)
    year = [int(v['IsolateYear']) for i, v in df.iterrows() if v['IsolateYear'] and v['IsolateYear'] != 'NA']
    print(create_binnned_year(year))
    print('=' * 40)

    country = count_number(
        [v for i, v in df.iterrows()], 'Country')
    print('Countries')
    print(country)
    print('=' * 40)

    country = count_number(
        [v for i, v in df.iterrows()], 'Country',
        translater=translate_country)
    print('Countries W/WO')
    print(country)
    print('=' * 40)

    genes = count_number(
        [v for i, v in df.iterrows()], 'Genes',
        translater=translate_gene)
    print('Genes')
    print(genes)
    print('=' * 40)

    aligns = [int(v['align_len']) for i, v in df.iterrows()]
    print('AlignLens')
    print(create_binned_seq_lens(aligns))
    print('=' * 40)

    num_na = [int(v['NumNA']) for i, v in df.iterrows()]
    print('NA length')
    print(create_binned_seq_lens(num_na))
    print('=' * 40)

    num_aa = [int(v['NumAA']) for i, v in df.iterrows()]
    print('AA length')
    print(create_binned_seq_lens(num_aa))
    print('=' * 40)

    pcnt_ident = [float(v['pcnt_id']) for i, v in df.iterrows()]
    print('PcntIDs')
    print(create_binned_pcnts(pcnt_ident))
    print('=' * 40)

    print('\n\n', '*' * 40, '\n\n')
