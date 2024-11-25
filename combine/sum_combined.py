from .utils import count_number
from .utils import merge_genbank_list_columns

def summarize_combined_data(combined):
    matches = combined[(combined['match'] == 'Yes')]

    publish_year = count_number([v for i, v in matches.iterrows()], 'Year')
    print('Publish Year')
    print(publish_year)
    print('=' * 40)

    journals = count_number([v for i, v in matches.iterrows()], 'Journal')
    print('Journals')
    print(journals)
    print('=' * 40)

    num_seq = sum([int(v['NumSeqs (GB)']) for i, v in matches.iterrows()])
    print('NumSeq', num_seq)
    print('=' * 40)

    countries = merge_genbank_list_columns(
        [v for i, v in matches.iterrows()], 'Countries (GB)')
    print('Countries')
    print(countries)
    print('=' * 40)

    year = merge_genbank_list_columns(
        [v for i, v in matches.iterrows()], 'SampleYr (GB)')
    print("Sample year")
    print(year)
    print('=' * 40)

    host = merge_genbank_list_columns(
        [v for i, v in matches.iterrows()], 'Hosts (GB)')
    print('Host')
    print(host)
    print('=' * 40)

    specimens = merge_genbank_list_columns(
        [v for i, v in matches.iterrows()], 'Specimen (GB)')
    print('Specimens')
    print(specimens)
    print('=' * 40)

    genes = merge_genbank_list_columns(
        [v for i, v in matches.iterrows()], 'Genes (GB)')
    print('Genes')
    print(genes)
    print('=' * 40)

    aligns = merge_genbank_list_columns(
        [v for i, v in matches.iterrows()], 'AlignLens (GB)')
    print('AlignLens')
    print(aligns)
    print('=' * 40)

    pcnt_ident = merge_genbank_list_columns(
        [v for i, v in matches.iterrows()], 'PcntIDs (GB)')
    print('PcntIDs')
    print(pcnt_ident)
    print('=' * 40)

    methods = count_number(
        [v for i, v in matches.iterrows()], 'SeqMethod (PM)')
    print('Seq method')
    print(methods)
    print('=' * 40)

    print('\n\n', '*' * 40, '\n\n')




