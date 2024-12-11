from .utils import count_number
from .utils import merge_genbank_list_columns
from .utils import int_sorter
from .sum_genbank import summarize_genbank_by_seq, summarize_genbank_full_genome
from .utils import split_value_by_comma
from .utils import get_values_of_value_count_list
from Utilities import create_binnned_year


def summarize_combined_data(combined, genbank):
    print('Matched')
    matches = combined[(combined['match'] == 'Yes')]

    publish_year = count_number([v for i, v in matches.iterrows()], 'Year', sorter=int_sorter)
    print('Publish Year')
    print(publish_year)
    publish_year = [
        int(v['Year']) for i, v in matches.iterrows()
        if v['Year'] and v['Year'] != 'NA']
    print(create_binnned_year(publish_year))
    print('=' * 40)

    journals = count_number([v for i, v in matches.iterrows()], 'Journal')
    print('Journals')
    print(journals)
    print('=' * 40)

    methods = count_number(
        [v for i, v in matches.iterrows()], 'SeqMethod (PM)')
    print('Seq method')
    print(methods)
    print('=' * 40)

    accessions = set([
        j.strip()
        for i, v in matches.iterrows()
        for j in v['GenBank (GB)'].split(',')
        ])
    num_seq = len(accessions)
    print('NumSeq', num_seq)
    print('=' * 40)

    genbank['Accession'] = genbank['Accession'].apply(lambda x: x.strip().split('.')[0])
    genbank = genbank[genbank['Accession'].isin(list(accessions))]

    summarize_genbank_by_seq(genbank)

    print('Similar virus')
    print(summarize_similarity(combined, 'Viruses'))

    print('Similar hosts')
    print(summarize_similarity(combined, 'Hosts'))

    print('Similar Specimens')
    print(summarize_similarity(combined, 'Specimen'))

    print('Similar countries')
    print(summarize_similarity(combined, 'Countries'))

    print('Similar Genes')
    print(summarize_similarity(combined, 'Genes'))


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

