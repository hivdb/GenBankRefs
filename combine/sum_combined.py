from .utils import count_number
from .utils import int_sorter
from .sum_genbank import summarize_genbank_by_seq
from .utils import get_values_of_value_count_list
from Utilities import create_binnned_year
import re


def summarize_combined_data(combined, features, genes, logger):
    summarize_report = []

    section = ['Summarize PubMed GenBank Match']
    summarize_report.append(section)

    genbank_only_pubmed = combined[
        (combined['match'] != 'Yes') &
        (combined['Reviewer(s) Seq'] == '') &
        (combined['PMID'] != '')
        ]
    section = ['GenBank only PMID']
    section.append(len(genbank_only_pubmed))
    summarize_report.append(section)

    matches = combined[(combined['match'] == 'Yes')]

    section = ['Publish Year']
    publish_year = count_number([
        v for i, v in matches.iterrows()], 'Year', sorter=int_sorter)
    section.append(publish_year)

    publish_year = [
        int(v['Year']) for i, v in matches.iterrows()
        if v['Year'] and v['Year'] != 'NA']
    section.append(create_binnned_year(publish_year))
    summarize_report.append(section)

    section = ['Journals']
    journals = count_number([v for i, v in matches.iterrows()], 'Journal')
    section.append(journals)
    summarize_report.append(section)

    section = ['Seq method']
    methods = count_number(
        [v for i, v in matches.iterrows()], 'SeqMethod (PM)')
    section.append(methods)
    summarize_report.append(section)

    section = ['After matching, NumSeq from GenBank']
    accessions = set([
        j.strip()
        for i, v in matches.iterrows()
        for j in v['GenBank (GB)'].split(',')
        ])
    num_seq = len(accessions)
    section.append(num_seq)
    summarize_report.append(section)

    features = features[features['Accession'].isin(list(accessions))]

    genes = genes[genes['Accession'].isin(list(accessions))]



    section = ['Pubmed Supplement GenBank']
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

        section.append((name, count))
    summarize_report.append(section)

    # section = ['Similar virus']
    # section.append(summarize_similarity(combined, 'Viruses'))

    # section = ['Similar hosts']
    # section.append(summarize_similarity(combined, 'Hosts'))

    # section = ['Similar Specimens']
    # section.append(summarize_similarity(combined, 'Specimen'))

    # section = ['Similar countries']
    # section.append(summarize_similarity(combined, 'Countries'))

    # section = ['Similar Genes']
    # section.append(summarize_similarity(combined, 'Genes'))

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

    section = ['# Matched pubmed genbank']
    summarize_genbank_by_seq(features, genes, logger)


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
