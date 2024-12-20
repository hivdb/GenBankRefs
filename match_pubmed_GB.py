import pandas as pd
from collections import defaultdict
import Levenshtein
import pandas as pd
import re

from Utilities import count_number
from Utilities import int_sorter
from Utilities import get_values_of_value_count_list
from Utilities import create_binnned_year

from summarize_genbank import summarize_genbank_by_ref
from summarize_pubmed import summarize_pubmed_data
from summarize_genbank import summarize_genbank_full_genome
from summarize_genbank import summarize_genbank_by_seq

from openpyxl import load_workbook
from openpyxl.styles import Alignment
from DataFrameLogic import merge_feature_rows


def match_pubmed_GB(
        pubmed, genbank_ref, genbank_feature, genbank_genes, virus_obj):

    logger = virus_obj.get_logger('compare_matched')
    pubmed_match, genbank_match, pubmed_unmatch, genbank_unmatch = match(
        pubmed, genbank_ref, logger)

    # pubmed_unmatch.to_excel('pubmed_unmatch.xlsx')
    # pd.DataFrame(genbank_unmatch).to_excel('genbank_unmatch.xlsx')

    logger = virus_obj.get_logger('compare_pubmed_only')
    logger.info('Pumbed only Literatures:', len(pubmed_unmatch))
    logger.report(summarize_pubmed_data(pubmed_unmatch))

    logger = virus_obj.get_logger('compare_genbank_only')
    logger.info('GenBank only References:', len(genbank_unmatch))
    logger.report(summarize_genbank_by_ref(genbank_unmatch))

    combined, columns = combine_file(
        pubmed_match, pubmed_unmatch, genbank_unmatch,
        genbank_feature, genbank_genes)

    combined.to_excel(str(virus_obj.pubmed_genbank_combined),
                      index=False, columns=columns)
    format_table(str(virus_obj.pubmed_genbank_combined))

    logger = virus_obj.get_logger('compare_matched')
    summarize_combined_data(
        combined, genbank_feature, genbank_genes, logger)

    logger.report(summarize_genbank_full_genome(
        genbank_match, genbank_feature,
        virus_obj.GENES))

    return pubmed, pubmed_match


def match(pubmed, genbank, logger):

    match_by_title_list = []
    match_by_pmid_list = []
    match_by_acc_list = []
    matched_pubmed_indices = []
    genbank_unmatch_list = []

    for index, row in genbank.iterrows():
        pmid = row['PMID']

        result = pubmed[pubmed['PMID'] == pmid]

        match_by_acc = None
        match_by_pmid = None
        match_by_title = None

        if not result.empty:
            matched_pubmed_indices.extend(result.index.tolist())
            match_by_pmid = [row, result, index]
            match_by_pmid_list.append(match_by_pmid)
            # continue

        title = row['Title'].replace('Direct Submission,', '').replace(', Direct Submission', '').strip()

        # Pubmed title always exists
        result = pubmed[pubmed['Title'].apply(
            lambda x: Levenshtein.distance(x.lower(), title.lower()) < 5)]

        if not result.empty:
            matched_pubmed_indices.extend(result.index.tolist())
            match_by_title = [row, result, index]
            match_by_title_list.append(match_by_title)

        accession_list = row['accession']
        accession_prefix_list = set([
            a.strip()[:6]
            for a in accession_list.split(',')
        ])

        result = search_access_prefix(pubmed, accession_prefix_list)
        if not result.empty:
            matched_pubmed_indices.extend(result.index.tolist())
            match_by_acc = [row, result, index]
            match_by_acc_list.append(match_by_acc)

        if match_by_acc or match_by_pmid or match_by_title:
            pass
        else:
            genbank_unmatch_list.append(row)

    genbank_match_list = match_by_title_list + match_by_pmid_list + match_by_acc_list

    logger.info("Genbank match by pmid:", len(set(i[-1] for i in match_by_pmid_list)))
    logger.info("Genbank match by title:", len(set(i[-1] for i in match_by_title_list)))
    logger.info("Genbank match by acc:", len(set(i[-1] for i in match_by_acc_list)))
    logger.info('Genbank match total:', len(set(i[-1] for i in genbank_match_list)))
    logger.info('-' * 80)

    pubmed_match = match_pubmed2genbank(genbank_match_list)

    genbank_match = {}
    for row, result, index in genbank_match_list:
        genbank_match[index] = row

    genbank_match = genbank_match.values()

    logger.info("Pubmed match by pmid:", len(match_pubmed2genbank(match_by_pmid_list)))
    logger.info("Pubmed match by title:", len(match_pubmed2genbank(match_by_title_list)))
    logger.info("Pubmed match by acc:", len(match_pubmed2genbank(match_by_acc_list)))
    logger.info('Pubmed match total:', len(pubmed_match))

    pubmed_unmatch = pubmed.drop(index=matched_pubmed_indices)

    return pubmed_match, pd.DataFrame(genbank_match), pubmed_unmatch, pd.DataFrame(genbank_unmatch_list)


def search_access_prefix(pubmed, accession_prefix_list):
    found = pd.DataFrame()

    for acc_prefix in accession_prefix_list:
        result = pubmed[pubmed['GenBank'].apply(
            lambda x: acc_prefix.lower() in str(x).lower())]
        if not result.empty:
            found = pd.concat([found, result])

    return found


def match_pubmed2genbank(genbank_match):
    pubmed_match = defaultdict(dict)
    for g, pubmed_list, index in genbank_match:
        for r, i in pubmed_list.iterrows():
            pubmed_match[r]['pubmed'] = i
            if 'genbank' not in pubmed_match[r]:
                pubmed_match[r]['genbank'] = []

            pubmed_match[r]['genbank'].append(g)

    pubmed_match = [
        [v['pubmed'], v['genbank']]
        for v in pubmed_match.values()
    ]

    return pubmed_match


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

    # section = ['Journals']
    # journals = count_number([v for i, v in matches.iterrows()], 'Journal')
    # section.append(journals)
    # summarize_report.append(section)

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

    section = ['# Matched pubmed genbank']
    summarize_report.append(section)

    logger.report(summarize_report)

    logger.report(summarize_genbank_by_seq(features, genes))


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


def combine_file(
        pubmed_match, pubmed_unmatch, genbank_unmatch,
        features_df, genes_df
        ):

    result = []
    for pubmed, genbank_list in pubmed_match:
        accessions = set([
             j.strip()
             for i in genbank_list
             for j in i['accession'].split(',')
        ])

        features = features_df[features_df['Accession'].isin(
            accessions)]

        genes = genes_df[genes_df['Accession'].isin(accessions)]

        features_stat = merge_feature_rows(features, genes)

        row = {
            'Authors': pubmed['Authors'],
            'Title': pubmed['Title'],
            'Journal': pubmed['Journal'],
            'Year': pubmed['Year'],
            'PMID': pubmed['PMID'],

            'Reviewer(s) Seq': pubmed['Reviewer(s) Seq'],
            'GPT seq (Y/N)': pubmed['GPT seq (Y/N)'],
            'Resolve Title': pubmed['Resolve Title'],
            'match': 'Yes',

            'Viruses (PM)': pubmed['Viruses'],
            'NumSeqs (PM)': pubmed['NumSeqs'],
            'Hosts (PM)': pubmed['Host'],
            'Specimen (PM)': pubmed['Specimen'],

            'SampleYr (PM)': pubmed['SampleYr'],
            'Countries (PM)': pubmed['Country'],

            'Genes (PM)': pubmed['Gene'],
            'SeqMethod (PM)': pubmed['SeqMethod'],
            'CloneMethod (PM)': pubmed['CloneMethod'],
            'GenBank (PM)': pubmed['GenBank'],

            'Viruses (GB)': features_stat['Organisms'],
            'NumSeqs (GB)': len(accessions),
            'Hosts (GB)': features_stat['Hosts'],
            'Specimen (GB)': features_stat['Specimens'],

            'SampleYr (GB)':  features_stat['IsolateYears'],
            'Countries (GB)': features_stat['Countries'],

            'Genes (GB)': features_stat['Gene'],
            'SeqMethod (GB)': '',
            'CloneMethod (GB)': '',
            'GenBank (GB)': ', '.join(sorted(list(accessions))),

            'NumSubSeqs': features_stat['NumSubSeqs'],
            'AlignLens (GB)': features_stat['AlignLens'],
            'PcntIDs (GB)': features_stat['PcntIDs'],
        }

        result.append(row)

    for r, pubmed in pubmed_unmatch.iterrows():
        row = {
            'Authors': pubmed['Authors'],
            'Title': pubmed['Title'],
            'Journal': pubmed['Journal'],
            'Year': pubmed['Year'],
            'PMID': pubmed['PMID'],

            'Reviewer(s) Seq': pubmed['Reviewer(s) Seq'],
            'GPT seq (Y/N)': pubmed['GPT seq (Y/N)'],
            'Resolve Title': pubmed['Resolve Title'],

            'Viruses (PM)': pubmed['Viruses'],
            'NumSeqs (PM)': pubmed['NumSeqs'],
            'Hosts (PM)': pubmed['Host'],
            'Specimen (PM)': pubmed['Specimen'],
            'SampleYr (PM)': pubmed['SampleYr'],
            'Countries (PM)': pubmed['Country'],
            'Genes (PM)': pubmed['Gene'],
            'SeqMethod (PM)': pubmed['SeqMethod'],
            'CloneMethod (PM)': pubmed['CloneMethod'],
            'GenBank (PM)': pubmed['GenBank'],
        }

        result.append(row)

    for row, genbank in genbank_unmatch.iterrows():

        accessions = set(genbank['accession'].split(','))
        features = features_df[features_df['Accession'].isin(
            accessions)]

        genes = genes_df[genes_df['Accession'].isin(accessions)]

        features_stat = merge_feature_rows(features, genes)

        row = {
            'Authors': genbank['Authors'],
            'Title': genbank['Title'],
            'Journal': genbank['Journal'],
            'Year': genbank['Year'],
            'PMID': genbank['PMID'],
            'Viruses (GB)': features_stat['Organisms'],
            'NumSeqs (GB)': len(accessions),
            'Hosts (GB)': features_stat['Hosts'],
            'Specimen (GB)': features_stat['Specimens'],
            'SampleYr (GB)': features_stat['IsolateYears'],
            'Countries (GB)': features_stat['Countries'],
            'Genes (GB)': features_stat['Gene'],
            'SeqMethod (GB)': '',
            'CloneMethod (GB)': '',
            'GenBank (GB)': ', '.join(sorted(list(accessions))),
            'NumSubSeqs': features_stat['NumSubSeqs'],
            'AlignLens (GB)': features_stat['AlignLens'],
            'PcntIDs (GB)': features_stat['PcntIDs'],
        }

        result.append(row)

    columns = [
        'Authors',
        'Title',
        'Journal',
        'Year',
        'PMID',

        'Reviewer(s) Seq',
        'GPT seq (Y/N)',
        'Resolve Title',
        'match',

        'Viruses (PM)',
        'Viruses (GB)',

        'NumSeqs (PM)',
        'NumSeqs (GB)',

        'Hosts (PM)',
        'Hosts (GB)',

        'Specimen (PM)',
        'Specimen (GB)',

        'SampleYr (PM)',
        'SampleYr (GB)',

        'Countries (PM)',
        'Countries (GB)',

        'Genes (PM)',
        'Genes (GB)',

        'SeqMethod (PM)',
        'SeqMethod (GB)',

        'CloneMethod (PM)',
        'CloneMethod (GB)',

        'GenBank (PM)',
        'GenBank (GB)',

        'NumSubSeqs',
        'AlignLens (GB)',
        'PcntIDs (GB)',
    ]

    for i in result:
        for c in columns:
            if c not in i:
                i[c] = ''

    return pd.DataFrame(result), columns


def format_table(excel_file):
    wb = load_workbook(excel_file)
    ws = wb.active

    for row in ws.iter_rows(
            min_row=1, max_row=ws.max_row, min_col=1, max_col=ws.max_column):
        for cell in row:
            cell.alignment = Alignment(
                horizontal="left", vertical="top", wrap_text=True)

    for col in ws.columns:
        col_letter = col[0].column_letter
        ws.column_dimensions[col_letter].width = 20

    wb.save(excel_file)
