from Utilities import count_number
from Utilities import split_value_by_comma
from Utilities import median_year
from Utilities import int_sorter
from Utilities import create_binnned_year, format_counts_and_percentages
import numpy as np
import pandas as pd
from collections import defaultdict
from Utilities import load_csv
from Utilities import dump_csv


def summarize_pubmed(pubmed, virus_obj):

    summarize_pubmed_reviewer_gpt(pubmed, virus_obj.get_logger('pubmed_workflow'))

    likely = pubmed[
        (
            (pubmed['Reviewer1  (Y/N)'].str.lower().isin(('likely', 'unsure')))
            |
            (pubmed['GPT (Y/N)'].str.lower().isin(('likely', 'unsure')))
        ) &
        (pubmed['Resolve Seq'].str.lower() != 'no') &
        (
            (pubmed['Reviewer(s) Seq'].str.lower() == 'yes') |
            (pubmed['GPT seq (Y/N)'].str.lower() == 'yes')
        )
    ]

    both_unlikely = pubmed[
        (
            (pubmed['Reviewer1  (Y/N)'].str.lower() == 'unlikely') &
            (pubmed['GPT (Y/N)'].str.lower() == 'unlikely')
        )
    ]

    pubmed = pd.concat([likely, both_unlikely])
    logger = virus_obj.get_logger('pubmed_workflow')
    logger.info('Pubmed Literatures:', len(pubmed))
    logger.info('Pubmed Literatures, likely:', len(likely))
    logger.info('Pubmed Literatures, both unlikely:', len(both_unlikely))
    pubmed = virus_obj.process_pubmed(pubmed)

    virus_obj.get_logger('pubmed').report(summarize_pubmed_data(pubmed))

    if virus_obj.pubmed_additional_from_gb:
        additional_pubmed = pd.read_excel(
            virus_obj.pubmed_additional_from_gb, dtype=str).fillna('')

        additional_pubmed['ref_source'] = 'GB reference'

        additional_pubmed = virus_obj.process_pubmed(additional_pubmed)
        logger.info(
            'Only from GenBank Search:',
            len(additional_pubmed))
        pubmed = pd.concat([pubmed, additional_pubmed], ignore_index=True)
        logger.info(
            'Pubmed Literature with additional Literature from GenBank Only:',
            len(pubmed))

    if virus_obj.pubmed_search_missing:
        pubmed_missing = pd.read_excel(
            virus_obj.pubmed_search_missing, dtype=str).fillna('')

        pubmed_missing['ref_source'] = 'Not found in PubMed or GenBank Search'

        pubmed_missing = virus_obj.process_pubmed(pubmed_missing)
        logger.info(
            'Not found in PubMed or GenBank Search:',
            len(pubmed_missing))
        pubmed = pd.concat([pubmed, pubmed_missing], ignore_index=True)

    # IF inlude additional PMID from GenBank
    # if summrize:
    #     summarize_pubmed_data(pubmed, virus_obj.get_logger('pubmed_from_GB'))

    # pubmed['PubID'] = pubmed.index + 1
    get_fixed_Pub_ID(virus_obj, pubmed)
    pubmed.to_excel(virus_obj.pubmed_with_index, index=False)

    return pubmed


def recursive_defaultdict():
    return defaultdict(recursive_defaultdict)


def summarize_pubmed_reviewer_gpt(df, logger):
    summary = recursive_defaultdict()

    summary['Title/Abstract:']['df'] = df

    summary['Title/Abstract:']['R1 likely']['df'] = df[
        (df['Reviewer1  (Y/N)'].str.lower().isin(('likely', 'unsure')))
        ]
    summary['Title/Abstract:']['R1 unlikely']['df'] = df[
        (df['Reviewer1  (Y/N)'].str.lower() == 'unlikely')
        ]
    summary['Title/Abstract:']['GPT likely']['df'] = df[
        (df['GPT (Y/N)'].str.lower().isin(('likely', 'unsure')))
        ]
    summary['Title/Abstract:']['GPT unlikely']['df'] = df[
        (df['GPT (Y/N)'].str.lower() == 'unlikely')]

    summary['Title/Abstract:']['R1 likely']['GPT likely']['df'] = df[
        (df['Reviewer1  (Y/N)'].str.lower().isin(('likely', 'unsure'))) &
        (df['GPT (Y/N)'].str.lower().isin(('likely', 'unsure')))
        ]
    summarize_pubmed_reviewer_gpt_to_full_text(
        summary['Title/Abstract:']['R1 likely']['GPT likely'])

    summarize_pubmed_reviewer_gpt_to_full_text(
        summary['Title/Abstract:']['GPT likely'])

    summary['Title/Abstract:']['R1 likely']['GPT unlikely']['df'] = df[
        (df['Reviewer1  (Y/N)'].str.lower().isin(('likely', 'unsure'))) &
        (df['GPT (Y/N)'].str.lower() == 'unlikely')
        ]
    summary['Title/Abstract:']['R1 likely']['GPT unlikely'][
        'R2 likely']['df'] = df[
            (df['Reviewer1  (Y/N)'].str.lower().isin(('likely', 'unsure'))) &
            (df['GPT (Y/N)'].str.lower() == 'unlikely') &
            (df['Resolve Title'].str.lower().isin(('likely', 'unsure')))
        ]
    summarize_pubmed_reviewer_gpt_to_full_text(
        summary['Title/Abstract:']['R1 likely']['GPT unlikely'][
            'R2 likely'])

    summary['Title/Abstract:']['R1 likely']['GPT unlikely'][
        'R2 unlikely']['df'] = df[
            (df['Reviewer1  (Y/N)'].str.lower().isin(('likely', 'unsure'))) &
            (df['GPT (Y/N)'].str.lower() == 'unlikely') &
            (df['Resolve Title'].str.lower().isin(('unlikely', '')))
    ]
    summarize_pubmed_reviewer_gpt_to_full_text(
        summary['Title/Abstract:']['R1 likely']['GPT unlikely'][
            'R2 unlikely'])

    summary['Title/Abstract:']['R1 unlikely']['GPT likely']['df'] = df[
        (df['Reviewer1  (Y/N)'].str.lower() == 'unlikely') &
        (df['GPT (Y/N)'].str.lower().isin(('likely', 'unsure')))
    ]
    summary['Title/Abstract:']['R1 unlikely']['GPT likely'][
        'R2 likely']['df'] = df[
            (df['Reviewer1  (Y/N)'].str.lower() == 'unlikely') &
            (df['GPT (Y/N)'].str.lower().isin(('likely', 'unsure'))) &
            (df['Resolve Title'].str.lower().isin(('likely', 'unsure')))
    ]
    summarize_pubmed_reviewer_gpt_to_full_text(
        summary['Title/Abstract:']['R1 unlikely']['GPT likely'][
            'R2 likely'])

    summary['Title/Abstract:']['R1 unlikely']['GPT likely'][
        'R2 unlikely']['df'] = df[
            (df['Reviewer1  (Y/N)'].str.lower() == 'unlikely') &
            (df['GPT (Y/N)'].str.lower().isin(('likely', 'unsure'))) &
            (df['Resolve Title'].str.lower().isin(('unlikely', '')))
    ]
    summarize_pubmed_reviewer_gpt_to_full_text(
        summary['Title/Abstract:']['R1 unlikely']['GPT likely'][
            'R2 unlikely'])

    summary['Title/Abstract:']['R1 unlikely']['GPT unlikely']['df'] = df[
        (df['Reviewer1  (Y/N)'].str.lower() == 'unlikely') &
        (df['GPT (Y/N)'].str.lower() == 'unlikely')
        ]

    summarize_pubmed_reviewer_gpt_to_full_text(
        summary['Title/Abstract:']['R1 unlikely']['GPT unlikely'])

    calc_GPT_title_abstract_accuracy(summary)
    calc_R1_title_abstract_accuracy(summary)

    basic_summary(logger, summary)
    logger.info('=' * 80)

    summary = recursive_defaultdict()
    summary['df'] = df[
        (df['Reviewer1  (Y/N)'].str.lower().isin(('likely', 'unsure'))) |
        (df['GPT (Y/N)'].str.lower().isin(('likely', 'unsure')))
    ]
    summarize_pubmed_reviewer_gpt_to_full_text(summary)

    calc_GPT_full_text_accuracy(summary)
    calc_R1_full_text_accuracy(summary)

    basic_summary(logger, summary)
    logger.info('=' * 80)


def calc_GPT_title_abstract_accuracy(summary):

    true_pos = summary['Title/Abstract:']['GPT']['true pos'] = (
        len(summary['Title/Abstract:']['R1 likely']['GPT likely']['df']) +
        len(summary['Title/Abstract:']['R1 unlikely']['GPT likely']['R2 likely']['df'])
    )
    false_neg = summary['Title/Abstract:']['GPT']['false neg'] = (
        len(summary['Title/Abstract:']['R1 likely']['GPT unlikely']['R2 likely']['df'])
        # + summary['Title/Abstract:']['R1 unlikely']['GPT unlikely']['R2 likely']
    )
    false_pos = summary['Title/Abstract:']['GPT']['false pos'] = (
        len(summary['Title/Abstract:']['R1 unlikely']['GPT likely']['R2 unlikely']['df'])
        # summary['Title/Abstract:']['R1 likely']['GPT likely']['R2 unlikely']
    )

    summary['Title/Abstract:']['GPT']['precision'] = 100 * true_pos / (true_pos + false_pos)
    summary['Title/Abstract:']['GPT']['recall'] = 100 * true_pos / (true_pos + false_neg)


def calc_R1_title_abstract_accuracy(summary):

    true_pos = summary['Title/Abstract:']['R1']['true pos'] = (
        len(summary['Title/Abstract:']['R1 likely']['GPT likely']['df']) +
        len(summary['Title/Abstract:']['R1 likely']['GPT unlikely']['R2 likely']['df'])
    )
    false_neg = summary['Title/Abstract:']['R1']['false neg'] = (
        len(summary['Title/Abstract:']['R1 unlikely']['GPT likely']['R2 likely']['df'])
    )
    false_pos = summary['Title/Abstract:']['R1']['false pos'] = (
        len(summary['Title/Abstract:']['R1 likely']['GPT unlikely']['R2 unlikely']['df'])
        # summary['Title/Abstract:']['R1 likely']['GPT likely']['R2 unlikely']
    )

    summary['Title/Abstract:']['R1']['precision'] = 100 * true_pos / (true_pos + false_pos)
    summary['Title/Abstract:']['R1']['recall'] = 100 * true_pos / (true_pos + false_neg)


def calc_GPT_full_text_accuracy(summary):

    true_pos = summary['Full text with seq:']['GPT']['true pos'] = (
        len(summary['Full text with seq:']['R1 yes']['GPT yes']['df']) +
        len(summary['Full text with seq:']['R1 no']['GPT yes']['R2 yes']['df'])
    )
    false_neg = summary['Full text with seq:']['GPT']['false neg'] = (
        len(summary['Full text with seq:']['R1 yes']['GPT no']['R2 yes']['df'])
    )
    false_pos = summary['Full text with seq:']['GPT']['false pos'] = (
        len(summary['Full text with seq:']['R1 no']['GPT yes']['R2 no']['df'])
    )

    summary['Full text with seq:']['GPT']['precision'] = 100 * true_pos / (true_pos + false_pos)
    summary['Full text with seq:']['GPT']['recall'] = 100 * true_pos / (true_pos + false_neg)


def calc_R1_full_text_accuracy(summary):

    true_pos = summary['Full text with seq:']['R1']['true pos'] = (
        len(summary['Full text with seq:']['R1 yes']['GPT yes']['df']) +
        len(summary['Full text with seq:']['R1 yes']['GPT no']['R2 yes']['df'])
    )
    false_neg = summary['Full text with seq:']['R1']['false neg'] = (
        len(summary['Full text with seq:']['R1 no']['GPT yes']['R2 yes']['df'])
    )
    false_pos = summary['Full text with seq:']['R1']['false pos'] = (
        len(summary['Full text with seq:']['R1 yes']['GPT no']['R2 no']['df'])
    )

    summary['Full text with seq:']['R1']['precision'] = 100 * true_pos / (true_pos + false_pos)
    summary['Full text with seq:']['R1']['recall'] = 100 * true_pos / (true_pos + false_neg)


def summarize_pubmed_reviewer_gpt_to_full_text(summary):

    df = summary['df']

    summary['Full text with seq:']['R1 yes']['df'] = df[
        (df['Reviewer(s) Seq'].str.lower() == 'yes')
    ]

    summary['Full text with seq:']['R1 yes']['GPT yes']['df'] = df[
        (df['Reviewer(s) Seq'].str.lower() == 'yes') &
        (df['GPT seq (Y/N)'].str.lower() == 'yes')
        ]

    summary['Full text with seq:']['R1 yes']['GPT no']['df'] = df[
        (df['Reviewer(s) Seq'].str.lower() == 'yes') &
        (df['GPT seq (Y/N)'].str.lower() == 'no')
        ]

    summary['Full text with seq:']['R1 yes']['GPT no'][
            'R2 yes']['df'] = df[
        (df['Reviewer(s) Seq'].str.lower() == 'yes') &
        (df['GPT seq (Y/N)'].str.lower() == 'no') &
        (df['Resolve Seq'].str.lower() == 'yes')
        ]

    summary['Full text with seq:']['R1 yes']['GPT no'][
            'R2 no']['df'] = df[
        (df['Reviewer(s) Seq'].str.lower() == 'yes') &
        (df['GPT seq (Y/N)'].str.lower() == 'no') &
        (df['Resolve Seq'].str.lower().isin(('no', '')))
        ]

    summary['Full text with seq:']['R1 no']['df'] = df[
        (df['Reviewer(s) Seq'].str.lower() == 'no')
    ]

    summary['Full text with seq:']['R1 no']['GPT yes']['df'] = df[
        (df['Reviewer(s) Seq'].str.lower() == 'no') &
        (df['GPT seq (Y/N)'].str.lower() == 'yes')
        ]

    summary['Full text with seq:']['R1 no']['GPT yes'][
            'R2 yes']['df'] = df[
        (df['Reviewer(s) Seq'].str.lower() == 'no') &
        (df['GPT seq (Y/N)'].str.lower() == 'yes') &
        (df['Resolve Seq'].str.lower() == 'yes')
        ]

    summary['Full text with seq:']['R1 no']['GPT yes'][
            'R2 no']['df'] = df[
        (df['Reviewer(s) Seq'].str.lower() == 'no') &
        (df['GPT seq (Y/N)'].str.lower() == 'yes') &
        (df['Resolve Seq'].str.lower().isin(('no', '')))
        ]

    summary['Full text with seq:']['R1 no']['GPT no']['df'] = df[
        (df['Reviewer(s) Seq'].str.lower() == 'no') &
        (df['GPT seq (Y/N)'].str.lower() == 'no')
        ]


def basic_summary(logger, summary, prefix=[]):

    for k, v in summary.items():
        if isinstance(v, pd.DataFrame):
            logger.info(f"{' ' * 4 * len(prefix)}{len(v)}")
            logger.info('-' * 80)
        elif isinstance(v, int):
            logger.info(f"{' ' * 4 * len(prefix)}{k}:{v}")
            logger.info('-' * 80)
        elif isinstance(v, float):
            logger.info(f"{' ' * 4 * len(prefix)}{k}:{v}")
            logger.info('-' * 80)
        elif isinstance(v, str):
            print(v)
        else:
            logger.info(f"{' ' * 4 * len(prefix)}{k}")
            basic_summary(logger, v, prefix + [k])


def summarize_pubmed_data(df):
    summarize_report = []

    section = ["Summarize PubMed"]
    summarize_report.append(section)

    section = ['# Sequences']
    section.append(sum([int(v['NumSeqs']) if v['NumSeqs'] else 0 for i, v in df.iterrows()]))
    summarize_report.append(section)

    section = ["Publish Year"]
    year = count_number([v for i, v in df.iterrows()], 'Year', sorter=int_sorter)
    publish_year = [
        int(v['Year']) for i, v in df.iterrows()
        if v['Year'] and v['Year'] != 'NA']

    counts_formatted, percentages_formatted = format_counts_and_percentages(year)
    section.append(f"Counts:\n{counts_formatted}\n")
    section.append(f"Percentages:\n{percentages_formatted}\n")

    section.append(
        create_binnned_year(publish_year) if publish_year else ''
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
    counts_formatted, percentages_formatted = format_counts_and_percentages(hosts)
    section.append(f"Counts:\n{counts_formatted}\n")
    section.append(f"Percentages:\n{percentages_formatted}\n")
    summarize_report.append(section)

    section = ["Specimen"]
    specimen = count_number([v for i, v in df.iterrows()], 'Specimen')
    counts_formatted, percentages_formatted = format_counts_and_percentages(specimen)
    section.append(f"Counts:\n{counts_formatted}\n")
    section.append(f"Percentages:\n{percentages_formatted}\n")
    summarize_report.append(section)

    section = ["Median of Sample Year"]
    df.loc[:, 'MedianYear'] = df['SampleYr'].apply(median_year)
    year = count_number([v for i, v in df.iterrows()], 'MedianYear', sorter=int_sorter)
    counts_formatted, percentages_formatted = format_counts_and_percentages(year)
    section.append(f"Counts:\n{counts_formatted}\n")
    section.append(f"Percentages:\n{percentages_formatted}\n")

    year = [int(v['MedianYear']) for i, v in df.iterrows() if v['MedianYear'] and v['MedianYear'] != 'NA']
    if year:
        section.append(('Median IQR', np.percentile(year, 25), np.percentile(year, 50), np.percentile(year, 75)))
        section.append(create_binnned_year(year))
    else:
        section.append(('Median IQR', ))
        section.append('')
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
    counts_formatted, percentages_formatted = format_counts_and_percentages(country)
    section.append(f"Counts:\n{counts_formatted}\n")
    section.append(f"Percentages:\n{percentages_formatted}\n")
    section.append(country)
    summarize_report.append(section)

    # section = ["Country W/WO"]
    # country = count_number(
    #     [v for i, v in df.iterrows()], 'Country', translater=with_country)
    # section.append(country)
    # summarize_report.append(section)

    section = ["Gene"]
    gene_list = split_value_by_comma(df, 'Gene')
    gene = count_number(gene_list)
    counts_formatted, percentages_formatted = format_counts_and_percentages(gene)
    section.append(f"Counts:\n{counts_formatted}\n")
    section.append(f"Percentages:\n{percentages_formatted}\n")
    summarize_report.append(section)

    section = ["Seq method"]
    methods = count_number([v for i, v in df.iterrows()], 'SeqMethod')
    counts_formatted, percentages_formatted = format_counts_and_percentages(methods)
    section.append(f"Counts:\n{counts_formatted}\n")
    section.append(f"Percentages:\n{percentages_formatted}\n")
    summarize_report.append(section)

    section = ["End of Report"]
    summarize_report.append(section)

    return summarize_report


def get_fixed_Pub_ID(virus, references):
    if virus.fixed_pub_id_file.exists():
        fixed_ref = load_csv(virus.fixed_pub_id_file)
    else:
        fixed_ref = []

    def find_fixed_ref_id(ref):
        for prev_ref in fixed_ref:
            if prev_ref['PMID'].strip() and ref['PMID'].strip() and prev_ref['PMID'].strip() == ref['PMID'].strip():
                return int(prev_ref['PubID'])

        return None

    max_ref_id = max([int(r['PubID']) for r in fixed_ref]) if fixed_ref else 0

    for idx, ref in references.iterrows():

        fixed_ref_id = find_fixed_ref_id(ref)
        if not fixed_ref_id:
            max_ref_id += 1
            fixed_ref_id = max_ref_id
            fixed_ref.append({
                'PubID': fixed_ref_id,
                'PMID': ref['PMID'],
            })

        references.at[idx, 'PubID'] = fixed_ref_id

    references["PubID"] = references["PubID"].astype(int)

    dump_csv(virus.fixed_pub_id_file, fixed_ref)

    return references
