from .utils import count_number
from .utils import split_value_by_comma
from .translate_value import median_year
from .utils import int_sorter
from Utilities import create_binnned_year
import numpy as np
import pandas as pd
from collections import defaultdict


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


def summarize_pubmed(df, logger):
    summarize_report = []

    section = ["Summarize PubMed"]
    summarize_report.append(section)

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
