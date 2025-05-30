from Utilities import median_year
from Utilities import with_country
from Utilities import count_number
from Utilities import int_sorter
from Utilities import create_binned_pcnts
from Utilities import create_binned_seq_lens
from Utilities import create_binnned_year, count_rev_sorter
from Utilities import format_counts_and_percentages
from collections import defaultdict, Counter
from Bio import SeqIO, pairwise2
import pandas as pd

import re


def summarize_genbank(genbank_ref, genbank_feature, genbank_genes, virus_obj):

    genbank_ref = genbank_ref.fillna('')

    genbank_feature = genbank_feature.fillna('')

    genbank_genes = genbank_genes.fillna('')

    logger = virus_obj.get_logger('genbank')
    logger.report(summarize_genbank_by_ref(genbank_ref))

    logger.report(summarize_genbank_by_seq(virus_obj, genbank_feature, genbank_genes))

    logger.report(summarize_genbank_full_genome(
        genbank_ref, genbank_feature,
        full_gene_set=virus_obj.GENES))


def summarize_genbank_by_ref(df):
    summarize_report = []

    section = ['Summarize Genbank By Ref']
    summarize_report.append(section)

    df['MedianPublishYear'] = df['Year'].apply(median_year)
    publish_year = count_number([
        v for i, v in df.iterrows()], 'MedianPublishYear', sorter=int_sorter)
    section = ['Publish Year']
    publish_year = '\n'.join([
        f'{k} ({v})'
        for k, v in count_rev_sorter(publish_year.items())
    ])
    section.append(publish_year)

    publish_year = [
        int(v['MedianPublishYear']) for i, v in df.iterrows()
        if v['MedianPublishYear'] and v['MedianPublishYear'] != 'NA']
    section.append(create_binnned_year(publish_year))

    summarize_report.append(section)

    # Journal information not included
    # journal_values = [row['Journal'].split(',')[0].strip()
    #                   for _, row in df.iterrows() if 'Journal' in row and pd.notnull(row['Journal'])]
    # cleaned_entries = [remove_parenthesis(entry) for entry in journal_values]
    # journals = count_number([{'Journal': value}
    #                         for value in cleaned_entries], 'Journal')
    # section = ['Journals']
    # section.append(journals)

    section = ['References Num Isolates (Sequences)']
    df['NumSeq (GB)'] = df['accession'].apply(lambda x: len(x.split(',')))
    num_seqs = count_number(
        [v for i, v in df.iterrows()], 'NumSeq (GB)', sorter=int_sorter)
    num_seqs = '\n'.join([
        f'{k} ({v})'
        for k, v in count_rev_sorter(num_seqs.items())
    ])
    section.append(num_seqs)

    section.append((
        'Total', len(set([
            j.strip()
            for i, v in df.iterrows() if v['NumSeq (GB)']
            for j in v['accession'].split(',')
            if j.strip()
        ]))))
    summarize_report.append(section)

    section = ["End of Report"]
    summarize_report.append(section)

    return summarize_report


def summarize_genbank_by_seq(virus_obj, df, genes_df, file_prefix=''):
    summarize_report = []

    if virus_obj.name == 'CCHF':
        # filter out non_CCHF
        df = df[df['Description'].str.contains('orthonairovirus|crimean-congo', case=False, na=False)]

    section = [f'Total after dropping non-clinical and non_virus: {df.shape[0]}']
    summarize_report.append(section)

    section = ['Summarize Genbank By Seq']
    summarize_report.append(section)

    print('# special accesions', virus_obj.special_accessions)

    section = ['Host']
    hosts = count_number([v for i, v in df.iterrows()], 'Host')
    counts_formatted, percentages_formatted = format_counts_and_percentages(hosts)

    section.append(f"Counts:\n{counts_formatted}\n")
    section.append(f"Percentages:\n{percentages_formatted}\n")
    summarize_report.append(section)

    section = ['Specimens']
    specimen = count_number([v for i, v in df.iterrows()], 'isolate_source')
    counts_formatted, percentages_formatted = format_counts_and_percentages(specimen)

    section.append(f"Counts:\n{counts_formatted}\n")
    section.append(f"Percentages:\n{percentages_formatted}\n")
    summarize_report.append(section)

    section = ['RecordYears']
    year = count_number(
        [v for i, v in df.iterrows()], 'RecordYear', sorter=int_sorter)
    counts_formatted, percentages_formatted = format_counts_and_percentages(year)

    section.append(f"Counts:\n{counts_formatted}\n")
    section.append(f"Percentages:\n{percentages_formatted}\n")

    year = [int(v['RecordYear']) for i, v in df.iterrows() if v['RecordYear']]
    section.append(create_binnned_year(year))
    summarize_report.append(section)

    section = ['Sample Years']
    year = count_number([
        int(v['IsolateYear'])
        if v['IsolateYear'] else '' for i, v in df.iterrows()],
        sorter=int_sorter)
    counts_formatted, percentages_formatted = format_counts_and_percentages(year)

    section.append(f"Counts:\n{counts_formatted}\n")
    section.append(f"Percentages:\n{percentages_formatted}\n")

    year = [int(v['IsolateYear']) for i, v in df.iterrows() if v['IsolateYear'] and v['IsolateYear'] != 'NA']
    section.append(create_binnned_year(year))
    summarize_report.append(section)

    section = ['Countries']
    country = count_number(
        [v for i, v in df.iterrows()], 'Country')
    counts_formatted, percentages_formatted = format_counts_and_percentages(country)

    section.append(f"Counts:\n{counts_formatted}\n")
    section.append(f"Percentages:\n{percentages_formatted}\n")
    summarize_report.append(section)

    section = ['Countries W/WO']
    country = count_number(
        [v for i, v in df.iterrows()], 'Country',
        translater=with_country)
    section.append(country)
    section.append('=' * 40)

    section = ['Genes']

    df[df['Genes'] == 'NA'].to_excel(f"OutputData/{virus_obj.name}/excels/isolates_without_genes.xlsx", index=False)
    genes = Counter()
    for row in df['Genes']:
        unique_genes = set(row.split(', '))
        genes.update(unique_genes)

    total_count = df[df['Genes'] != 'NA'].shape[0]
    counts_formatted, percentages_formatted = format_counts_and_percentages(genes, total=total_count)

    # get number of entries that has 1 gene, 2 gene, 3 gene
    num_gene_counts = df.loc[df['Genes'] != 'NA', 'Genes'].apply(
        lambda x: len(set(x.split(','))) if isinstance(x, str) else 0
    )
    num_gene_distribution = Counter(num_gene_counts)
    num_gene_distribution_f, num_gene_percent_f = format_counts_and_percentages(num_gene_distribution, total=total_count)

    section.append(f"Counts:\n{counts_formatted}\n")
    section.append(f"Percentages:\n{percentages_formatted}\n")
    section.append(f"Counts # of genes:\n{num_gene_distribution_f}\n")
    section.append(f"Percentages # of genes:\n{num_gene_percent_f}\n")

    final_df = merge_segements(virus_obj, df, file_prefix)

    total_count_after_na = final_df[final_df['Genes'] != 'NA'].shape[0]
    combined_num_gene_counts = final_df.loc[final_df['Genes'] != 'NA', 'Genes'].apply(
        lambda x: len(set(x.split(','))) if isinstance(x, str) else 0
    )
    combined_num_gene_distribution = Counter(combined_num_gene_counts)
    combined_num_gene_distribution_f, combined_num_gene_percent_f = format_counts_and_percentages(combined_num_gene_distribution, total=total_count_after_na)

    combined_genes = Counter()
    for row in final_df['Genes']:
        unique_genes = set(row.split(', '))
        combined_genes.update(unique_genes)
    counts_formatted_combined, percentages_formatted_combined = format_counts_and_percentages(combined_genes, total=total_count_after_na)

    # # Check merge is correct, isolateName often the same for different isolates - incorrect for Nipah
    # merged_groups = df_to_merge.groupby(['country_region', 'collection_date', 'Host', 'IsolateName']).filter(lambda x: len(x) > 1)

    section.append(f"Counts genes after combining isolate:\n{counts_formatted_combined}\n")
    section.append(f"Percentages genes after combining isolate:\n{percentages_formatted_combined}\n")
    section.append(f"Counts # of genes after combining isolate:\n{combined_num_gene_distribution_f}\n")
    section.append(f"Percentages # of genes after combining isolate:\n{combined_num_gene_percent_f}\n")

    summarize_report.append(section)

    for i in virus_obj.GENES:
        summarize_report.append([f"Gene {i}"])

        sub_gene_df = genes_df[genes_df['Gene'] == i]

        section = ['AlignLens']
        aligns = [int(v['align_len']) for i, v in sub_gene_df.iterrows()]
        section.append(create_binned_seq_lens(aligns))
        summarize_report.append(section)

        section = ['NA length']
        num_na = [int(v['NA_length']) for i, v in sub_gene_df.iterrows()]
        section.append(create_binned_seq_lens(num_na))
        summarize_report.append(section)

        section = ['NA start']
        num_na = [int(v['NA_start']) for i, v in sub_gene_df.iterrows()]
        section.append(create_binned_seq_lens(num_na))
        summarize_report.append(section)

        section = ['NA stop']
        num_na = [int(v['NA_stop']) for i, v in sub_gene_df.iterrows()]
        section.append(create_binned_seq_lens(num_na))
        summarize_report.append(section)

        section = ['AA length']
        num_aa = [int(v['AA_length']) for i, v in sub_gene_df.iterrows()]
        section.append(create_binned_seq_lens(num_aa))
        summarize_report.append(section)

        section = ['PcntIDs']
        pcnt_ident = [float(v['pcnt_id']) for i, v in sub_gene_df.iterrows()]
        section.append(create_binned_pcnts(pcnt_ident))
        summarize_report.append(section)

    section = ['End of report']
    summarize_report.append(section)

    return summarize_report


def summarize_genbank_full_genome(
        ref_df, features_df, full_gene_set):

    summarize_report = []

    section = ['Summarize references with full genomes']
    summarize_report.append(section)

    num_ref = 0
    # num_seq = 0
    for i, row in ref_df.iterrows():
        accessions = set([
            a.strip()
            for a in row['accession'].split(',')
        ])

        features = features_df[
            features_df['Accession'].isin(accessions)]

        gene_count = defaultdict(int)
        for gene_list in features['Genes']:
            for gene in gene_list.split(','):
                gene = gene.strip()
                gene_count[gene] += 1

        genome = 0
        if 'genome' in gene_count:
            genome += gene_count['genome']
            del gene_count['genome']

        if set(gene_count.keys()) == set(full_gene_set):
            genome += min(gene_count.values())

        if genome:
            # num_seq += genome
            num_ref += 1

    section = ['Number of References with Full genome:']
    section.append(num_ref)
    # section.append('Number of full genome seq', num_seq)

    summarize_report.append(section)

    section = ['End of report']
    summarize_report.append(section)

    return summarize_report


def merge_segements(virus_obj, df, file_prefix):
    # Try to merge on isolateName, the fields should all be the same for same isolate
    # Filter out rows where "IsolateName" is not empty and perform merging

    for idx, row in df.iterrows():
        accession_prefix = str(row['Accession'])[:7]
        df.loc[idx, 'Accession_prefix'] = accession_prefix

    df_no_merge = df[df['IsolateName'] == ""].copy()
    # df_no_merge.to_excel(f"OutputData/{virus}/excels/{file_prefix}_G1.xlsx")

    # Count occurrences of each group, if >3 don't merge
    # maybe 6 for Nipah, 3 for cchf?
    # Gene should also be unique
    def count_genes(genes):
        unique_genes = set(
            g.strip()
            for v in genes
            for g in v.split(',')
        )
        return len(unique_genes)

    df_to_merge = df[df['IsolateName'] != ""].copy()
    df_to_merge['count'] = df_to_merge.groupby([
        'Accession_prefix', 'country_region', 'collection_date', 'Host', 'IsolateName'
    ])['Genes'].transform(count_genes)

    valid_groups = df_to_merge[(df_to_merge['count'] > 1)]
    valid_groups.to_excel(f"OutputData/{virus_obj.name}/excels/{virus_obj.name}_{file_prefix}_merge.xlsx")
    invalid_groups = df_to_merge[(df_to_merge['count'] == 1)]
    # Combine invalid groups back into df_no_merge
    df_no_merge = pd.concat([df_no_merge, invalid_groups], ignore_index=True)
    invalid_groups.to_excel(f"OutputData/{virus_obj.name}/excels/{virus_obj.name}_{file_prefix}_no_merge.xlsx", index=False)

    # Merge the valid_groups back with the original DataFrame to retain all columns

    merged_valid_df = (
        valid_groups.groupby(['Accession_prefix', 'country_region', 'collection_date', 'Host', 'IsolateName'], as_index=False)
        .agg({
            'Genes': lambda x: ', '.join(sorted(set(g.strip() for v in x.dropna() for g in v.split(',')))),
            'Accession': lambda x: ', '.join(sorted(set(x.dropna()))),  # Combine unique Accessions
            'cds': lambda x: ', '.join(sorted(set(g.strip() for v in x.dropna() for g in v.split(',')))),  # Combine unique cds
            'isolate_source2': lambda x: ', '.join(sorted(set(x.dropna()))),
            **{col: 'first' for col in valid_groups.columns if col not in
                ['country_region', 'collection_date', 'Host', 'IsolateName',
                'Genes', 'Accession', 'cds', 'isolate_source2',
                'Description', 'record_date', 'organism', 'segment_source']}
        })
    )

    # Combine both datasets
    final_df = pd.concat([merged_valid_df, df_no_merge], ignore_index=True)

    # print(df.shape[0], final_merged_df.shape[0])
    if file_prefix:
        final_df.to_excel(f"OutputData/{virus_obj.name}/excels/{virus_obj.name}_{file_prefix}_merge_result.xlsx")
    else:
        final_df.to_excel(virus_obj.isolate_file)

    return final_df
