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

    logger.report(summarize_genbank_by_seq(
        genbank_feature, genbank_genes))

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


def match_genes(pattern, description, gene_dict_keys):
    """
    Only return genes in gene_dict_keys (listt)
    """
    matches = pattern.findall(description)
    matched_genes = set()

    for match in matches:
        matched_terms = [m for m in match if m and m in gene_dict_keys]  # Filter only valid genes
        matched_genes.update(matched_terms)

    return ', '.join(matched_genes) if matched_genes else "NA"


def local_align_genes(seq, description, virus_name, acc):

    gene_dict = {}  # Dictionary to store gene names and sequences

    with open(f"ReferenceData/{virus_name}/{virus_name}_RefNAs.fasta", "r") as fasta_file:
        for record in SeqIO.parse(fasta_file, "fasta"):
            gene_dict[record.id] = str(record.seq)

    matched_genes = []
    for gene, ref_seq in gene_dict.items():
        align_score = pairwise2.align.localxx(seq, ref_seq, score_only=True)
        if align_score > len(ref_seq) * 0.8:  # 80% similarity threshold
            # if acc == 'OQ791499':
            #     print('new', gene, align_score)
            #     for i in pairwise2.align.localms(seq, ref_seq, 2, -3, -5, -2):
            #         print(pairwise2.format_alignment(*i))
            matched_genes.append(gene)

    genes = ', '.join(matched_genes) if matched_genes else 'NA'

    # for those still NA, try finding protein name from submission description
    if genes == 'NA':
        unique_genes = list(gene_dict.keys())
        pattern = re.compile(
            r'\b(' + '|'.join(map(re.escape, unique_genes)) + r')\b(?=\s+(protein|gene))'
            r'|\(\b(' + '|'.join(map(re.escape, unique_genes)) + r')\b\)\s+gene',
            re.IGNORECASE
        )

        genes = match_genes(pattern, description, gene_dict.keys())

    # if still not found
    if genes == 'NA':
        matched_genes = []
        description = description.lower()
        if virus_name == 'CCHF':
            # filter out other viruses that are not CCHF
            if not ("orthonairovirus" in description or "crimean-congo" in description):
                return genes
            # attempt matching by gene name in description
            segment_s_keywords = ['nucleocapsid', 'nucleoprotein', 'segment: s', 'segment s']
            if any(w in description for w in segment_s_keywords):
                matched_genes.append("S")
            if 'glycoprotein' in description or "segment m" in description:
                matched_genes.append("M")
            if "rdrp" in description or 'rna polymerase' in description:
                matched_genes.append("L")

        elif virus_name == 'Lassa':
            if not ("lassa" in description):
                return genes
            if 'polymerase' in description:
                matched_genes.append("L")
            if 'nucleoprotein' in description or 'nucleocapsid' in description:
                matched_genes.append("N")
            if 'glycoprotein' in description or 'gpc' in description:
                matched_genes.append("G")
            if 'matrix' in description or '(z)' in description:
                matched_genes.append("Z")

        genes = ', '.join(matched_genes) if matched_genes else 'NA'

    return genes


def summarize_genbank_by_seq(df, genes_df):
    summarize_report = []

    virus = df.iloc[0]['organism']
    if virus == 'CCHF':
        # filter out non_CCHF
        df = df[df['Description'].str.contains('orthonairovirus|crimean-congo', case=False, na=False)]

    section = [f'Total after dropping non-clinical and non_virus: {df.shape[0]}']
    summarize_report.append(section)

    section = ['Summarize Genbank By Seq']
    summarize_report.append(section)

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

    # align sequences where gene is empty to get gene
    for index, row in df.iterrows():
        if row['Genes'] == "": # if we comment this out, more accurate, only detect genes present in seq, currently has dup
            # should move to an earlier step when generating feature_df? some accessions not in genes_df
            # combining by isolateName & removing dup solves the problem partially
            new_genes = local_align_genes(row['Seq'], row['Description'], virus, row['Accession'])
            if (row['Accession'] == 'OQ791499'):
                print(new_genes)
            df.at[index, 'Genes'] = new_genes

    df[df['Genes'] == 'NA'].to_excel(f"OutputData/{virus}/excels/gene_missing.xlsx")
    genes = Counter()
    for row in df['Genes']:
        unique_genes = set(row.split(', '))
        genes.update(unique_genes)

    total_count = df.shape[0]
    counts_formatted, percentages_formatted = format_counts_and_percentages(genes, total=total_count)


    # get number of entries that has 1 gene, 2 gene, 3 gene
    num_gene_counts = df.loc[df['Genes'] != 'NA', 'Genes'].apply(lambda x: len(set(x.split(', '))) if isinstance(x, str) else 0)
    num_gene_distribution = Counter(num_gene_counts)
    num_gene_distribution_f, num_gene_percent_f = format_counts_and_percentages(num_gene_distribution, total=total_count)

    section.append(f"Counts:\n{counts_formatted}\n")
    section.append(f"Percentages:\n{percentages_formatted}\n")
    section.append(f"Counts # of genes:\n{num_gene_distribution_f}\n")
    section.append(f"Percentages # of genes:\n{num_gene_percent_f}\n")

    # Try to merge on isolateName, the fields should all be the same for same isolate
    df_no_merge = df[df['IsolateName'] == ""].copy()

    # Filter out rows where "IsolateName" is not empty and perform merging
    df_to_merge = df[df['IsolateName'] != ""].copy()


    # Count occurrences of each group
    group_counts = df_to_merge.groupby(['country_region', 'collection_date', 'Host', 'IsolateName']).size().reset_index(name='count')

    # Separate groups: ones with <= 3 occurrences (to merge) and ones with > 3 (to keep separate, not merge)
    valid_groups = group_counts[group_counts['count'] <= 3].drop(columns=['count'])
    invalid_groups = group_counts[group_counts['count'] > 3].drop(columns=['count'])

    # Process valid groups (merge)
    filtered_df = df_to_merge.merge(valid_groups, on=['country_region', 'collection_date', 'Host', 'IsolateName'])
    merged_df = (filtered_df.groupby(['country_region', 'collection_date', 'Host', 'IsolateName'], as_index=False)
             .agg({
                 'Genes': lambda x: ', '.join(sorted(set(g.strip() for v in x for g in v.split(',')))),
                 'Accession': lambda x: ', '.join(sorted(set(x))),  # Combine Accessions
                 'cds': lambda x: ', '.join(sorted(set(x))),  # Combine cds
                 'isolate_source2': lambda x: ', '.join(sorted(set(x))),
                 **{col: 'first' for col in df_to_merge.columns if col not in ['Description', 'record_date', 'organism', 'segment_source',]}  # Keep first value
             })
           )
    # Keep invalid groups as they are
    invalid_df = df_to_merge.merge(invalid_groups, on=['country_region', 'collection_date', 'Host', 'IsolateName'])

    # Combine both datasets
    final_df = pd.concat([merged_df, invalid_df], ignore_index=True)

    # Concatenate merged rows with non-merged rows
    final_merged_df = pd.concat([final_df, df_no_merge], ignore_index=True)
    # print(df.shape[0], final_merged_df.shape[0])
    final_merged_df.to_excel("tmp_gene_after.xlsx")

    combined_num_gene_counts = final_merged_df.loc[final_merged_df['Genes'] != 'NA', 'Genes'].apply(lambda x: len(set(x.split(', '))) if isinstance(x, str) else 0)
    combined_num_gene_distribution = Counter(combined_num_gene_counts)
    combined_num_gene_distribution_f, combined_num_gene_percent_f = format_counts_and_percentages(combined_num_gene_distribution, total=final_merged_df.shape[0])

    combined_genes = Counter()
    for row in final_merged_df['Genes']:
        unique_genes = set(row.split(', '))
        combined_genes.update(unique_genes)
    counts_formatted_combined, percentages_formatted_combined = format_counts_and_percentages(combined_genes, total=final_merged_df.shape[0])

    # # Check merge is correct, isolateName often the same for different isolates - incorrect for Nipah
    # merged_groups = df_to_merge.groupby(['country_region', 'collection_date', 'Host', 'IsolateName']).filter(lambda x: len(x) > 1)

    section.append(f"Counts genes after combining isolate:\n{counts_formatted_combined}\n")
    section.append(f"Percentages genes after combining isolate:\n{percentages_formatted_combined}\n")
    section.append(f"Counts # of genes after combining isolate:\n{combined_num_gene_distribution_f}\n")
    section.append(f"Percentages # of genes after combining isolate:\n{combined_num_gene_percent_f}\n")

    summarize_report.append(section)

    section = ['AlignLens']
    aligns = [int(v['align_len']) for i, v in genes_df.iterrows()]
    section.append(create_binned_seq_lens(aligns))
    summarize_report.append(section)

    section = ['NA length']
    num_na = [int(v['NA_length']) for i, v in genes_df.iterrows()]
    section.append(create_binned_seq_lens(num_na))
    summarize_report.append(section)

    section = ['AA length']
    num_aa = [int(v['AA_length']) for i, v in genes_df.iterrows()]
    section.append(create_binned_seq_lens(num_aa))
    summarize_report.append(section)

    section = ['PcntIDs']
    pcnt_ident = [float(v['pcnt_id']) for i, v in genes_df.iterrows()]
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
