import pandas as pd
import numpy as np
import re
import os
from collections import Counter
import Levenshtein

from Utilities import (get_pcnt_authors_overlap,
                       get_pcnt_shared_accessions,
                       get_pcnt_shared_stems,
                       calc_year_dif,
                       count_unique_elements,
                       dict_to_sorted_string,
                       convert_dict_to_list_of_sets,
                       combine_items_in_different_lists,
                       create_binned_seq_lens,
                       create_binned_pcnts)

from GenBankFunctions import is_reference_genome


def compare_authors_titles_year_accession_overlap(row_i, row_j):
    # !!Change the order of comparison will change the final result

    authors_i = row_i['Authors']
    authors_j = row_j['Authors']
    title_i = row_i['Title']
    title_j = row_j['Title']
    year_i = str(row_i['Year'])
    year_j = str(row_j['Year'])
    accessions_i = row_i['accession']
    accessions_j = row_j['accession']

    def log_match_reference():
        with open('MatchingReferences.txt', 'a') as file:
            title_distance = Levenshtein.distance(title_i, title_j)
            max_year_dif = calc_year_dif(year_i, year_j)
            pcnt_authors_overlap = get_pcnt_authors_overlap(authors_i, authors_j)
            pcnt_shared_accessions = get_pcnt_shared_accessions(accessions_i, accessions_j)
            pcnt_shared_stems = get_pcnt_shared_stems(accessions_i, accessions_j, 3)

            file.write(f'Title_i:{title_i}\nTitle_j:{title_j}\nTitle_distance:{title_distance}\n')
            file.write(f'Authors_i: {authors_i}\nAuthors_j: {authors_j}\n')
            file.write(f'Authors_overlap:{pcnt_authors_overlap}\n')
            file.write(f'Year_i:{year_i} Year_j:{year_j} Max_year_dif:{max_year_dif}\n')
            file.write(f'Accessions_i:{accessions_i}\nAccessions_j:{accessions_j}\n')
            file.write(f'Pcnt_shared_accessions:{pcnt_shared_accessions} Pcnt_shared_stems:{pcnt_shared_stems}\n\n')

    if (accessions_i == accessions_j) and is_reference_genome(accessions_i):
        return 0

    # Title
    title_distance = Levenshtein.distance(title_i, title_j)
    if title_i != 'Direct Submission' and title_distance < 5:
        return 1

    # Ref seq
    if is_reference_genome(accessions_i[0]) or is_reference_genome(accessions_j[0]):
        return 0

    # if row_i['PMID'] and row_j['PMID'] and row_i['PMID'] == row_j['PMID']:
    #     match = 1
    #     log_match_reference()
    #     return match

    # Accesion only
    pcnt_shared_accessions = get_pcnt_shared_accessions(accessions_i, accessions_j)
    if pcnt_shared_accessions > 0.8:
        return 1

    # Title, author, accession stem
    max_year_dif = calc_year_dif(year_i, year_j)
    pcnt_authors_overlap = get_pcnt_authors_overlap(authors_i, authors_j)
    pcnt_shared_stems = get_pcnt_shared_stems(accessions_i, accessions_j, 3)
    if (title_i == 'Direct Submission') | (title_j == 'Direct Submission') \
            and authors_i != 'NCBI' \
            and pcnt_authors_overlap >= 0.75 \
            and max_year_dif <= 1 \
            and pcnt_shared_stems > 0.75:
        return 1

    return 0


def process_authors_titles(df):
    close_lists = {}
    for i, row_i in df.iterrows():

        close_matches = []
        for j, row_j in df.iterrows():
            if i >= j:
                continue

            score = compare_authors_titles_year_accession_overlap(row_i, row_j)
            if score == 1:
                close_matches.append(j)

        if len(close_matches) >= 1:
            close_lists[i] = close_matches

    list_of_merged_indexes, complete_list_of_merged_indexes = convert_dict_to_list_of_sets(close_lists)
    # print(len(list_of_merged_indexes), len(complete_list_of_merged_indexes))
    # print("Close lists:", close_lists)
    # print(f'''No with shared author_titles: {len(list_of_sets_w_shared_indexes)}: {list_of_sets_w_shared_indexes}''')
    # print(f'''To be dropped: {len(complete_list_of_shared_indexes)}: {complete_list_of_shared_indexes}''')

    list_of_new_rows = []
    for item in list_of_merged_indexes:
        new_row = merge_rows(df, list(item))
        list_of_new_rows.append(new_row)

    df = df.drop(complete_list_of_merged_indexes)
    list_of_new_rows.insert(0, df)
    list_of_new_rows = pd.concat(list_of_new_rows, ignore_index=True)
    return list_of_new_rows


def merge_rows(df, merged_indexes):
    new_row = {}
    authors_list = df.loc[merged_indexes, 'Authors'].tolist()
    titles_list = df.loc[merged_indexes, 'Title'].tolist()
    titles_list = [i for i in titles_list if i != 'Direct Submission']
    journal_list = df.loc[merged_indexes, 'Journal'].tolist()
    pmid_list = df.loc[merged_indexes, 'PMID'].tolist()
    year_list = df.loc[merged_indexes, 'Year'].tolist()
    accession_list = df.loc[merged_indexes, 'accession'].tolist()
    new_row['Authors'] = combine_items_in_different_lists(authors_list)
    new_row['Year'] = combine_items_in_different_lists(year_list)
    new_row['Title'] = combine_items_in_different_lists(titles_list)
    new_row['PMID'] = combine_items_in_different_lists(pmid_list)
    new_row['Journal'] = combine_items_in_different_lists(journal_list)
    new_row['accession'] = combine_items_in_different_lists(accession_list)
    # new_row['merged_indexes'] = ';'.join([str(i) for i in merged_indexes])
    new_row_df = pd.DataFrame(new_row, index=[0])
    return new_row_df


def merge_feature_rows(df):
    df = df.copy()
    df = df.replace("", "NA")
    new_row = {}

    unique_organisms = count_unique_elements(df['organism'].tolist())
    new_row['Organisms'] = dict_to_sorted_string(unique_organisms)

    unique_record_years = count_unique_elements(df['RecordYear'].tolist())
    new_row['RecordYears'] = dict_to_sorted_string(unique_record_years)

    unique_hosts = count_unique_elements(df['Host'].tolist())
    new_row['Hosts'] = dict_to_sorted_string(unique_hosts)

    unique_countries = count_unique_elements(df['Country'].tolist())
    new_row['Countries'] = dict_to_sorted_string(unique_countries)

    unique_isolate_years = count_unique_elements(df['IsolateYear'].tolist())
    new_row['IsolateYears'] = dict_to_sorted_string(unique_isolate_years)

    unique_specimens = count_unique_elements(df['isolate_source'].tolist())
    new_row['Specimens'] = dict_to_sorted_string(unique_specimens)

    unique_gene = count_unique_elements(df['segment_source'].tolist())
    new_row['SegmentSource'] = dict_to_sorted_string(unique_gene)

    unique_gene = count_unique_elements(df['Genes'].tolist())
    new_row['Gene'] = dict_to_sorted_string(unique_gene)

    unique_cds = count_unique_elements(df['cds'].tolist())
    new_row['CDS'] = dict_to_sorted_string(unique_cds)

    new_row['NumNA'] = create_binned_seq_lens(df['NumNA'].tolist())
    new_row['NumAA'] = create_binned_seq_lens(df['NumAA'].tolist())
    new_row['AlignLens'] = create_binned_seq_lens(df['align_len'].tolist())
    new_row['PcntIDs'] = create_binned_pcnts(df['pcnt_id'].tolist())
    return new_row


def combine_refs_and_features(ref_df, features_df):
    combined_df = ref_df.copy()
    feature_columns = ['Organisms', 'RecordYears',  'Hosts', 'Countries', 'Gene',
                       'IsolateYears', 'Specimens', 'CDS', 'NumNA', 'NumAA', 'AlignLens', 'PcntIDs']
    combined_df[feature_columns] = 'None'
    count = 0
    for index, row in combined_df.iterrows():
        count += 1
        accession_string = row['accession']

        accession_list = accession_string.split(',')
        accession_list = [i.strip() for i in accession_list]
        features_rows = features_df[features_df['Accession'].isin(
            accession_list)]

        new_dict = merge_feature_rows(features_rows)

        combined_df.at[index, 'Organisms'] = new_dict['Organisms']
        combined_df.at[index, 'RecordYears'] = new_dict['RecordYears']
        combined_df.at[index, 'Hosts'] = new_dict['Hosts']
        combined_df.at[index, 'Countries'] = new_dict['Countries']
        combined_df.at[index, 'IsolateYears'] = new_dict['IsolateYears']
        combined_df.at[index, 'Specimens'] = new_dict['Specimens']
        combined_df.at[index, 'Gene'] = new_dict['Gene']
        combined_df.at[index, 'CDS'] = new_dict['CDS']
        combined_df.at[index, 'NumNA'] = new_dict['NumNA']
        combined_df.at[index, 'NumAA'] = new_dict['NumAA']
        combined_df.at[index, 'AlignLens'] = new_dict['AlignLens']
        combined_df.at[index, 'PcntIDs'] = new_dict['PcntIDs']
        # remove .x after accession ID
        row["accession"] = re.sub(r'\.\d+', '', row['accession'])
    return combined_df


def compare_output_files(saved_df, new_df):
    saved_df = saved_df.sort_values(by='Authors')
    new_df = new_df.sort_values(by='Authors')

    print(f'Number of rows: Saved file:{len(saved_df)} New file:{len(new_df)}')
    if (saved_df.columns == new_df.columns).all():
        print("The DataFrames have the same columns in the same order.")
    else:
        print("The DataFrames do not have the same columns or order.")

    for (index_i, row_i), (index_j, row_j) in zip(saved_df.fillna('').iterrows(), new_df.fillna('').iterrows()):
        for col in saved_df.columns:
            if row_i[col] == row_j[col] or (pd.isna(row_i[col]) and pd.isna(row_j[col])):
                continue
            else:
                print(f'Column:{col}: {row_i[col]}\n{row_j[col]}\n')
                input('pause')
