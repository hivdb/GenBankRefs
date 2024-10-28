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
                       extract_year_from_date_fields,
                       create_binned_seq_lens,
                       create_binned_pcnts)

from GenBankFunctions import is_reference_genome

def compare_authors_titles(row_i, row_j):
    authors_i = row_i['authors']
    authors_j = row_j['authors']
    title_i = row_i['title']
    title_j = row_j['title']
    year_i = str(row_i['year'])
    year_j = str(row_j['year'])
    accessions_i = row_i['accession']
    accessions_j = row_j['accession']

    pcnt_authors_overlap = get_pcnt_authors_overlap(authors_i, authors_j)
    title_distance = Levenshtein.distance(title_i, title_j)
    max_year_dif = calc_year_dif(year_i, year_j)
    pcnt_shared_accessions = get_pcnt_shared_accessions(accessions_i, accessions_j)
    pcnt_shared_stems = get_pcnt_shared_stems(accessions_i, accessions_j)

    if title_i != 'Direct Submission' and title_distance < 5:
        match = 1
    elif (title_i == 'Direct Submission') | (title_j == 'Direct Submission') \
        and authors_i != 'NCBI' \
        and pcnt_authors_overlap >= 0.75 \
        and max_year_dif <= 1 \
        and pcnt_shared_stems >0.75:
            match = 1
    else:
        match = 0

    if match == 1:
        with open('MatchingReferences.txt', 'a') as file:
            file.write(f'Title_i:{title_i}\nTitle_j:{title_j}\nTitle_distance:{title_distance}\n')
            file.write(f'Authors_i: {authors_i}\nAuthors_j: {authors_j}\n')
            file.write(f'Authors_overlap:{pcnt_authors_overlap}\n')
            file.write(f'Year_i:{year_i} Year_j:{year_j} Max_year_dif:{max_year_dif}\n')
            file.write(f'Accessions_i:{accessions_i}\nAccessions_j:{accessions_j}\n')
            file.write(f'Pcnt_shared_accessions:{pcnt_shared_accessions} Pcnt_shared_stems:{pcnt_shared_stems}\n\n')
    return match


def process_authors_titles(df):
    close_lists = {}
    for i, row_i in df.iterrows():
        close_matches = []
        for j, row_j in df.iterrows():
            if i >= j:
                continue
            score = compare_authors_titles(row_i, row_j)
            if score == 1:
                close_matches.append(j)
            if len(close_matches) >=1:
                close_lists[i] = close_matches

    list_of_sets_w_shared_indexes, complete_list_of_shared_indexes = convert_dict_to_list_of_sets(close_lists)
    #print("Close lists:", close_lists)
    #print(f'''No with shared author_titles: {len(list_of_sets_w_shared_indexes)}: {list_of_sets_w_shared_indexes}''')
    #print(f'''To be dropped: {len(complete_list_of_shared_indexes)}: {complete_list_of_shared_indexes}''')

    list_of_new_rows = []
    for item in list_of_sets_w_shared_indexes:
        new_row = merge_rows(df, list(item))
        list_of_new_rows.append(new_row)

    df = df.drop(complete_list_of_shared_indexes)
    list_of_new_rows.insert(0, df)
    list_of_new_rows = pd.concat(list_of_new_rows, ignore_index=True)
    return list_of_new_rows


def process_accession_lists(df):
    accession_list = df['accession']
    close_lists = {}
    for i, item1 in enumerate(accession_list):
        close_matches = []
        for j, item2 in enumerate(accession_list):
            if i >= j or is_reference_genome(item1[0]) == True:
                continue
            pcnt_shared_accessions = get_pcnt_shared_accessions(item1, item2)
            if pcnt_shared_accessions >= 0.9:
                close_matches.append(j)
        if len(close_matches) >=1:
            close_lists[i] = close_matches

    (list_of_sets_w_shared_indexes, complete_list_of_shared_indexes) = convert_dict_to_list_of_sets(close_lists)
    #print(f'''No with shared accessions: {len(list_of_sets_w_shared_indexes)}: {list_of_sets_w_shared_indexes}''')
    #print(f'''To be dropped: {len(complete_list_of_shared_indexes)}: {complete_list_of_shared_indexes}''')

    list_of_new_rows = []
    for item in list_of_sets_w_shared_indexes:
        new_row = merge_rows(df, list(item))
        list_of_new_rows.append(new_row)

    df = df.drop(complete_list_of_shared_indexes)
    list_of_new_rows.insert(0, df)
    list_of_new_rows = pd.concat(list_of_new_rows, ignore_index=True)
    return list_of_new_rows


def merge_rows(df, shared_indexes):
    new_row = {}
    authors_list = df.loc[shared_indexes, 'authors'].tolist()
    titles_list = df.loc[shared_indexes, 'title'].tolist()
    journal_list = df.loc[shared_indexes, 'journal'].tolist()
    pmid_list = df.loc[shared_indexes, 'pmid'].tolist()
    year_list = df.loc[shared_indexes, 'year'].tolist()
    accession_list = df.loc[shared_indexes, 'accession'].tolist()
    new_accessions = combine_items_in_different_lists(accession_list)
    new_years = combine_items_in_different_lists(year_list)
    new_titles = combine_items_in_different_lists(titles_list)
    new_authors = combine_items_in_different_lists(authors_list)
    new_pmids = combine_items_in_different_lists(pmid_list)
    new_journals = combine_items_in_different_lists(journal_list)
    new_row['authors'] = new_authors
    new_row['year'] = new_years
    new_row['title'] = new_titles
    new_row['pmid'] = new_pmids
    new_row['journal'] = new_journals
    new_row['accession'] = new_accessions
    new_row_df = pd.DataFrame(new_row, index = [0])
    return new_row_df


def merge_feature_rows(df):
    df = df.copy()
    df = df.replace("", "NA")
    new_row = {}
    unique_organisms = count_unique_elements(df['organism'].tolist())
    new_row['Organisms'] = dict_to_sorted_string(unique_organisms)
    #df['record_year'] = df['record_date'].astype(str).str[:4]
    #df['record_year'] = extract_year_from_date_fields(df['record_date'].astype(str))
    df['record_year'] = df['record_date'].apply(extract_year_from_date_fields)
    #df['record_year'] = df['record_date']
    unique_record_years = count_unique_elements(df['record_year'].tolist())
    new_row['RecordYears'] =  dict_to_sorted_string(unique_record_years)
    unique_hosts = count_unique_elements(df['host'].tolist())
    new_row['Hosts'] =  dict_to_sorted_string(unique_hosts)
    unique_countries = count_unique_elements(df['country_region'].tolist())
    new_row['Countries'] = dict_to_sorted_string(unique_countries)
    #df['isolate_year'] = df['collection_date'].astype(str).str[:4]
    #df['isolate_year'] = extract_year_from_date_fields(df['collection_date'].astype(str))
    df['isolate_year'] = df['collection_date'].apply(extract_year_from_date_fields)
    #df['isolate_year'] = df['collection_date']
    unique_isolate_years = count_unique_elements(df['isolate_year'].tolist())
    new_row['IsolateYears'] = dict_to_sorted_string(unique_isolate_years)
    unique_cds = count_unique_elements(df['cds'].tolist())
    new_row['CDS'] = dict_to_sorted_string(unique_cds)
    new_row['SeqLens'] = create_binned_seq_lens(df['seq_len'].tolist())
    new_row['AlignLens'] = create_binned_seq_lens(df['align_len'].tolist())
    new_row['PcntIDs'] = create_binned_pcnts(df['pcnt_id'].tolist())
    return new_row


def combine_refs_and_features(ref_df, features_df):
    combined_df = ref_df.copy()
    feature_columns = ['Organisms', 'RecordYears',  'Hosts', 'Countries', 
                      'IsolateYears', 'CDS', 'SeqLens', 'AlignLens', 'PcntIDs']
    combined_df[feature_columns] = 'None'
    
    count = 0
    for index, row in combined_df.iterrows():
        count += 1
        accession_string = row['accession']
        accession_list = accession_string.split(', ')
        print("\n", index)
        print(accession_list)
        features_rows = features_df[features_df['acc_num'].isin(accession_list)]
        print("features_rows_dates: ", features_rows['record_date'], " ", features_rows['collection_date'], "\n")
        new_dict = merge_feature_rows(features_rows)
        print(new_dict)
        combined_df.at[index, 'Organisms'] = new_dict['Organisms']
        combined_df.at[index, 'RecordYears'] = new_dict['RecordYears']
        combined_df.at[index, 'Hosts'] = new_dict['Hosts']
        combined_df.at[index, 'Countries'] = new_dict['Countries']
        combined_df.at[index, 'IsolateYears'] = new_dict['IsolateYears']
        combined_df.at[index, 'CDS'] = new_dict['CDS']
        combined_df.at[index, 'SeqLens'] = new_dict['SeqLens']
        combined_df.at[index, 'AlignLens'] = new_dict['AlignLens']
        combined_df.at[index, 'PcntIDs'] = new_dict['PcntIDs']     

    return combined_df


