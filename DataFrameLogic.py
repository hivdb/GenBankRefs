import pandas as pd
import re
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


def aggregate_references(references, virus_obj):

    # references.to_excel(virus_obj.genbank_raw_ref_file)

    # Aggregation
    grouped_ref = references.groupby(
        ['Authors', 'Title', 'Journal', 'PMID', 'Year'])[
            'accession'].apply(list).reset_index()

    grouped_ref['accession'] = grouped_ref['accession'].apply(
        lambda x: ', '.join(x))
    print("Number of entries following aggregation by exact matches: ",
          len(grouped_ref))

    grouped_ref.to_excel(virus_obj.genbank_ref_file)

    # merge rows that are dups
    merged_ref = merge_by_author_title_acc(grouped_ref)
    merged_ref['RefID'] = merged_ref.index + 1

    for idx, row in merged_ref.iterrows():
        authors = row['Authors']
        if ',' in authors:
            authors = authors.split(',')
        else:
            authors = authors.split(';')
        first_author_surname = ''
        if authors:
            first_author = authors[0]
            first_author_name_list = first_author.split()
            if len(first_author_name_list) == 1:
                first_author_surname = first_author_name_list[0]
            else:
                first_author_surname = ' '.join(first_author_name_list[:-1])
        merged_ref.at[idx, 'FirstAuthorSurname'] = first_author_surname

        PMID = row['PMID']
        years = [int(y) for y in str(row['Year']).split(',')]
        year = years[len(years) // 2]
        if PMID:
            short_name = f"{first_author_surname} ({year}, {PMID})"
        else:

            short_name = f"{first_author_surname} ({year})"

        merged_ref.at[idx, 'ShortName'] = short_name

    print("Number of entries following aggregation by similarity: ",
          len(merged_ref))

    merged_ref.to_excel(str(virus_obj.merged_ref_file), index=False)

    return merged_ref


def merge_by_author_title_acc(df):
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
    # print(close_lists)
    # for i in list_of_merged_indexes:
    #     if 39 in i:
    #         print(i)
    # print(list_of_merged_indexes, len(complete_list_of_merged_indexes))
    # print("Close lists:", close_lists)
    # print(f'''No with shared author_titles: {len(list_of_sets_w_shared_indexes)}: {list_of_sets_w_shared_indexes}''')
    # print(f'''To be dropped: {len(complete_list_of_shared_indexes)}: {complete_list_of_shared_indexes}''')

    list_of_new_rows = []
    for indexes in list_of_merged_indexes:
        indexes = list(indexes)
        rows = df.loc[indexes]

        indexes_with_pmid = [
            i
            for i, row in rows.iterrows()
            if row['PMID'].strip()
        ]

        if len(indexes_with_pmid) < 2:
            new_row = merge_rows(df, indexes)
            list_of_new_rows.append(new_row)
        else:
            first_index = indexes_with_pmid[0]
            new_indexes = [first_index] + [
                _i
                for _i in indexes
                if _i not in indexes_with_pmid
            ]
            new_row = merge_rows(df, new_indexes)
            list_of_new_rows.append(new_row)

            list_of_new_rows.append(df.loc[indexes_with_pmid[1:]])

    df = df.drop(complete_list_of_merged_indexes)
    list_of_new_rows.insert(0, df)
    list_of_new_rows = pd.concat(list_of_new_rows, ignore_index=True)
    return list_of_new_rows


def compare_authors_titles_year_accession_overlap(row_i, row_j):
    # !!Change the order of comparison will change the final result

    authors_i = row_i['Authors']
    authors_j = row_j['Authors']
    title_i = row_i['Title']
    title_j = row_j['Title']
    year_i = str(row_i['Year'])
    year_j = str(row_j['Year'])

    accessions_i = [
        r.strip()
        for r in row_i['accession'].split(',')
        if r.strip()
    ]
    accessions_j = [
        r.strip()
        for r in row_j['accession'].split(',')
        if r.strip()
    ]

    accessions_i = [s for s in accessions_i if not is_reference_genome(s)]
    accessions_j = [s for s in accessions_j if not is_reference_genome(s)]

    # def log_match_reference():
    #     with open('MatchingReferences.txt', 'a') as file:
    #         title_distance = Levenshtein.distance(title_i, title_j)
    #         max_year_dif = calc_year_dif(year_i, year_j)
    #         pcnt_authors_overlap = get_pcnt_authors_overlap(authors_i, authors_j)
    #         pcnt_shared_accessions = get_pcnt_shared_accessions(accessions_i, accessions_j)
    #         pcnt_shared_stems = get_pcnt_shared_stems(accessions_i, accessions_j, 3)

    #         file.write(f'Title_i:{title_i}\nTitle_j:{title_j}\nTitle_distance:{title_distance}\n')
    #         file.write(f'Authors_i: {authors_i}\nAuthors_j: {authors_j}\n')
    #         file.write(f'Authors_overlap:{pcnt_authors_overlap}\n')
    #         file.write(f'Year_i:{year_i} Year_j:{year_j} Max_year_dif:{max_year_dif}\n')
    #         file.write(f'Accessions_i:{accessions_i}\nAccessions_j:{accessions_j}\n')
    #         file.write(f'Pcnt_shared_accessions:{pcnt_shared_accessions} Pcnt_shared_stems:{pcnt_shared_stems}\n\n')

    # PMID
    if row_i['PMID'] and row_j['PMID']:
        if row_i['PMID'] == row_j['PMID']:
            return 1
        else:
            return 0

    # Accession
    if accessions_i and accessions_j:
        if (set(accessions_i).issubset(set(accessions_j)) or set(accessions_j).issubset(set(accessions_i))):
            return 1

    pcnt_shared_accessions = get_pcnt_shared_accessions(accessions_i, accessions_j)
    if pcnt_shared_accessions > 0.8:
        return 1

    # Title - when both not Direct Submission
    title_distance = Levenshtein.distance(title_i, title_j)
    if title_i != 'Direct Submission' and (title_j != 'Direct Submission') and title_distance < 5:
        return 1

    if (title_i != 'Direct Submission') and (title_j != 'Direct Submission'):
        return 0

    # if any author not NCBI (empty), skip comparison
    if (authors_i and authors_i != 'NCBI') or (authors_j and authors_j != 'NCBI'):
        return 0

    # Title, author, accession stem
    max_year_dif = calc_year_dif(year_i, year_j)
    pcnt_authors_overlap = get_pcnt_authors_overlap(authors_i, authors_j)
    pcnt_shared_stems = get_pcnt_shared_stems(accessions_i, accessions_j, 3)
    if pcnt_authors_overlap >= 0.75 \
            and max_year_dif <= 1 \
            and pcnt_shared_stems > 0.75:
        return 1
    return 0


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
    return pd.DataFrame([new_row])


def merge_feature_rows(df, genes_df):
    df = df.copy()
    df = df.replace("", "NA")
    new_row = {}

    unique_organisms = count_unique_elements(df['organism'].tolist())
    new_row['Organisms'] = dict_to_sorted_string(unique_organisms)

    unique_record_years = count_unique_elements(df['RecordYear'].tolist())
    new_row['RecordYears'] = dict_to_sorted_string(unique_record_years)

    unique_hosts = count_unique_elements(df['Host'].tolist())
    new_row['Hosts'] = dict_to_sorted_string(unique_hosts)

    unique_specimens = count_unique_elements(df['isolate_source'].tolist())
    new_row['Specimens'] = dict_to_sorted_string(unique_specimens)

    unique_gene = count_unique_elements(df['Genes'].tolist())
    new_row['Genes'] = dict_to_sorted_string(unique_gene)

    unique_countries = count_unique_elements(df['Country'].tolist())
    new_row['Countries'] = dict_to_sorted_string(unique_countries)

    unique_isolate_years = count_unique_elements(df['IsolateYear'].tolist())
    new_row['IsolateYears'] = dict_to_sorted_string(unique_isolate_years)

    unique_specimens = count_unique_elements(df['isolate_source'].tolist())
    new_row['Specimens'] = dict_to_sorted_string(unique_specimens)

    unique_gene = count_unique_elements(df['Genes'].tolist())
    new_row['Gene'] = dict_to_sorted_string(unique_gene)

    new_row['NumSubSeqs'] = create_binned_seq_lens(df['NumSubSeqs'].tolist())

    new_row['NumNA'] = create_binned_seq_lens(genes_df['NA_length'].tolist())
    new_row['NumAA'] = create_binned_seq_lens(genes_df['AA_length'].tolist())

    new_row['AlignLens'] = create_binned_seq_lens(genes_df['align_len'].tolist())
    new_row['PcntIDs'] = create_binned_pcnts(genes_df['pcnt_id'].tolist())
    return new_row


def combine_refs_and_features(ref_df, features_df, genes_df):
    combined_df = ref_df.copy()
    feature_columns = [
        'Organisms', 'RecordYears',  'Hosts', 'Countries', 'Gene',
        'IsolateYears', 'Specimens', 'NumNA', 'NumAA',
        'AlignLens', 'PcntIDs']
    combined_df[feature_columns] = 'None'
    count = 0
    for index, row in combined_df.iterrows():
        count += 1
        accession_string = row['accession']
        accession_list = accession_string.split(',')
        accession_list = [i.strip() for i in accession_list]
        features_rows = features_df[features_df['Accession'].isin(
            accession_list)]

        genes_row = genes_df[genes_df['Accession'].isin(accession_list)]

        new_dict = merge_feature_rows(features_rows, genes_row)

        combined_df.at[index, 'Organisms'] = new_dict['Organisms']
        combined_df.at[index, 'RecordYears'] = new_dict['RecordYears']
        combined_df.at[index, 'Hosts'] = new_dict['Hosts']
        combined_df.at[index, 'Countries'] = new_dict['Countries']
        combined_df.at[index, 'IsolateYears'] = new_dict['IsolateYears']
        combined_df.at[index, 'Specimens'] = new_dict['Specimens']
        combined_df.at[index, 'Gene'] = new_dict['Gene']
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
