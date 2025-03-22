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
from Utilities import load_csv
from Utilities import dump_csv

from GenBankFunctions import is_reference_genome


def aggregate_references(references, virus_obj, save_data=False):

    # references.to_excel(virus_obj.genbank_raw_ref_file)

    # Aggregation
    grouped_ref = references.groupby(
        ['Authors', 'Title', 'Journal', 'PMID', 'Year'])[
            'accession'].apply(list).reset_index()

    grouped_ref['accession'] = grouped_ref['accession'].apply(
        lambda x: ', '.join(x))
    print("Number of Submission sets following aggregation by exact matches: ",
          len(grouped_ref))

    grouped_ref['RowID'] = grouped_ref.index + 1
    grouped_ref.to_excel(virus_obj.genbank_ref_file)

    # merge rows that are dups
    merged_ref = merge_by_author_title_acc(grouped_ref)

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
        years = [int(y) for y in str(row['Year']).split(',') if y]
        if years:
            year = years[len(years) // 2]
        else:
            year = ''
        if PMID:
            short_name = f"{first_author_surname} ({year}, {PMID})"
        else:

            short_name = f"{first_author_surname} ({year})"

        merged_ref.at[idx, 'ShortName'] = short_name

    print("Number of Submission sets following aggregation by similarity: ",
          len(merged_ref))

    # merged_ref['RefID'] = merged_ref.index + 1
    merged_ref, fixed_ref = get_fixed_Ref_ID(virus_obj, merged_ref)

    if save_data:
        dump_csv(virus_obj.fixed_ref_id_file, fixed_ref)
        merged_ref.to_excel(str(virus_obj.merged_ref_file), index=False)

    return merged_ref


def remove_no_pmid_ref_by_linked_accession(virus, ref):
    ref_with_pmid = ref[ref['PMID'] != '']
    acc_in_ref_with_pmid = set([
        acc.strip()
        for idx, row in ref_with_pmid.iterrows()
        for acc in row['accession'].split(',')
    ])
    # print(acc_in_ref_with_pmid[:3], len(acc_in_ref_with_pmid))
    ref_without_pmid = ref[ref['PMID'] == '']
    # print(len(ref_with_pmid), len(ref_without_pmid))

    remove_ref = []
    keep_ref = []
    for idx, row in ref_without_pmid.iterrows():
        acc_list = [
            acc.strip()
            for acc in row['accession'].split(',')
        ]
        acc_checked = [
            acc
            for acc in acc_list
            if acc in acc_in_ref_with_pmid
        ]

        if (len(acc_checked) / len(acc_list)) == 1:
            remove_ref.append(row)
        else:
            keep_ref.append(row)

    keep_ref = pd.concat([ref_with_pmid, pd.DataFrame(keep_ref)])
    print('Number of Submission sets after remove duplicated Acc:', len(keep_ref))
    pd.DataFrame(remove_ref).to_excel(virus.output_excel_dir / f'{virus.name}_remove_submissions.xlsx')
    return keep_ref


def get_fixed_Ref_ID(virus, references):
    if virus.fixed_ref_id_file.exists():
        fixed_ref = load_csv(virus.fixed_ref_id_file)
    else:
        fixed_ref = []

    def metadata_match(prev_ref, ref):
        for key in ['Title', 'Authors', 'Journal', 'PMID', 'Year', 'accession']:
            if str(prev_ref[key]) != str(ref[key]):
                return False

        return True

    def find_fixed_ref_id(ref):
        for prev_ref in fixed_ref:
            if prev_ref['PMID'].strip() and ref['PMID'].strip() and prev_ref['PMID'].strip() == ref['PMID'].strip():
                return int(prev_ref['RefID'])
            elif metadata_match(prev_ref, ref):
                return int(prev_ref['RefID'])

        return None

    max_ref_id = max([int(r['RefID']) for r in fixed_ref]) if fixed_ref else 0

    for idx, ref in references.iterrows():

        fixed_ref_id = find_fixed_ref_id(ref)
        if not fixed_ref_id:
            max_ref_id += 1
            fixed_ref_id = max_ref_id
            fixed_ref.append({
                'RefID': fixed_ref_id,
                'PMID': ref['PMID'],
                'Title': ref['Title'],
                'Authors': ref['Authors'],
                'Journal': ref['Journal'],
                'Year': ref['Year'],
                'accession': ref['accession']
            })

        references.at[idx, 'RefID'] = fixed_ref_id

    references["RefID"] = references["RefID"].astype(int)

    return references, fixed_ref


def merge_by_author_title_acc(df):
    close_lists = {}

    for _, row_i in df.iterrows():

        close_matches = []
        for _, row_j in df.iterrows():
            if row_i['RowID'] >= row_j['RowID']:
                continue

            score = is_same_submission_set(row_i, row_j)
            if score == 1:
                close_matches.append(row_j['RowID'])

        if len(close_matches) >= 1:
            close_lists[row_i['RowID']] = close_matches

    # for i, j in close_lists.items():
    #     check_list = [367, 369, 464, 475]
    #     if set(check_list) & set([i] + j):
    #         print(i, j)
    list_of_merged_rows, complete_list_of_merged_rows = convert_dict_to_list_of_sets(
        df, close_lists)
    # print(close_lists)
    # for i in list_of_merged_indexes:
    #     if 39 in i:
    #         print(i)
    # print(list_of_merged_indexes, len(complete_list_of_merged_indexes))
    # print("Close lists:", close_lists)
    # print(f'''No with shared author_titles: {len(list_of_merged_rows)}: {list_of_merged_rows}''')
    # print(f'''To be dropped: {len(complete_list_of_shared_indexes)}: {complete_list_of_merged_rows}''')

    list_of_new_rows = []
    for rowIDs in list_of_merged_rows:
        rowIDs = list(rowIDs)

        new_row = merge_rows(df, rowIDs)
        list_of_new_rows.append(new_row)

        # rows = df[df['RowID'].isin(rowIDs)]
        # rows_with_pmid = [
        #     i
        #     for i, row in rows.iterrows()
        #     if row['PMID'].strip()
        # ]

        # if len(rows_with_pmid) < 2:
        #     new_row = merge_rows(df, rowIDs)
        #     list_of_new_rows.append(new_row)
        # else:
        #     first_rowID = rows_with_pmid[0]
        #     new_rowIDs = [first_rowID] + [
        #         _i
        #         for _i in rowIDs
        #         if _i not in rows_with_pmid
        #     ]
        #     new_row = merge_rows(df, new_rowIDs)
        #     list_of_new_rows.append(new_row)

        #     list_of_new_rows.append(
        #         df[df['RowID'].isin(rows_with_pmid[1:])]
        #     )

    list_of_new_rows.insert(
        0,
        df[~df['RowID'].isin(complete_list_of_merged_rows)])
    list_of_new_rows = pd.concat(list_of_new_rows, ignore_index=True)
    return list_of_new_rows


def is_same_submission_set(row_i, row_j):
    # !!Change the order of comparison will change the final result

    # PMID
    if row_i['PMID'] and row_j['PMID']:
        if str(row_i['PMID']).strip() == str(row_j['PMID']).strip():
            return 1
        else:
            return 0

    accessions_i = [
        r.strip()
        for r in row_i['accession'].split(',')
        if r.strip()
    ]
    accessions_i = [s for s in accessions_i if not is_reference_genome(s)]

    accessions_j = [
        r.strip()
        for r in row_j['accession'].split(',')
        if r.strip()
    ]
    accessions_j = [s for s in accessions_j if not is_reference_genome(s)]

    author_score = closed_author(row_i, row_j)
    title_score = closed_title(row_i, row_j)
    journal_score = closed_journal(row_i, row_j)
    acc_score = (
        closed_accession_stem(accessions_i, accessions_j)
        or
        closed_accession_group(accessions_i, accessions_j)
    )
    year_score = closed_year(row_i, row_j)

    valid_column = get_valid_fuzzy_column(row_i) & get_valid_fuzzy_column(row_j)

    if 'Title' in valid_column:
        multiple_scores = {
            ('Authors', 'Title', 'Year'): author_score and title_score and year_score,
            ('Authors', 'Title', 'Accession'): author_score and title_score and acc_score,
            ('Title', 'Journal', 'Year'): title_score and journal_score and year_score,
            ('Title', 'Journal', 'Accession'): title_score and journal_score and acc_score,
            ('Authors', 'Title', 'Journal'): author_score and title_score and journal_score,
            # journal year, accession
        }
    else:
        multiple_scores = {
            ('Authors', 'Journal', 'Year'): author_score and journal_score and year_score,
            ('Authors', 'Year', 'Accession'): author_score and year_score and acc_score,
            ('Authors', 'Journal', 'Accession'): author_score and journal_score and acc_score,
            ('Authors', 'Accession'): author_score and closed_accession_group(accessions_i, accessions_j, 1) and (len(accessions_i) >= 10)
            # journal year, accession
        }

    title_list = [
        'Emerging viruses are an underestimated cause of undiagnosed febrile illness in uganda'
    ]
    if row_i['Title'] in title_list and row_j['Title'] in title_list:
        print(multiple_scores)
        print(valid_column)

    multiple_scores = [
        v
        for k, v in multiple_scores.items()
        if set(k).issubset(valid_column)
    ]
    if any(multiple_scores):
        return 1

    return 0


def get_valid_fuzzy_column(row):
    columns = []
    if row['Authors'] and row['Authors'] != 'NCBI':
        columns.append('Authors')
    if row['Title'] and row['Title'] != 'Direct Submission':
        columns.append('Title')
    if row['Journal']:
        columns.append('Journal')

    columns.append('Year')
    columns.append('Accession')

    return set(columns)


def closed_accession_stem(accessions_i, accessions_j):
    if not accessions_i or not accessions_j:
        return 0

    pcnt_shared_stems = get_pcnt_shared_stems(accessions_i, accessions_j, 3)
    if pcnt_shared_stems > 0.75:
        return 1
    else:
        return 0


def closed_accession_group(accessions_i, accessions_j, threshold=0.8):
    # Accession
    if accessions_i and accessions_j:
        if (set(accessions_i).issubset(set(accessions_j)) or set(accessions_j).issubset(set(accessions_i))):
            return 1

    # idealy should be sharing 100%
    pcnt_shared_accessions = get_pcnt_shared_accessions(accessions_i, accessions_j)
    if pcnt_shared_accessions >= threshold:
        return 1

    return 0


def closed_author(row_i, row_j):

    authors_i = row_i['Authors'].upper()
    authors_j = row_j['Authors'].upper()

    if not authors_i or not authors_j:
        return 0

    if (authors_i == 'NCBI') or (authors_j == 'NCBI'):
        return 0

    pcnt_authors_overlap = get_pcnt_authors_overlap(authors_i, authors_j)
    if pcnt_authors_overlap > 0.66:
        return 1
    else:
        return 0


def closed_title(row_i, row_j):
    title_i = row_i['Title'].upper()
    title_j = row_j['Title'].upper()

    if not title_i or not title_j:
        return 0

    if (title_i == 'Direct Submission'.upper()) or (title_j == 'Direct Submission'.upper()):
        return 0

    if title_i in title_j:
        return 1
    elif title_j in title_i:
        return 1

    title_distance = Levenshtein.distance(title_i, title_j)
    if title_distance < 5:
        return 1
    else:
        return 0


def closed_year(row_i, row_j):
    year_i = str(row_i['Year'])
    year_j = str(row_j['Year'])

    if not year_i and not year_j:
        return 1
    elif year_i and not year_j:
        return 1
    elif year_j and not year_i:
        return 1

    max_year_dif = calc_year_dif(year_i, year_j)
    if max_year_dif <= 1:
        return 1
    else:
        return 0


def closed_journal(row_i, row_j):
    journal_i = row_i['Journal']
    journal_j = row_j['Journal']

    if not journal_i or not journal_j:
        return 0

    # if (journal_i.lower() == 'Unpublished'.lower()) or (journal_j.lower() == 'Unpublished'.lower()):
    #     return 0

    # if journal_i.lower() in journal_i.lower():
    #     return 1
    # elif journal_j.lower() in journal_i.lower():
    #     return 1

    distance = Levenshtein.distance(journal_i, journal_j)
    if distance < 5:
        return 1
    else:
        return 0


def merge_similar_title(title_list):
    new_title_list = []
    for idx, t1 in enumerate(title_list):
        for jdx, t2 in enumerate(title_list):
            if idx >= jdx:
                continue
            title_distance = Levenshtein.distance(t1, t2)
            if title_distance < 5:
                if t1 in new_title_list or t2 in new_title_list:
                    continue
                else:
                    new_title_list.append(t1)

    return new_title_list


def merge_rows(df, merged_rowID):
    new_row = {}
    rows = df[df['RowID'].isin(merged_rowID)]

    authors_list = rows['Authors'].tolist()
    new_row['Authors'] = combine_items_in_different_lists(authors_list, spliter=',')

    titles_list = rows['Title'].tolist()
    titles_list = [i.capitalize() for i in titles_list if i != 'Direct Submission']
    if titles_list:
        # titles_list = merge_similar_title(title_list)
        new_row['Title'] = combine_items_in_different_lists(titles_list)
    else:
        new_row['Title'] = 'Direct Submission'

    journal_list = rows['Journal'].tolist()
    new_row['Journal'] = combine_items_in_different_lists(journal_list)

    pmid_list = rows['PMID'].tolist()
    new_row['PMID'] = combine_items_in_different_lists(pmid_list)

    year_list = rows['Year'].tolist()
    new_row['Year'] = combine_items_in_different_lists(year_list)

    accession_list = rows['accession'].tolist()
    new_row['accession'] = combine_items_in_different_lists(accession_list, spliter=',')

    new_row['merged_indexes'] = ';'.join([str(i) for i in sorted(merged_rowID)])
    new_row['num_merged_indexes'] = len(merged_rowID)
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
