from Bio import SeqIO
from Bio import Entrez
import pandas as pd
import re
import Levenshtein

Entrez.email = "rshafer.stanford.edu"

def extract_year_from_journal(text):
    match = re.search(r'\((\d{4})\)', text)
    if match:
        return match.group(1)
    match = re.search (r'\d{2}-[A-Z]{3}-(\d{4})', text)
    if match:
        return match.group(1)
    return ''

# def is_incorrect_accession(accession):
#     if accession == 'NM_010185.4':
#         return True
#     else: 
#         return False

def is_reference_genome(acc):
    prefix_list = ['NC', 'NG', 'NM', 'NR']
    for item in prefix_list:
        if acc.startswith(item):
            return True
    return False


def process_author_field(names):
    if names == '':
        return 'NCBI'
    names = names.replace(" and", ",")
    name_list = names.split(", ")
    processed_names = []
    for name in name_list:
        parts = name.split(",")
        if len(parts) == 2:
            last_name = parts[0]  # The last name is the first part
            last_name = last_name.capitalize()
            initials = parts[1]  # The initials are the second part
            first_initial = initials[0]  # Take only the first character of the initials
            processed_name = last_name + ' ' + first_initial + '.'
            processed_names.append(processed_name)  
    processed_names = ', '.join(processed_names)
    return (processed_names)


#Convert {key: [values], key: [values]] to a list of sets comprising {key, values}
#In this program, each key is the index of a row and the values contain
# one or more indexes of rows that share some property
#The subset flag code ensures sets do not share items with one another
def convert_dict_to_list_of_sets(dict):
    list_of_sets = []
    list_of_items = []
    for key, values in dict.items():
        key_plus_values = set()
        key_plus_values.add(key)
        key_plus_values.update(set(values))    
        subset_flag = False
        for items in list_of_sets:
            if key_plus_values.issubset(items):
                subset_flag = True
        if subset_flag == True:
            continue
        list_of_sets.append(key_plus_values)
    
    list_of_items = [
        j
        for i in list_of_sets
        for j in i
    ]

    return (list_of_sets, list_of_items)


def get_pcnt_authors_overlap(authors1, authors2):
    set1 = set(authors1.split(', '))
    set2 = set(authors2.split(', '))
    shared_set = set1 & set2
    combined_set = set1 | set2
    pcnt = len(shared_set) / len(combined_set)
    return pcnt


def get_pcnt_shared_accessions(stringlist1, stringlist2):
    list1 = stringlist1.split(', ')
    list2 = stringlist2.split(', ')
    pcnt_shared_accessions = len(set(list1) & set(list2)) / len(set(list1) | set(list2))
    return pcnt_shared_accessions


def get_pcnt_shared_stems(list1, list2):
    acc_num_stem_list1 = []
    acc_num_stem_list2 = []
    for acc_num in list1:
        acc_num_stem_list1.append(acc_num[:3])
    for acc_num in list2:
        acc_num_stem_list2.append(acc_num[:3])
    set_acc_num_stem1 = set(acc_num_stem_list1)
    set_acc_num_stem2 = set(acc_num_stem_list2)
    pcnt_shared_stems = len(set_acc_num_stem1 & set_acc_num_stem2) / len((set_acc_num_stem1 | set_acc_num_stem2))
    return pcnt_shared_stems


def calc_year_dif(list1, list2):
    years1 = list1.split(", ")
    years2 = list2.split(", ")
    years1 = [int(year) for year in years1 if year]
    years2 = [int(year) for year in years2 if year]
    max_dif = 0
    for year1 in years1:
        for year2 in years2:
            if abs(year1 - year2) > max_dif:
                max_dif = abs(year1 - year2)
    return max_dif


# The score is 1 if the author list is the same and 1 of the following criteria is met
# (1) One or both of the titles are 'Direct Submission" and the years are within 1. If the years
# are not present they are treated as being the same
# (2) The titles are the same (and not both 'Direct Submission')
def compare_authors_titles(row_i, row_j):
    authors_i = row_i['authors']
    authors_j = row_j['authors']
    title_i = row_i['title']
    title_j = row_j['title']
    year_i = str(row_i['year'])
    year_j = str(row_j['year'])
    pmid_i = str(row_i['pmid'])
    pmid_j = str(row_j['pmid'])
    accessions_i = row_i['accession']
    accessions_j = row_j['accession']

    pcnt_authors_overlap = get_pcnt_authors_overlap(authors_i, authors_j)
    title_distance = Levenshtein.distance(title_i, title_j)
    max_year_dif = calc_year_dif(year_i, year_j)
    pcnt_shared_accessions = get_pcnt_shared_accessions(accessions_i, accessions_j)
    pcnt_shared_stems = get_pcnt_shared_stems(accessions_i, accessions_j)

    if title_i != 'Direct Submission' and title_distance < 5:
        with open('SameTitles.txt', 'a') as file:
            file.write(f'Title_i:{title_i}\nTitle_j:{title_j}\nTitle_distance:{title_distance}\n')     
            file.write(f'Authors_i: {authors_i}\nAuthors_j: {authors_j}\n')
            file.write(f'Authors_overlap:{pcnt_authors_overlap}\n')
            file.write(f'Year_i:{year_i} Year_j:{year_j} Max_year_dif:{max_year_dif}\n')
            file.write(f'Accessions_i:{accessions_i}\nAccessions_j:{accessions_j}\n')
            file.write(f'Pcnt_shared_accessions:{pcnt_shared_accessions} Pcnt_shared_stems:{pcnt_shared_stems}\n\n')   
        return 1 
    elif (title_i == 'Direct Submission') | (title_j == 'Direct Submission') \
        and authors_i != 'NCBI' \
        and pcnt_authors_overlap >= 0.75 \
        and max_year_dif <= 1 \
        and pcnt_shared_stems >0.75:
        with open('SameAuthors.txt', 'a') as file:
            file.write(f'Index_i: {row_i.name} Index_j: {row_j.name}\n')
            file.write(f'Title_i:{title_i}\nTitle_j:{title_j}\nTitle_distance:{title_distance}\n')     
            file.write(f'Authors_i: {authors_i}\nAuthors_j: {authors_j}\n')
            file.write(f'Authors_overlap:{pcnt_authors_overlap}\n')
            file.write(f'Year_i:{year_i} Year_j:{year_j} Max_year_dif:{max_year_dif}\n')
            file.write(f'Accessions_i:{accessions_i}\nAccessions_j:{accessions_j}\n')
            file.write(f'Pcnt_shared_accessions:{pcnt_shared_accessions} Pcnt_shared_stems:{pcnt_shared_stems}\n\n')                  
        return 1
    else: 
        return 0


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
    print("Close lists:", close_lists)

    list_of_sets_w_shared_indexes, complete_list_of_shared_indexes = convert_dict_to_list_of_sets(close_lists)  
    print(f'''No with shared author_titles: {len(list_of_sets_w_shared_indexes)}: {list_of_sets_w_shared_indexes}''')
    print(f'''To be dropped: {len(complete_list_of_shared_indexes)}: {complete_list_of_shared_indexes}''')      

    list_of_new_rows = []
    for item in list_of_sets_w_shared_indexes:
        new_row = merge_refs_sharing_author_title(df, list(item))
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
        new_row = merge_refs_sharing_accessions(df, list(item))
        list_of_new_rows.append(new_row)

    df = df.drop(complete_list_of_shared_indexes)
    list_of_new_rows.insert(0, df)
    list_of_new_rows = pd.concat(list_of_new_rows, ignore_index=True)
    return list_of_new_rows


# Each list is a string in which items are separated by ', '
def combine_items_in_different_lists(lists):
    unique_dict = {}
    for string in lists:
        if isinstance(string, int):
            string = str(string)
        items = string.split(', ')
        for item in items:
            unique_dict[item] = 1
    list_unique_items = list(unique_dict.keys())
    list_unique_items = sorted(list_unique_items)
    unique_items = ', '.join(list_unique_items)
    if unique_items.startswith(","):
        unique_items = unique_items[1:]
    unique_items = unique_items.strip()
    return unique_items


def merge_refs_sharing_accessions(df, shared_indexes):
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


def merge_refs_sharing_author_title(df, shared_indexes):
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

    #new_row['authors'] = authors[:1]
    new_row['authors'] = new_authors
    #new_row['year'] = [';'.join(str(i) for i in year_list if i)]
    new_row['year'] = new_years
    new_row['title'] = new_titles
    new_row['pmid'] = new_pmids
    new_row['journal'] = new_journals
     #new_row['title'] = ['\n'.join(i for i in titles_list)]
    #new_row['pmid'] = ['; '.join(i for i in pmid_list if i)]
    #new_row['journal'] = ['; '.join(journal_list)]

    new_row['accession'] = new_accessions
    new_row_df = pd.DataFrame(new_row, index=[0])
    return new_row_df

